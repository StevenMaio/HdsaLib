#ifndef HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP

#include "HDSA_MD_Opt_Prob_Interface.hpp"
#include "HDSA_Tpetra_Vector.hpp"
#include "HDSA_Std_Vector.hpp"
#include "HDSA_Solver_Interface_MrHyDE.hpp"

template <class RealT>
class MD_Opt_Prob_Interface_MrHyDE : public HDSA::MD_Opt_Prob_Interface<RealT>
{

private:
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> solver_;
  HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> postproc_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  
  HDSA::Ptr<Solver_Interface_MrHyDE<RealT>> solver_interface_;
  HDSA::Ptr<HDSA::Vector<RealT>> grad_nom_;

public:
  MD_Opt_Prob_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> &postproc, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface)
  {
    postproc_ = postproc;
    solver_ = solver;
    params_ = params;

    solver_interface_ = HDSA::makePtr<Solver_Interface_MrHyDE<RealT>>(solver_, params_);

    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set = 0; set < solver_->setnames.size(); set++)
    {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode>>(solver_->settings);
    }

    grad_nom_ = data_interface->get_z_opt()->Clone();
    gradient(*grad_nom_, *data_interface->get_z_opt());
  }

  virtual ~MD_Opt_Prob_Interface_MrHyDE()
  {
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Implementation of base class pure virtual functions:
  // Apply_Solution_Operator_z_Jacobian_Transpose
  // Apply_RS_Hessian
  // Misfit_Gradient
  // Apply_Misfit_Hessian
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &z) const
  {
    writedata_solopt(u_in);
    do_solop(true);
    gradient(z_out, z);
    z_out.Scale(-1.0);
  }

  void Apply_RS_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z) const
  {
    do_solop(false);
    HDSA::Ptr<HDSA::Vector<RealT>> z_pert = z.Clone();
    z_pert->Set(z);
    RealT h = 1.e-4;
    z_pert->axpy(h, z_in);
    gradient(z_out, *z_pert);
    z_out.axpy(-1.0, *grad_nom_);
    z_out.Scale(1.0 / h);
  }

  void Misfit_Gradient(HDSA::Vector<RealT> &u_grad, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {
    if (solver_->isTransient)
    {
      const HDSA::Transient_Vector<RealT> &eu = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(u);
      HDSA::Transient_Vector<RealT> &eu_grad = dynamic_cast<HDSA::Transient_Vector<RealT> &>(u_grad);
      int n_t = solver_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;

      if (postproc_->objectives[0].type == "integrated control")
      { // only works for one objective term
        eu_grad[0]->Zeros();
        for (int i = 0; i < n_t - 1; i++)
        { // exludes initial condition
          const HDSA::Tpetra_Vector<RealT> &eu_i = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(*eu[i + 1]);
          HDSA::Tpetra_Vector<RealT> &eu_grad_i = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(*eu_grad[i + 1]);
          RealT currenttime = solver_->initial_time + (double)i * solver_->deltat;

          // the gradient should be a non-overlapping vector, but the state should be an overlapping vector
          HDSA::Ptr<Tpetra::MultiVector<RealT>> eu_grad_i_tpetra = eu_grad_i.getVector();
          HDSA::Ptr<Tpetra::MultiVector<RealT>> ui_over = solver_->linalg->getNewOverlappedVector(0);
          solver_->linalg->importVectorToOverlapped(0, ui_over, eu_i.getVector());

          postproc_->setTimeIndex(i);
          solver_->assembler->updateStage(0, currenttime, solver_->deltat);
          postproc_->computeObjectiveGradState(0, ui_over, currenttime, solver_->deltat, eu_grad_i_tpetra);
          if (i == 0)
          {
            eu_grad_i.Scale(solver_->deltat);
          }
        }
        u_grad.Scale(-1.0);
      }
    }
    else
    {
      const HDSA::Tpetra_Vector<RealT> &eu = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(u);
      HDSA::Tpetra_Vector<RealT> &eu_grad = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(u_grad);

      // the gradient should be a non-overlapping vector, but the state should be an overlapping vector
      HDSA::Ptr<Tpetra::MultiVector<RealT>> grad = eu_grad.getVector();
      HDSA::Ptr<Tpetra::MultiVector<RealT>> u_over = solver_->linalg->getNewOverlappedVector(0);
      solver_->linalg->importVectorToOverlapped(0, u_over, eu.getVector());
      postproc_->computeObjectiveGradState(0, u_over, 0.0, solver_->deltat, grad);

      u_grad.Scale(-1.0);
    }
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> ugrad_nom = u_out.Clone();
    Misfit_Gradient(*ugrad_nom, u, z);

    HDSA::Ptr<HDSA::Vector<RealT>> u_pert = u_out.Clone();
    u_pert->Set(u);
    RealT h = 1.0e-4;
    u_pert->axpy(h, u_in);
    Misfit_Gradient(u_out, *u_pert, z);

    u_out.axpy(-1.0, *ugrad_nom);
    u_out.Scale(1.0 / h);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Functions that aid in the implementation of the base class pure virtual function
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void writedata_solopt(const HDSA::Vector<RealT> &u) const
  {
    if (solver_->isTransient)
    {
      const HDSA::Transient_Vector<RealT> &eu = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(u);
      int n_t = solver_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
      for (int i = 0; i < n_t; i++)
      {
        const HDSA::Tpetra_Vector<RealT> &eu_i = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(*eu[i]);
        RealT currenttime = solver_->initial_time + (double)i * solver_->deltat;
        postproc_->hdsa_solop_data[0]->store(eu_i.getVector(), currenttime, 0);
      }
    }
    else
    {
      const HDSA::Tpetra_Vector<RealT> &eu = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(u);
      postproc_->hdsa_solop_data[0]->store(eu.getVector(), 0.0, 0);
    }
  }

  void do_solop(bool solop_flag) const
  {
    postproc_->hdsa_solop = solop_flag;
  }

  void gradient(HDSA::Vector<RealT> &grad_z, const HDSA::Vector<RealT> &z) const
  {
    bool new_z = checkNewParams(z);
    if (new_z)
    {
      HDSA::Ptr<MrHyDE_OptVector> z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(z);
      MrHyDE_OptVector curr_z = params_->getCurrentVector();
      ROL::Ptr<ROL::Vector<RealT>> z_tmp = curr_z.clone();
      MrHyDE_OptVector ez_tmp = Teuchos::dyn_cast<MrHyDE_OptVector>(dynamic_cast<ROL::Vector<RealT> &>(*z_tmp));
      ez_tmp.set(*z_rol);

      params_->updateParams(ez_tmp);
      ScalarT val = 0.0;
      solver_->forwardModel(val);
    }

    HDSA::Ptr<MrHyDE_OptVector> grad_z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(grad_z);
    grad_z_rol->zero();

    solver_->adjointModel(*grad_z_rol);
  }

  bool checkNewParams(const HDSA::Vector<RealT> &z) const
  {
    HDSA::Ptr<MrHyDE_OptVector> z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(z);
    MrHyDE_OptVector curr_z = params_->getCurrentVector();

    ROL::Ptr<ROL::Vector<RealT>> diff = curr_z.clone();
    diff->zero();
    diff->set(curr_z);
    diff->axpy(-1.0, *z_rol);

    ScalarT dnorm = diff->norm();
    ScalarT refnorm = curr_z.norm();
    dnorm = dnorm / refnorm;
    ScalarT reltol = 1.0e-12;

    bool new_z = false;
    if (dnorm > reltol)
    {
      new_z = true;
    }
    return new_z;
  }

};
#endif
