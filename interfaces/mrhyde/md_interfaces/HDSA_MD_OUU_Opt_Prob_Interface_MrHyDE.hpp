/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_OUU_OPT_PROB_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OUU_OPT_PROB_INTERFACE_MRHYDE_HPP

#include "HDSA_MD_OUU_Opt_Prob_Interface.hpp"

template <class RealT>
class MD_OUU_Opt_Prob_Interface_MrHyDE : public HDSA::MD_OUU_Opt_Prob_Interface<RealT>
{

private:
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> solver_;
  HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> postproc_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  HDSA::Ptr<ROL::SampleGenerator<RealT>> sampler_;

  HDSA::Ptr<Solver_Interface_MrHyDE<RealT>> solver_interface_;
  std::vector<std::vector<std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>>>> state_soln_;
  std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> grad_nom_;

public:
  MD_OUU_Opt_Prob_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> &postproc, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, HDSA::Ptr<ROL::SampleGenerator<RealT>> &sampler, std::vector<RealT> &ens_weights) : HDSA::MD_OUU_Opt_Prob_Interface<RealT>(ens_weights), solver_(solver), postproc_(postproc), params_(params), sampler_(sampler)
  {
    solver_interface_ = HDSA::makePtr<Solver_Interface_MrHyDE<RealT>>(solver_, params_);

    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set = 0; set < solver_->setnames.size(); set++)
    {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode>>(solver_->settings);
    }

    int ens_size = ens_weights.size();
    Initialize_State_Solution(*data_interface->Get_z_opt(), ens_size);
    grad_nom_.resize(ens_size);
    for (int s = 0; s < ens_size; s++)
    {
      std::vector<RealT> pt_s = sampler_->getMyPoint(s);
      params_->updateParams(pt_s, "stochastic");
      grad_nom_[s] = data_interface->Get_z_opt()->Clone();
      RS_Gradient(*grad_nom_[s], *data_interface->Get_z_opt(), s, false);
    }
  }

  virtual ~MD_OUU_Opt_Prob_Interface_MrHyDE()
  {
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &z, int s) const
  {
    std::vector<RealT> pt_s = sampler_->getMyPoint(s);
    params_->updateParams(pt_s, "stochastic");
    Write_Data_Solution_Operator(u_in);
    Do_Solution_Operator(true);
    RS_Gradient(z_out, z, s, false);
    z_out.Scale(-1.0);
  }

  void Apply_RS_Hessian_Per_Sample(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, int s) const
  {
    std::vector<RealT> pt_s = sampler_->getMyPoint(s);
    params_->updateParams(pt_s, "stochastic");
    Do_Solution_Operator(false);
    HDSA::Ptr<HDSA::Vector<RealT>> z_pert = z.Clone();
    z_pert->Set(z);
    RealT h = 1.e-4;
    z_pert->Scaled_Plus(h, z_in);
    RS_Gradient(z_out, *z_pert, s, true);
    z_out.Scaled_Plus(-1.0, *grad_nom_[s]);
    z_out.Scale(1.0 / h);
  }

  void Misfit_Gradient_Per_Sample(HDSA::Vector<RealT> &u_grad, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, int s) const
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

  void Apply_Misfit_Hessian_Per_Sample(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, int s) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> ugrad_nom = u_out.Clone();
    Misfit_Gradient_Per_Sample(*ugrad_nom, u, z, s);

    HDSA::Ptr<HDSA::Vector<RealT>> u_pert = u_out.Clone();
    u_pert->Set(u);
    RealT h = 1.0e-4;
    u_pert->Scaled_Plus(h, u_in);
    Misfit_Gradient_Per_Sample(u_out, *u_pert, z, s);

    u_out.Scaled_Plus(-1.0, *ugrad_nom);
    u_out.Scale(1.0 / h);
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Functions that aid in the implementation of the base class pure virtual function
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void Write_Data_Solution_Operator(const HDSA::Vector<RealT> &u) const
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

  void Do_Solution_Operator(bool solop_flag) const
  {
    postproc_->hdsa_solop = solop_flag;
  }

  void Initialize_State_Solution(const HDSA::Vector<RealT> &z, int ens_size) 
  {
    Instantiate_Solution_Storage();
    
    HDSA::Ptr<MrHyDE_OptVector> z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(z);
    MrHyDE_OptVector curr_z = params_->getCurrentVector();
    ROL::Ptr<ROL::Vector<RealT>> z_tmp = curr_z.clone();
    MrHyDE_OptVector ez_tmp = Teuchos::dyn_cast<MrHyDE_OptVector>(dynamic_cast<ROL::Vector<RealT> &>(*z_tmp));
    ez_tmp.set(*z_rol);
    params_->updateParams(ez_tmp);

    for (int s = 0; s < ens_size; s++)
    {
      std::vector<RealT> pt_s = sampler_->getMyPoint(s);
      params_->updateParams(pt_s, "stochastic");

      ScalarT val = 0.0;
      solver_->forwardModel(val);

      Cache_State(s);
    }
  }

  void RS_Gradient(HDSA::Vector<RealT> &grad_z, const HDSA::Vector<RealT> &z, int s, bool force_fwd_solve = false) const
  {
    if (force_fwd_solve)
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
    else
    {
      Set_State(s);
    }

    HDSA::Ptr<MrHyDE_OptVector> grad_z_rol = solver_interface_->Map_HDSA_Vector_to_MrHyDE_OptVector(grad_z);
    grad_z_rol->zero();

    solver_->adjointModel(*grad_z_rol);
  }

  void Instantiate_Solution_Storage(void)
  {
    int num_samples = sampler_->numMySamples();
    state_soln_.resize(num_samples);
    for (int sample_id = 0; sample_id < num_samples; sample_id++)
    {
      state_soln_[sample_id].resize(postproc_->soln.size());
      for (int set = 0; set < state_soln_[sample_id].size(); set++)
      {
        int numFwdSteps = postproc_->soln[set]->getTotalTimes(0);
        state_soln_[sample_id][set].resize(numFwdSteps);
        for (int time_index = 0; time_index < numFwdSteps; time_index++)
        {
          state_soln_[sample_id][set][time_index] = solver_->linalg->getNewOverlappedVector(set);
        }
      }
    }
  }

  void Cache_State(int sample_id) const
  {
    for (int set = 0; set < postproc_->soln.size(); set++)
    {
      std::vector<std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>>> vecs = postproc_->soln[set]->extractAllData();
      int numFwdSteps = vecs[0].size();
      for (int time_index = 0; time_index < numFwdSteps; time_index++)
      {
        state_soln_[sample_id][set][time_index]->update(1.0, *vecs[0][time_index], 0.0);
      }
    }
  }

  void Set_State(int sample_id) const
  {
    for (int set = 0; set < postproc_->soln.size(); set++)
    {
      postproc_->soln[set]->reset();
      int numFwdSteps = state_soln_[sample_id][set].size();
      for (int time_index = 0; time_index < numFwdSteps; time_index++)
      {
        RealT current_time = solver_->initial_time + (double)time_index * solver_->deltat;
        postproc_->soln[set]->store(state_soln_[sample_id][set][time_index], current_time, 0);
      }
    }
  }
};

#endif
