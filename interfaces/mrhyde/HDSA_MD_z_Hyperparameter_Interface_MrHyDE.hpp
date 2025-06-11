#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_MRHYDE_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_MRHYDE_HPP

template <class RealT>
class MD_z_Hyperparameter_Interface_MrHyDE : public HDSA::MD_z_Hyperparameter_Interface<RealT>
{

private:
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> solver_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;

public:
  std::vector<std::vector<RealT>> Spatial_Domain_Bounds(void) const override
  {
    // This is a placeholder. TODO: Generalize to read nodal data
    std::vector<std::vector<RealT>> vec; // vec.size() = spatial dimension, e.g. 1,2, or 3, [ vec[i][0],vec[i][1] ] is an interval bounding the ith spatial coordinate
    vec.resize(1);
    vec[0].resize(2);
    vec[0][0] = 0.0;
    vec[0][1] = 1.0;
    return vec;
  }

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const override
  {

    const HDSA_Tpetra_Vector<RealT> &ez = dynamic_cast<const HDSA_Tpetra_Vector<RealT> &>(z);
    HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> fvec = ez.getVector();
    ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<RealT>>(0);
    ROL::Ptr<MrHyDE_OptVector> z_rol = ROL::makePtr<MrHyDE_OptVector>(fvec, svec);

    params_->updateParams(*z_rol);
    ScalarT val = 0.0;
    solver_->forwardModel(val);

    if (solver_->isTransient)
    {
      Transient_Vector<RealT> &u_trans = dynamic_cast<Transient_Vector<RealT> &>(u);
      int n_t = solver_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
      for (int i = 0; i < n_t; i++)
      {
        HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_vec;
        solver_->postproc->soln[0]->extract(u_vec, i);
        HDSA_Tpetra_Vector<RealT> &eu_i = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(*u_trans[i]);
        eu_i.getVector()->update(1.0, *u_vec, 0.0);
      }
    }
    else
    {
      HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_vec;
      solver_->postproc->soln[0]->extract(u_vec, 0);
      HDSA_Tpetra_Vector<RealT> &eu = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(u);
      eu.getVector()->update(1.0, *u_vec, 0.0);
    }

  }

  MD_z_Hyperparameter_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(z_type, num_state_solves)
  {
    solver_ = solver;
    params_ = params;
  }

  MD_z_Hyperparameter_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, int seed, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(seed, z_type, num_state_solves)
  {
    solver_ = solver;
    params_ = params;
  }

  MD_z_Hyperparameter_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, z_type, num_state_solves)
  {
    solver_ = solver;
    params_ = params;
  }
};

#endif
