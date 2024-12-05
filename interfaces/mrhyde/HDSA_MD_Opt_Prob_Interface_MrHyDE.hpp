#ifndef HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP

template <class RealT>
class MD_Opt_Prob_Interface_MrHyDE : public HDSA::MD_Opt_Prob_Interface<RealT> {

private:

  HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode> > postproc_;
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > solver_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode> > params_;
  
public:
  MD_Opt_Prob_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > & solver, HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode> > & postproc, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode> > & params,const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator)
  {  
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<State_Vector_MrHyDE<RealT> > (solver,random_number_generator);
        
    postproc_ = postproc;
    solver_ = solver;
    params_ = params;

    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set=0; set<solver_->setnames.size(); set++) {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode> >(solver_->settings);
    }
  }

  virtual ~MD_Opt_Prob_Interface_MrHyDE()
  { }
  
  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const {      
    writedata_solopt(u_in);
    do_solop(true);
    gradient(z_out,z);
  }

  void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const {
    do_solop(false);
    HDSA::Ptr<HDSA::Vector<RealT> > g = z.clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_pert = z.clone();
    gradient(*g,z);
    z_pert->set(z);
    RealT h = 1.e-4;
    z_pert->axpy(h,z_in);
    gradient(z_out,*z_pert);
    z_out.axpy(-1.0,*g);
    z_out.scale(1/h);
  }

  void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const{
    const State_Vector_MrHyDE<RealT> &eu = dynamic_cast<const State_Vector_MrHyDE<RealT>&>(u);
    State_Vector_MrHyDE<RealT> &eu_grad = dynamic_cast<State_Vector_MrHyDE<RealT>&>(u_grad);
    for (int set = 0; set<eu.mrhyde_state_vec.size(); ++set) {
      for (int i = 0; i<solver_->numsteps[set]; i++) {
	RealT currenttime = solver_->initial_time + (double)i*solver_->deltat;
	postproc_->computeObjectiveGradState(set,eu.mrhyde_state_vec[set][i], currenttime,solver_->deltat,eu_grad.mrhyde_state_vec[set][i]);
      }
    }
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const {
    HDSA::Ptr<HDSA::Vector<RealT> > ugrad_nom = u_out.clone();
    Misfit_Gradient(*ugrad_nom,u,z);
    HDSA::Ptr<HDSA::Vector<RealT> > u_pert = u_out.clone();
    u_pert->set(u);
    RealT h = 1.0E-4;
    u_pert->axpy(h,u_in);
    Misfit_Gradient(u_out,*u_pert,z);
    u_out.axpy(-1.0,*ugrad_nom);
    u_out.scale(1.0/h);
  }

  void writedata_solopt(const HDSA::Vector<RealT> &u) const {
    const State_Vector_MrHyDE<RealT> &eu = dynamic_cast<const State_Vector_MrHyDE<RealT>&>(u);
    for (int set = 0; set<eu.mrhyde_state_vec.size(); ++set) {
      for (int i = 0; i<solver_->numsteps[set]; i++) { 
        RealT currenttime = solver_->initial_time + (double)i*solver_->deltat;
        postproc_->hdsa_solop_data[set]->store(eu.mrhyde_state_vec[set][i],currenttime,0);
      }
    }
  }

  void gradient(HDSA::Vector<RealT> & grad_z, const HDSA::Vector<RealT> & z) const {      
    bool new_z = checkNewParams(z);
    
    if (new_z) {
      MrHyDE_OptVector curr_z = params_->getCurrentVector();
      ROL::Ptr<ROL::Vector<RealT> > z_tmp = curr_z.clone();
      MrHyDE_OptVector ez_tmp = Teuchos::dyn_cast<MrHyDE_OptVector >(dynamic_cast<ROL::Vector<RealT> &>(*z_tmp));
      const Opt_Vector_MrHyDE<RealT> &ez = dynamic_cast<const Opt_Vector_MrHyDE<RealT>&>(z);
      ez_tmp.set(*ez.mrhyde_vec);
      
      params_->updateParams(ez_tmp);
      DFAD val = 0.0;
      solver_->forwardModel(val);
    }
    Opt_Vector_MrHyDE<RealT> e_grad_z = Teuchos::dyn_cast<Opt_Vector_MrHyDE<RealT> >(dynamic_cast<HDSA::Vector<RealT> &>(grad_z));
    e_grad_z.zeros();
    MrHyDE_OptVector ee_grad_z = Teuchos::dyn_cast<MrHyDE_OptVector >(dynamic_cast<ROL::Vector<RealT> &>(*e_grad_z.mrhyde_vec));
    
    solver_->adjointModel(ee_grad_z);
  }
  
  bool checkNewParams(const HDSA::Vector<RealT> &z) const {
    MrHyDE_OptVector curr_z = params_->getCurrentVector();
    ROL::Ptr<ROL::Vector<RealT> > diff = curr_z.clone();
    MrHyDE_OptVector ediff = Teuchos::dyn_cast<MrHyDE_OptVector >(const_cast<ROL::Vector<RealT> &>(*diff));
    Opt_Vector_MrHyDE<RealT> ez = Teuchos::dyn_cast<Opt_Vector_MrHyDE<RealT> >(const_cast<HDSA::Vector<RealT> &>(z));
    ediff.zero();
    ediff.set(curr_z);
    ediff.axpy(-1.0,*ez.mrhyde_vec);
    ScalarT dnorm = ediff.norm();
    ScalarT refnorm = curr_z.norm();
    dnorm = dnorm/refnorm;
    ScalarT reltol = 1.0e-12;
    bool new_z = false;
    if (dnorm > reltol) {
        new_z = true;
    }
    return new_z;
  }

  void do_solop(bool solop_flag) const
  {
    postproc_->hdsa_solop = solop_flag;
  } 
};
#endif
