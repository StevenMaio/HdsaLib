#ifndef HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OPT_PROB_INTERFACE_MRHYDE_HPP

  template <class RealT>
  class MD_Opt_Prob_Interface_MrHyDE : public HDSA::MD_Opt_Prob_Interface<RealT> {

  private:

    Teuchos::RCP<HDSA::Objective_MrHyDE<RealT> > obj_;
    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > postproc_;
    HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > solver_;
  public:
    MD_Opt_Prob_Interface_MrHyDE(HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > & solver, Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > & postproc, Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > & params)
  {  
    obj_ = Teuchos::rcp( new HDSA::Objective_MrHyDE<RealT> (solver, postproc, params));
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<HDSA::Vector_MrHyDE_State<RealT> > (solver);
        
    postproc_ = postproc;
    solver_ = solver;

    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set=0; set<solver_->setnames.size(); set++) {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode> >(solver_->settings);
    }

    }

    virtual ~MD_Opt_Prob_Interface_MrHyDE()
    { }

    void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const {      
      writedata_solopt(u_in);
      obj_->do_solop(true);
      RealT tol = 1.0E-7;
      obj_->gradient(z_out,z,tol);
    }

    void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const {
      //    HDSA::Vector_Mrhyde<RealT> &ez_out = dynamic_cast<HDSA::Vector_Mrhyde<RealT>&>(z_out);
      //const HDSA::Vector_Mrhywhde<RealT> &ez_in = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(z_in);
      //const HDSA::Vector_Mrhyde<RealT> &ez = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(z);
      obj_->do_solop(false);
      RealT tol = 1E-8;
      obj_->hessVec(z_out, z_in, z, tol );
    }

    void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const{
      // J_u - need full space objective function gradient
      // use postproc.computeobjective to get gradient and access the gradient through Sacado access  (full space objective)
      const HDSA::Vector_MrHyDE_State<RealT> &eu = dynamic_cast<const HDSA::Vector_MrHyDE_State<RealT>&>(u);
      HDSA::Vector_MrHyDE_State<RealT> &eu_grad = dynamic_cast<HDSA::Vector_MrHyDE_State<RealT>&>(u_grad);
      for (int set=0; set<eu.mrhyde_state_vec.size(); ++set) {
        for (int i=0; i<solver_->numsteps[set]; i++) {
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
    const HDSA::Vector_MrHyDE_State<RealT> &eu = dynamic_cast<const HDSA::Vector_MrHyDE_State<RealT>&>(u);
    for (int set=0; set<eu.mrhyde_state_vec.size(); ++set) {
      for (int i=0; i<solver_->numsteps[set]; i++) { 
        RealT currenttime = solver_->initial_time + (double)i*solver_->deltat;
        postproc_->hdsa_solop_data[set]->store(eu.mrhyde_state_vec[set][i],currenttime,0);
      }
    }
  } 
                                                                                                                                                                                                     
  };

#endif
