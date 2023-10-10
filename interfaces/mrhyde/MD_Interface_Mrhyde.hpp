#ifndef HDSA_MD_INTERFACE_MRHYDE_HPP
#define HDSA_MD_INTERFACE_MRHYDE_HPP

template <class RealT>
class Model_Discrepancy_Interface_Mrhyde : public HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT> {

private:
  int m_; // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_; // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S_; // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M_; // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Minv_; // Mass matrix inverse
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > E_L_; // State elliptic operator
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma_; // Control weighting matrix
  Teuchos::RCP<HDSA::Objective_Mrhyde<RealT> > obj_;
  Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > postproc_;
  HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > solver_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode> > params_;
public:

  Model_Discrepancy_Interface_Mrhyde( HDSA::Ptr<MrHyDE::SolverManager<SolverNode> > & solver, Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> > & postproc, Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > & params, Teuchos::RCP<Teuchos::MpiComm<int> > &comm)
  {  
    obj_ = Teuchos::rcp( new HDSA::Objective_Mrhyde<RealT> (solver, postproc, params));
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<HDSA::Vector_Mrhyde_State<RealT> > (solver); 

    postproc_ = postproc;
    solver_ = solver;
    params_= params;
    postproc_->hdsa_solop_data.resize(solver_->setnames.size());
    for (int set=0; set<solver_->setnames.size(); set++) {
      postproc_->hdsa_solop_data[set] = HDSA::makePtr<MrHyDE::SolutionStorage<SolverNode> >();
    }

    Instantiate_Prior_Operators(comm);
    int num_sing_vals = 50;
    int oversampling = 1;
    int num_subspace_iters = 2;
    // HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<Std_Vector<RealT> >(m_);

    HDSA::Model_Discrepancy_Interface_Elliptic_Prior<RealT>::Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,*u_vec);
  }

  virtual ~Model_Discrepancy_Interface_Mrhyde()
  { }

  void Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int> > comm) {
    std::string input_file_name = "input-HDSA-prior.yaml";


    ////////////////////////////////////////////////////////////////////////////////
    // Import default and user-defined settings into a parameter list
    ////////////////////////////////////////////////////////////////////////////////
    
    //Teuchos::RCP<Teuchos::ParameterList> settings = HDSA::Mrhyde_UserInterface.Mrhyde_settings(input_file_name);

    HDSA::Ptr<HDSA::Mrhyde_UserInterface> mrhyde_user_interface = HDSA::makePtr<HDSA::Mrhyde_UserInterface>(input_file_name);

    //   Teuchos::RCP<Teuchos::ParameterList> settings = mrhyde_user_interface->Mrhyde_settings(input_file_name);    
    //    mrhyde_user_interface->settings_ = mrhyde_user_interface->Mrhyde_settings(input_file_name);    
    
    //    verbosity = settings->get<int>("verbosity",0);
    //  debug_level = settings->get<int>("debug level",0);
    // profile = settings->get<bool>("profile",false);
    
    ////////////////////////////////////////////////////////////////////////////////
    // Create the mesh
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(mrhyde_user_interface->settings_, comm) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the physics
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp( new MrHyDE::PhysicsInterface(mrhyde_user_interface->settings_, comm, mesh->stk_mesh) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Mesh only needs the variable names and types to finalize
    ////////////////////////////////////////////////////////////////////////////////
    
    mesh->finalize(physics);
    
    ////////////////////////////////////////////////////////////////////////////////
    // Define the discretization(s)
    ////////////////////////////////////////////////////////////////////////////////
        
    Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp( new MrHyDE::DiscretizationInterface(mrhyde_user_interface->settings_, comm,
                                                                                           mesh->stk_mesh,
                                                                                           physics) );
            
    ////////////////////////////////////////////////////////////////////////////////
    // Create the solver object
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > params = Teuchos::rcp( new MrHyDE::ParameterManager<SolverNode>(comm, mrhyde_user_interface->settings_,
                                                                                                        mesh->stk_mesh, physics, disc));
    
    Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode> > assembler = Teuchos::rcp( new MrHyDE::AssemblyManager<SolverNode>(comm, mrhyde_user_interface->settings_, mesh,
                                                                                                         disc, physics, params));
    
    mesh->setMeshData(assembler->groups,
                      assembler->boundary_groups);
    
    if (mrhyde_user_interface->settings_->get<bool>("enable memory purge",true)) {
      disc->purgeLIDs();
      if (!mrhyde_user_interface->settings_->sublist("Postprocess").get("write solution",false) && 
          !mrhyde_user_interface->settings_->sublist("Postprocess").get("create optimization movie",false)) {
        mesh->stk_mesh = Teuchos::null;
        mesh->mesh_factory = Teuchos::null;
        disc->mesh = Teuchos::null;
        params->mesh = Teuchos::null;
      }
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Set up the subgrid discretizations/models if using multiscale method
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp( new MrHyDE::MultiscaleManager(comm, mesh, mrhyde_user_interface->settings_,
                                                                                             assembler->groups,
                                                                                             assembler->function_managers_AD) );
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create the postprocessing object
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> >
      postproc = Teuchos::rcp( new MrHyDE::PostprocessManager<SolverNode>(comm, mrhyde_user_interface->settings_, mesh,
                                                                disc, physics, //assembler->function_managers_AD, 
                                                                multiscale_manager,
								  assembler, params) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the solver and finalize some objects
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solve = Teuchos::rcp( new MrHyDE::SolverManager<SolverNode>(comm, mrhyde_user_interface->settings_, mesh,
                                                                                                 disc, physics, assembler, params) );
    
    solve->multiscale_manager = multiscale_manager;
    assembler->multiscale_manager = multiscale_manager;
    solve->postproc = postproc;
    postproc->linalg = solve->linalg;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Purge Panzer memory before solving
    ////////////////////////////////////////////////////////////////////////////////

    if (mrhyde_user_interface->settings_->get<bool>("enable memory purge",true)) {

      disc->purgeMemory();
      mesh->purgeMemory();
      assembler->purgeMemory();
      params->purgeMemory();
      physics->purgeMemory();
           
    }
    
    assembler->allocateGroupStorage();

    solve->completeSetup();

    ////////////////////////////////////////////////////////////////////////////////
    // Finalize the function and multiscale managers
    ////////////////////////////////////////////////////////////////////////////////
    
    assembler->finalizeFunctions();

    solve->finalizeMultiscale();

    ////////////////////////////////////////////////////////////////////////////////
    // Perform the requested analysis (fwd solve, adj solve, dakota run, etc.)
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::AnalysisManager> analysis = Teuchos::rcp( new MrHyDE::AnalysisManager(comm, mrhyde_user_interface->settings_,
                                                                               solve, postproc, params) );
    
    // Make sure all processes are caught up at this point
    Kokkos::fence();
    comm->barrier();
    
    analysis->run();
    
    //    Kokkos::finalize();
  }

  void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    // state elliptic operator invers, need linear solve
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    // Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	b->Replace_Element(k,0,u_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	u_out_std.Replace_Element(k,(*x)(k,0));
    //   }
    u_out.set(u_in);
  }

  void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    // state elliptic operator invers, need linear solve, this function is needed if symmetric depending on BCs
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    // Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	b->Replace_Element(k,0,u_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*E_L_,*x,*b);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	u_out_std.Replace_Element(k,(*x)(k,0));
    //   }
    u_out.set(u_in);
  }

  void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    // expose FE interface to expose mass matrix
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    // Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	b->Replace_Element(k,0,u_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // M_->Multiply(*x,*b);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	u_out_std.Replace_Element(k,(*x)(k,0));
    //   }
    u_out.set(u_in);
  }

  void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
  {
    // execute linear with mass matrix
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT>& u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    // Std_Vector<RealT>& u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	b->Replace_Element(k,0,u_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // Minv_->Multiply(*x,*b);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	u_out_std.Replace_Element(k,(*x)(k,0));
    //   }
    u_out.set(u_in);
  }

  // Manipulate prior covariances using direct linear algebra since the dimension is small
  void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    // need logicals 
    // elliptic solve, need to determine how to hard code an elliptic solve instead of yaml
    // use case z on the mesh
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > b = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT>& z_in_std = dynamic_cast<const Std_Vector<RealT>&>(z_in);
    // Std_Vector<RealT>& z_out_std = dynamic_cast<Std_Vector<RealT>&>(z_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	b->Replace_Element(k,0,z_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*Gamma_,*x,*b); // Gamma is already multiplied by mass matrix - lines above this fucntion 
    // for(int k = 0; k < m_; k++)
    //   {
    // 	z_out_std.Replace_Element(k,(*x)(k,0));
    //   }
    z_out.set(z_in);
  }



  // Assume a constraint u = z^3 nodewise on the mesh defined by nodes in x_
  // Assume an objective (1/2)*(u-T)^t*M*(u-T) where T = (x_+1.0)^3 so that the optimal solution is u_opt=(x_+1.0)^3 and z_opt=x_+1.0
  // Assume a high-fidelity model u = z^3 + .2*z^2

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const 
  {
      //solutionStorage type for datagen_soln
      //create two postproc objective each with a different objective function
      //datagen = u_in;
    writedata_solopt(u_in);
    obj_->do_solop(true);
    RealT tol = 1.0E-7;
    obj_->gradient(z_out,z,tol);
  }

  // This implementation assumes that it is evaluated at the optimal z so that the adjoint=0, a more general implementation would include a term multiplied by the adjoint variable
  //bvbw  void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out,í ¼í·±í ¼í·° const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const
  void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
  {
    //    HDSA::Vector_Mrhyde<RealT> &ez_out = dynamic_cast<HDSA::Vector_Mrhyde<RealT>&>(z_out);
    //const HDSA::Vector_Mrhywhde<RealT> &ez_in = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(z_in);
    //const HDSA::Vector_Mrhyde<RealT> &ez = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(z);
    obj_->do_solop(false);
    RealT tol = 1E-8;
    obj_->hessVec(z_out, z_in, z, tol );
    
  }

  void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const 
  {
    // J_u - need full space objective function gradient
    // use postproc.computeobjective to get gradient and access the gradient through Sacado access  (full space objective)
    const HDSA::Vector_Mrhyde_State<RealT> &eu = dynamic_cast<const HDSA::Vector_Mrhyde_State<RealT>&>(u);
    HDSA::Vector_Mrhyde_State<RealT> &eu_grad = dynamic_cast<HDSA::Vector_Mrhyde_State<RealT>&>(u_grad);
    for (int set=0; set<eu.mrhyde_state_vec.size(); ++set) {
      for (int i=0; i<solver_->numsteps[set]; i++) {
	RealT currenttime = solver_->initial_time + (double)i*solver_->deltat;
	postproc_->computeObjectiveGradState(set,eu.mrhyde_state_vec[set][i], currenttime,solver_->deltat,eu_grad.mrhyde_state_vec[set][i]);
      }
    }

    //    u_grad.set(u);
    // need to implement our won misfit Gradient dependent of the objective function - this is the RHS adjoint solve
    // we traced it down but ran into access issues - Kokkos interfaces, multiscale interfaces....vector was not exposed

 // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
 //    const Std_Vector<RealT> u_std = dynamic_cast<const Std_Vector<RealT>&>(u);
 //    Std_Vector<RealT> u_grad_std = dynamic_cast<Std_Vector<RealT>&>(u_grad);
 //    for(int k = 0; k < m_; k++)
 //      {
 // 	v->Replace_Element(k,0,u_std(k)-std::pow((*x_)(k,0)+1.0,3.0));
 //      }
 //    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > grad = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
 //    M_->Multiply(*grad,*v);
 //    for(int k = 0; k < m_; k++)
 //      {
 // 	u_grad_std.Replace_Element(k,(*grad)(k,0));
 //      }

  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const
  {
    // J_uu - need full space objective function hessian
    // finite difference gradient
    // evaluate gradient u+U_in - u step size of one since quadratic, will need smaller step if non-quadratic
    // long term implement Hessians
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > v = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // const Std_Vector<RealT> u_in_std = dynamic_cast<const Std_Vector<RealT>&>(u_in);
    // Std_Vector<RealT> u_out_std = dynamic_cast<Std_Vector<RealT>&>(u_out);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	v->Replace_Element(k,0,u_in_std(k));
    //   }
    // HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Hv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1);
    // M_->Multiply(*Hv,*v);
    // for(int k = 0; k < m_; k++)
    //   {
    // 	u_out_std.Replace_Element(k,(*Hv)(k,0));
    //   }

    HDSA::Ptr<HDSA::Vector<RealT> > ugrad_nom = u_out.clone();
    Misfit_Gradient(*ugrad_nom,u,z);
    HDSA::Ptr<HDSA::Vector<RealT> > u_pert = u_out.clone();
    u_pert->set(u);
    RealT h = 1.0E-4;
    u_pert->axpy(h,u_in);
    Misfit_Gradient(u_out,*u_pert,z);
    u_out.axpy(-1.0,*ugrad_nom);
    u_out.scale(1.0/h);
    //    u_out.set(u_in);
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u( ) const 
  {
    // logic to read Exodus file that contains optimal state solution from MrHyDE
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<HDSA::Vector_Mrhyde_State<RealT> > (solver_); 
    u_opt->Test_Vector();
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z( ) const 
  {
    MrHyDE_OptVector vectmp = params_->getCurrentVector();
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt = HDSA::makePtr<HDSA::Vector_Mrhyde<RealT> >(vectmp);
    // logic to read Exodus file that contains optimal control solution from MrHyDE
    z_opt->Test_Vector();
    return z_opt;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data( ) const 
  {
    MrHyDE_OptVector vectmp = params_->getCurrentVector();
    HDSA::Ptr<HDSA::Vector<RealT> > z = HDSA::makePtr<HDSA::Vector_Mrhyde<RealT> >(vectmp);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*z);
    HDSA::Ptr<HDSA::Vector<RealT> > z0 = (*Z)[0];
    HDSA::Ptr<HDSA::Vector<RealT> > z1 = (*Z)[1];
    z0->randomize_standard_normal();
    z1->randomize_standard_normal();



    // HDSA::Ptr<HDSA::Vector<RealT> > z0 = (*Z)[0];
    // HDSA::Ptr<HDSA::Vector<RealT> > z1 = (*Z)[1];
    // Std_Vector<RealT> z0_std = dynamic_cast<Std_Vector<RealT>&>(*z0);
    // Std_Vector<RealT> z1_std = dynamic_cast<Std_Vector<RealT>&>(*z1);

    // for(int k = 0; k < m_; k++)
    //   {
    // 	z0_std.Replace_Element(k,(*x_)(k,0)+1.0);
    // 	z1_std.Replace_Element(k,(*x_)(k,0) + std::pow((*x_)(k,0),2.0));
    //   }

    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data( ) const 
  {
    // logic to choose which data target
    // use yaml specification, need to pass yaml reader to contructor
    // if else blocks

    HDSA::Ptr<HDSA::Vector<RealT> > u = HDSA::makePtr<HDSA::Vector_Mrhyde_State<RealT> > (solver_);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*u);
    HDSA::Ptr<HDSA::Vector<RealT> > y0 = (*Y)[0];
    HDSA::Ptr<HDSA::Vector<RealT> > y1 = (*Y)[1];
    y0->randomize_standard_normal();
    y1->randomize_standard_normal();

    // HDSA::Ptr<Std_Vector<RealT> > y = HDSA::makePtr<Std_Vector<RealT> >(m_);
    // HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*y);

    // HDSA::Ptr<HDSA::Vector<RealT> > y0 = (*Y)[0];
    // HDSA::Ptr<HDSA::Vector<RealT> > y1 = (*Y)[1];
    // Std_Vector<RealT> y0_std = dynamic_cast<Std_Vector<RealT>&>(*y0);
    // Std_Vector<RealT> y1_std = dynamic_cast<Std_Vector<RealT>&>(*y1);

    // for(int k = 0; k < m_; k++)
    //   {
    // 	y0_std.Replace_Element(k,0.2*std::pow((*x_)(k,0)+1.0,2.0));
    // 	y1_std.Replace_Element(k,0.2*std::pow((*x_)(k,0) + std::pow((*x_)(k,0),2.0),2.0));
    //   }

    return Y;
  }

  void writedata_solopt(const HDSA::Vector<RealT> &u) const {
    const HDSA::Vector_Mrhyde_State<RealT> &eu = dynamic_cast<const HDSA::Vector_Mrhyde_State<RealT>&>(u);
    for (int set=0; set<eu.mrhyde_state_vec.size(); ++set) {
      for (int i=0; i<solver_->numsteps[set]; i++) {
	RealT currenttime = solver_->initial_time + (double)i*solver_->deltat;
	postproc_->hdsa_solop_data[set]->store(eu.mrhyde_state_vec[set][i],currenttime,0);
      }
    }
  }
  // HDSA::Ptr<HDSA::Vector<RealT> > Load_Matlab_z_Update( ) const 
  // {
  //   HDSA::Ptr<Std_Vector<RealT> > z_opt = HDSA::makePtr<Std_Vector<RealT> >(m_);

  //   RealT val = 0.0;
  //   // read in data
  //    std::ifstream in("z_update_matlab_solution.txt");           
  //   // read the elements in the file into a vector  
  //   // test file open   
  //   if (in) {   
  //     for(int i = 0; i < m_; i++)
  // 	{
  // 	  in >> val;
  // 	  z_opt->Replace_Element(i,val);
  // 	}   
  //   }
  //   else
  //     {
  // 	std::cout << "Error loading the data from z_update_matlab_solution.txt" << std::endl;
  //     }  
  //   return z_opt;
  // }

};

#endif


