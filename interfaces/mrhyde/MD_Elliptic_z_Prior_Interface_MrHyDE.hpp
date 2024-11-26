#ifndef HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP
#define HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP


  template <class RealT,
            class LO=Tpetra::Map<>::local_ordinal_type,
            class GO=Tpetra::Map<>::global_ordinal_type,
            class Node=Tpetra::Map<>::node_type >
  class MD_Elliptic_z_Prior_Interface_MrHyDE : public HDSA::MD_Elliptic_z_Prior_Interface<RealT> {

  typedef Tpetra::CrsMatrix<ScalarT,LO,GO,Node>   LA_CrsMatrix;
  typedef Teuchos::RCP<LA_CrsMatrix>              matrix_RCP;

  private:

   std::vector<Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > > solve_;
   std::vector<Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > > params_;
   Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode> > assembler_;
    std::vector<matrix_RCP> Ez_; 
  public:
    MD_Elliptic_z_Prior_Interface_MrHyDE(RealT & alpha_z, Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat): HDSA::MD_Elliptic_z_Prior_Interface<RealT>(alpha_z)
    { 
    Instantiate_Prior_Operators(comm,settings,blockNames,solve_concat);
    }


    
    virtual ~MD_Elliptic_z_Prior_Interface_MrHyDE()
    { }

  void Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int> > comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat) {

    params_.resize(blockNames.size());
    solve_.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
    std::string input_file_name = "input-HDSA-prior-z.yaml";
    // Todo: loop over blockname and instantiate solve and params per block
    //       in apply_u_inverse: loop over blocknames, execute state solve per block

    ////////////////////////////////////////////////////////////////////////////////
    // Import default and user-defined settings into a parameter list
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::userInterface> UI = Teuchos::rcp(new MrHyDE::userInterface() );
    Teuchos::RCP<Teuchos::ParameterList> settings = UI->UserInterface(input_file_name);

    //Teuchos::RCP<Teuchos::ParameterList> settings = MrHyDE::UserInterface(input_file_name);
    
    int debug_level = settings->get<int>("debug level",0);
    
    ////////////////////////////////////////////////////////////////////////////////
    // Create the mesh
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(settings, comm) );
  
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the physics
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp( new MrHyDE::PhysicsInterface(settings, comm, 
                                                                                mesh->getBlockNames(),
                                                                                mesh->getSideNames(),
                                                                                mesh->getDimension()) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Mesh only needs the variable names and types to finalize
    ////////////////////////////////////////////////////////////////////////////////
    
    mesh->finalize(physics->getVarList(), physics->getVarTypes(), physics->getDerivedList());
    
    ////////////////////////////////////////////////////////////////////////////////
    // Define the discretization(s)
    ////////////////////////////////////////////////////////////////////////////////
        
    Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp( new MrHyDE::DiscretizationInterface(settings, comm,
                                                                                           mesh, physics) );
            
    ////////////////////////////////////////////////////////////////////////////////
    // Create the solver object
    ////////////////////////////////////////////////////////////////////////////////
    
    params_[i] = Teuchos::rcp( new MrHyDE::ParameterManager<SolverNode>(comm, settings, mesh, physics, disc));
    
    assembler_ = Teuchos::rcp( new MrHyDE::AssemblyManager<SolverNode>(comm, settings, mesh, disc, physics, params_[i]));
    
    assembler_->setMeshData();
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the subgrid discretizations/models if using multiscale method
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp( new MrHyDE::MultiscaleManager(comm, mesh, settings,
                                                                                             assembler_->groups,
                                                                                             #ifndef MrHyDE_NO_AD
                                                                                             assembler_->function_managers_AD) );
                                                                                             #else
                                                                                             assembler_->function_managers) );
                                                                                             #endif
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create the postprocessing object
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> >
    postproc = Teuchos::rcp( new MrHyDE::PostprocessManager<SolverNode>(comm, settings, mesh,
                                                                disc, physics, //assembler->function_managers_AD, 
                                                                multiscale_manager,
                                                                assembler_, params_[i]) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the solver and finalize some objects
    ////////////////////////////////////////////////////////////////////////////////
      
    solve_[i] = Teuchos::rcp( new MrHyDE::SolverManager<SolverNode>(comm, settings, mesh, disc, physics, assembler_, params_[i]) );
    

    solve_[i]->multiscale_manager = multiscale_manager;
    assembler_->multiscale_manager = multiscale_manager;
    solve_[i]->postproc = postproc;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Allocate most of the per-element memory usage
    ////////////////////////////////////////////////////////////////////////////////
    
    mesh->allocateMeshDataStructures();
    assembler_->allocateGroupStorage();

    ////////////////////////////////////////////////////////////////////////////////
    // Purge Panzer memory before solving
    ////////////////////////////////////////////////////////////////////////////////
      
    if (settings->get<bool>("enable memory purge",false)) {
      if (debug_level > 0 && comm->getRank() == 0) {
        std::cout << "******** Starting driver memory purge ..." << std::endl;
      }
      if (!settings->sublist("Postprocess").get("write solution",false) && 
          !settings->sublist("Postprocess").get("create optimization movie",false)) {
        mesh->purgeMesh();
        disc->mesh = Teuchos::null;
        params_[i]->mesh = Teuchos::null;
      }
      disc->purgeOrientations();
      disc->purgeLIDs();
      mesh->purgeMemory();
      assembler_->purgeMemory();
      params_[i]->purgeMemory();
      physics->purgeMemory();
      if (debug_level > 0 && comm->getRank() == 0) {
        std::cout << "******** Finished driver memory purge ..." << std::endl;
      } 
           
    }
    
    solve_[i]->completeSetup();
    postproc->linalg = solve_[i]->linalg;
    solve_[i]->setupExplicitMass();
    
    if (settings->get<bool>("enable memory purge",false)) {
      disc->purgeMemory();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Finalize the function and multiscale managers
    ////////////////////////////////////////////////////////////////////////////////
    
    assembler_->finalizeFunctions();

    solve_[i]->finalizeMultiscale();

    ////////////////////////////////////////////////////////////////////////////////
    // Perform the requested analysis (fwd solve, adj solve, dakota run, etc.)
    ////////////////////////////////////////////////////////////////////////////////
    
    // Make sure all processes are caught up at this point
    Kokkos::fence();
    comm->barrier();
    
    } // end for loop
    Ez_.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
      Ez_[i] = solve_[i]->linalg->getNewMatrix(i);
      std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_out;
      veci_out.resize(1);
      veci_out[0]=solve_[i]->linalg->getNewVector(i);
      matrix_RCP Ez_over = solve_[i]->linalg->getNewOverlappedMatrix(i);
	    
      assembler_->assembleJacRes(i,0,veci_out,veci_out,veci_out,veci_out,veci_out,veci_out,true,false,false,
				 veci_out[0],Ez_over,false,0,false,false,veci_out[0]->getGlobalLength(),veci_out[0],veci_out[0],false,0.0);
      solve_[i]->linalg->fillComplete(Ez_over);
      Ez_[i]->resumeFill();
      solve_[i]->linalg->exportMatrixFromOverlapped(i, Ez_[i], Ez_over);
      solve_[i]->linalg->fillComplete(Ez_[i]);
    }
  }

    void Apply_E_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
     const HDSA::Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const HDSA::Vector_MrHyDE<RealT>&>(z_in);  
     HDSA::Vector_MrHyDE<RealT> &ez_out = dynamic_cast<HDSA::Vector_MrHyDE<RealT>&>(z_out);

     const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
     MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);

     std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_in;
     std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_out;
     veci_in.resize(1);
     veci_out.resize(1);
     veci_in[0]=eez_in.getField()[0]->getVector();
     veci_out[0]=eez_out.getField()[0]->getVector();
     solve_[0]->stateSolve(veci_out,veci_in);
       
    }
    void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
      Apply_E_z_Inverse(z_out,z_in);
    } 

    void Apply_M_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {

      const HDSA::Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const HDSA::Vector_MrHyDE<RealT>&>(z_in);  
      HDSA::Vector_MrHyDE<RealT> &ez_out = dynamic_cast<HDSA::Vector_MrHyDE<RealT>&>(z_out);

      const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
      MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);

      //const HDSA::Vector_MrHyDE_Steady_State<RealT> &ez_in = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(z_in);
      //HDSA::Vector_MrHyDE_Steady_State<RealT> &ez_out = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(z_out);
      // for (int i=0; i<ez_in.mrhyde_steady_state_vec.size(); ++i){
      Teuchos::RCP<Tpetra::MultiVector<RealT,LO,GO,Node>> eez_tmp = eez_out.getField()[0]->getVector();
      if (solve_[0]->linalg->getHaveOverlapped()) {
	Teuchos::RCP<Tpetra::MultiVector<RealT,LO,GO,Node>> z_in_over=solve_[0]->linalg->getNewOverlappedVector(0);
	Teuchos::RCP<Tpetra::MultiVector<RealT,LO,GO,Node>> z_out_over=solve_[0]->linalg->getNewOverlappedVector(0);
	solve_[0]->linalg->importVectorToOverlapped(0,z_in_over,eez_in.getField()[0]->getVector());
	solve_[0]->assembler->applyMassMatrixFree(0, z_in_over, z_out_over);
	
	// solve_[0]->linalg->exportVectorFromOverlapped(0,eez_out.getField()[0]->getVector(),z_out_over);
	solve_[0]->linalg->exportVectorFromOverlapped(0,eez_tmp,z_out_over);
      } else {
	solve_[0]->assembler->applyMassMatrixFree(0,eez_in.getField()[0]->getVector(),eez_tmp);
      }
	  //}
    }

    // Compute samples from a mean zero Gaussian with covariance W_z^{-1}                                                                                                                                                              
    virtual void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> & samples) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "MD_Elliptic_z_Prior_Interface::Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms" << std::endl);
    }   
    
    virtual void Apply_E_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
    {
      const HDSA::Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const HDSA::Vector_MrHyDE<RealT>&>(z_in);  
      HDSA::Vector_MrHyDE<RealT> &ez_out = dynamic_cast<HDSA::Vector_MrHyDE<RealT>&>(z_out);

      const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
      MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);

      //const HDSA::Vector_MrHyDE_Steady_State<RealT> &ez_in = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(z_in);
      //HDSA::Vector_MrHyDE_Steady_State<RealT> &ez_out = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(z_out);
      //      for (int i=0; i<ez_in.mrhyde_steady_state_vec.size(); ++i){
      Ez_[0]->apply(*eez_in.getField()[0]->getVector(),*eez_out.getField()[0]->getVector());
	//}
    }

    virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      Apply_E_z(z_out,z_in);
    }

    virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "MD_Elliptic_z_Prior_Interface::Apply_M_z_Inverse must be implemented to use the Hessian GEVP" << std::endl);
         // identity
    }

  };

#endif

