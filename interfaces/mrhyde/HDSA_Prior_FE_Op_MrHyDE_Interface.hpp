#ifndef HDSA_PRIOR_FE_OP_MRHYDE_INTERFACE_HPP
#define HDSA_PRIOR_FE_OP_MRHYDE_INTERFACE_HPP
#include "userInterface.hpp"

  template <class RealT,
            class LO=Tpetra::Map<>::local_ordinal_type,
            class GO=Tpetra::Map<>::global_ordinal_type,
            class Node=Tpetra::Map<>::node_type >
  class HDSA_Prior_FE_Op_MrHyDE_Interface {

  typedef Tpetra::CrsMatrix<RealT,LO,GO,Node>   LA_CrsMatrix;
  typedef Teuchos::RCP<LA_CrsMatrix>              matrix_RCP;

  private:

   std::vector<Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > > solve_;
   std::vector<Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > > params_;
   Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode> > assembler_;

  public:
    std::vector<matrix_RCP> M;
    std::vector<matrix_RCP> S; 
    HDSA_Prior_FE_Op_MrHyDE_Interface(Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames)
    { 
    Teuchos::RCP<MrHyDE::userInterface> UI = Teuchos::rcp(new MrHyDE::userInterface() );
    std::string input_file_name = "input-HDSA-prior-z.yaml";
    Teuchos::RCP<Teuchos::ParameterList> settings_prior = UI->UserInterface(input_file_name);
    settings_prior->sublist("Functions").set("ellipticPrior diffusion","0.0");
    M = Instantiate_Prior_Operators(comm,settings_prior,blockNames);
    settings_prior->sublist("Functions").set("ellipticPrior diffusion","1.0");
    settings_prior->sublist("Functions").set("ellipticPrior reaction","0.0");
    S = Instantiate_Prior_Operators(comm,settings_prior,blockNames);
    }
    
    virtual ~HDSA_Prior_FE_Op_MrHyDE_Interface()
    { }

    std::vector<matrix_RCP> Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int> > comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames) {

    params_.resize(blockNames.size());
    solve_.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
    std::string input_file_name = "input-HDSA-prior-z.yaml";
    // Todo: loop over blockname and instantiate solve and params per block
    //       in apply_u_inverse: loop over blocknames, execute state solve per block

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
    std::vector<matrix_RCP> A;
    A.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
      A[i] = solve_[i]->linalg->getNewMatrix(i);
      std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_out;
      veci_out.resize(1);
      veci_out[0]=solve_[i]->linalg->getNewVector(i);
      matrix_RCP A_over = solve_[i]->linalg->getNewOverlappedMatrix(i);
	    
      assembler_->assembleJacRes(i,0,veci_out,veci_out,veci_out,veci_out,veci_out,veci_out,true,false,false,
				 veci_out[0],A_over,false,0,false,false,veci_out[0]->getGlobalLength(),veci_out[0],veci_out[0],false,0.0);
      solve_[i]->linalg->fillComplete(A_over);
      A[i]->resumeFill();
      solve_[i]->linalg->exportMatrixFromOverlapped(i, A[i], A_over);
      solve_[i]->linalg->fillComplete(A[i]);
    }
    return A;
  } 

  };

#endif

