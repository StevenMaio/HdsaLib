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

  // { 
  //   Teuchos::RCP<Teuchos::ParameterList> settings_prior = HDSA::makePtr<Teuchos::ParameterList>();

  //   settings_prior->setParameters(settings->sublist("Mesh"));
  //   settings_prior->setParameters(settings->sublist("Discretization"));
  //   //    settings_prior->set("Discretization",settings->sublist("Discretization"));
  //   Teuchos::RCP<Teuchos::ParameterList> settings_physics = HDSA::makePtr<Teuchos::ParameterList>();
  //   settings_prior->sub("Physics").set("modules","ellipticPrior");
  //   settings_prior->sublist("Solver").set("solver","steady-state");
  //   settings_prior->sublist("Solver").set("matrix free",true);
  //   settings_prior->sublist("Analysis").set("Analysis type","forward");
    
  //   settings_prior->sublist("Functions").set("ellipticPrior diffusion","0.0");
  //   settings_prior->sublist("Functions").set("ellipticPrior reaction","1.0");
  //   M = Instantiate_Prior_Operators(comm,settings_prior,blockNames);

  //   settings_prior->sublist("Functions").set("ellipticPrior diffusion","1.0");
  //   settings_prior->sublist("Functions").set("ellipticPrior reaction","0.0");
  //   S = Instantiate_Prior_Operators(comm,settings_prior,blockNames);
  // }
    
  virtual ~HDSA_Prior_FE_Op_MrHyDE_Interface()
  { }

  std::vector<matrix_RCP> Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int> > comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames) {

    std::vector<matrix_RCP> A;
    A.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
      
      int debug_level = settings->get<int>("debug level",0);
    
      Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(settings, comm) );
  
      Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp( new MrHyDE::PhysicsInterface(settings, comm, 
												  mesh->getBlockNames(),
												  mesh->getSideNames(),
												  mesh->getDimension()) );
    
      mesh->finalize(physics->getVarList(), physics->getVarTypes(), physics->getDerivedList());
    
      Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp( new MrHyDE::DiscretizationInterface(settings, comm,
													     mesh, physics) );
     
      Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > params = Teuchos::rcp( new MrHyDE::ParameterManager<SolverNode>(comm, settings, mesh, physics, disc));
    
      Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode> > assembler = Teuchos::rcp( new MrHyDE::AssemblyManager<SolverNode>(comm, settings, mesh, disc, physics, params));
    
      assembler->setMeshData();
    
      Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp( new MrHyDE::MultiscaleManager(comm, mesh, settings,
													       assembler->groups,
													       #ifndef MrHyDE_NO_AD
													       assembler->function_managers_AD) );
      #else
	 assembler->function_managers) );
      #endif
    
  Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> >
	postproc = Teuchos::rcp( new MrHyDE::PostprocessManager<SolverNode>(comm, settings, mesh,
									    disc, physics, //assembler->function_managers_AD, 
									    multiscale_manager,
									    assembler, params) );
    
      Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solve = Teuchos::rcp( new MrHyDE::SolverManager<SolverNode>(comm, settings, mesh, disc, physics, assembler, params) );
      
      
      solve->multiscale_manager = multiscale_manager;
      assembler->multiscale_manager = multiscale_manager;
      solve->postproc = postproc;
      
      mesh->allocateMeshDataStructures();
      assembler->allocateGroupStorage();
      
      if (settings->get<bool>("enable memory purge",false)) {
	if (debug_level > 0 && comm->getRank() == 0) {
	  std::cout << "******** Starting driver memory purge ..." << std::endl;
	}
	if (!settings->sublist("Postprocess").get("write solution",false) && 
	    !settings->sublist("Postprocess").get("create optimization movie",false)) {
	  mesh->purgeMesh();
	  disc->mesh = Teuchos::null;
	  params->mesh = Teuchos::null;
	}
	disc->purgeOrientations();
	disc->purgeLIDs();
	mesh->purgeMemory();
	assembler->purgeMemory();
	params->purgeMemory();
	physics->purgeMemory();
	if (debug_level > 0 && comm->getRank() == 0) {
	  std::cout << "******** Finished driver memory purge ..." << std::endl;
	} 
      }
      
      solve->completeSetup();
      postproc->linalg = solve->linalg;
      solve->setupExplicitMass();
      
      if (settings->get<bool>("enable memory purge",false)) {
	disc->purgeMemory();
      }

      assembler->finalizeFunctions();
      
      solve->finalizeMultiscale();

      Kokkos::fence();
      comm->barrier();

      A[i] = solve->linalg->getNewMatrix(i);
      std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_out;
      veci_out.resize(1);
      veci_out[0]=solve->linalg->getNewVector(i);
      matrix_RCP A_over = solve->linalg->getNewOverlappedMatrix(i);
	    
      assembler->assembleJacRes(i,0,veci_out,veci_out,veci_out,veci_out,veci_out,veci_out,true,false,false,
				 veci_out[0],A_over,false,0,false,false,veci_out[0]->getGlobalLength(),veci_out[0],veci_out[0],false,0.0);
      solve->linalg->fillComplete(A_over);
      A[i]->resumeFill();
      solve->linalg->exportMatrixFromOverlapped(i, A[i], A_over);
      solve->linalg->fillComplete(A[i]);
  
    }
    return A;
  } 
};
#endif

