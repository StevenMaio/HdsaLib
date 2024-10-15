#ifndef HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_MRHYDE_HPP
#define HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_MRHYDE_HPP

#include <algorithm>
#include <cstdlib>
#include <random>
#include "userInterface.hpp"

  template <class RealT,
            class LO=Tpetra::Map<>::local_ordinal_type,
            class GO=Tpetra::Map<>::global_ordinal_type,
            class Node=Tpetra::Map<>::node_type >
  class MD_Elliptic_u_Prior_Interface_MrHyDE : public HDSA::MD_Elliptic_u_Prior_Interface<RealT> {
  
  private:
    std::vector<Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > > solve_;
    std::vector<Teuchos::RCP<MrHyDE::ParameterManager<SolverNode> > > params_;
  public:
    MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat): HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u)
    { 
    Instantiate_Prior_Operators(comm,settings,blockNames,solve_concat);
    }

    MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, int & seed, Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat): HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u,seed)
    { 
    Instantiate_Prior_Operators(comm,settings,blockNames,solve_concat);
    }

    MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator, Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<Teuchos::ParameterList> & settings, std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat)
      : HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u,random_number_generator)
    {
    Instantiate_Prior_Operators(comm,settings,blockNames,solve_concat);
    }

    virtual ~MD_Elliptic_u_Prior_Interface_MrHyDE()
    { }

  void Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int> > comm, Teuchos::RCP<Teuchos::ParameterList> & settings,std::vector<string> & blockNames, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve_concat) {

    params_.resize(blockNames.size());
    solve_.resize(blockNames.size());
    for (int i=0; i<blockNames.size(); ++i){
    std::string input_file_name = "input-HDSA-prior.yaml";
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
    
    Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode> > assembler = Teuchos::rcp( new MrHyDE::AssemblyManager<SolverNode>(comm, settings, mesh,
                                                                                                         disc, physics, params_[i]));
    
    assembler->setMeshData();
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the subgrid discretizations/models if using multiscale method
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp( new MrHyDE::MultiscaleManager(comm, mesh, settings,
                                                                                             assembler->groups,
                                                                                             #ifndef MrHyDE_NO_AD
                                                                                             assembler->function_managers_AD) );
                                                                                             #else
                                                                                             assembler->function_managers) );
                                                                                             #endif
    
    ///////////////////////////////////////////////////////////////////////////////
    // Create the postprocessing object
    ////////////////////////////////////////////////////////////////////////////////
    
    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode> >
    postproc = Teuchos::rcp( new MrHyDE::PostprocessManager<SolverNode>(comm, settings, mesh,
                                                                disc, physics, //assembler->function_managers_AD, 
                                                                multiscale_manager,
                                                                assembler, params_[i]) );
    
    ////////////////////////////////////////////////////////////////////////////////
    // Set up the solver and finalize some objects
    ////////////////////////////////////////////////////////////////////////////////
      
    solve_[i] = Teuchos::rcp( new MrHyDE::SolverManager<SolverNode>(comm, settings, mesh, disc, physics, assembler, params_[i]) );
    

    solve_[i]->multiscale_manager = multiscale_manager;
    assembler->multiscale_manager = multiscale_manager;
    solve_[i]->postproc = postproc;
    
    ////////////////////////////////////////////////////////////////////////////////
    // Allocate most of the per-element memory usage
    ////////////////////////////////////////////////////////////////////////////////
    
    mesh->allocateMeshDataStructures();
    assembler->allocateGroupStorage();

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
      assembler->purgeMemory();
      params_[i]->purgeMemory();
      physics->purgeMemory();
      if (debug_level > 0 && comm->getRank() == 0) {
        std::cout << "******** Finished driver memory purge ..." << std::endl;
      } 
           
    }
    
    solve_[i]->completeSetup();
    postproc->linalg = solve_[i]->linalg;
    
    if (settings->get<bool>("enable memory purge",false)) {
      disc->purgeMemory();
    }

    ////////////////////////////////////////////////////////////////////////////////
    // Finalize the function and multiscale managers
    ////////////////////////////////////////////////////////////////////////////////
    
    assembler->finalizeFunctions();

    solve_[i]->finalizeMultiscale();

    ////////////////////////////////////////////////////////////////////////////////
    // Perform the requested analysis (fwd solve, adj solve, dakota run, etc.)
    ////////////////////////////////////////////////////////////////////////////////
    
    // Make sure all processes are caught up at this point
    Kokkos::fence();
    comm->barrier();
    
    } // end for loop

    int num_sing_vals = 50;
    int oversampling = 1;
    int num_subspace_iters = 2;
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = HDSA::makePtr<HDSA::Vector_MrHyDE_Steady_State<RealT> >(solve_concat);

    HDSA::MD_Elliptic_u_Prior_Interface<RealT>::Compute_E_u_Inverse_GSVD(num_sing_vals,oversampling,num_subspace_iters,*u_vec);

  }

    void Apply_E_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {
     const HDSA::Vector_MrHyDE_Steady_State<RealT> &eu_in = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(u_in);  
     HDSA::Vector_MrHyDE_Steady_State<RealT> &eu_out = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(u_out);  
     for (int i=0; i<eu_in.mrhyde_steady_state_vec.size(); ++i){
       std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_in;
       std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > veci_out;
       veci_in.resize(1);
       veci_out.resize(1);
       veci_in[0]=eu_in.mrhyde_steady_state_vec[i];
       veci_out[0]=eu_out.mrhyde_steady_state_vec[i];
       solve_[i]->stateSolve(veci_out,veci_in); 
     }
    }
    void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {
      Apply_E_u_Inverse(u_out, u_in);
    }

    void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {

      const HDSA::Vector_MrHyDE_Steady_State<RealT> &eu_in = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(u_in);
      HDSA::Vector_MrHyDE_Steady_State<RealT> &eu_out = dynamic_cast<HDSA::Vector_MrHyDE_Steady_State<RealT>&>(u_out);
      for (int i=0; i<eu_in.mrhyde_steady_state_vec.size(); ++i){
        solve_[i]->explicitMass[0]->apply(*eu_in.mrhyde_steady_state_vec[0],*eu_out.mrhyde_steady_state_vec[0]);
      }
    }
  };

#endif

