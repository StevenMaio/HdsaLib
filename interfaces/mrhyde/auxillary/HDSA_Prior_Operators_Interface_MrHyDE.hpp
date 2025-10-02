#ifndef HDSA_PRIOR_OPERATORS_INTERFACE_MRHYDE_HPP
#define HDSA_PRIOR_OPERATORS_INTERFACE_MRHYDE_HPP

template <class RealT,
          class LO = Tpetra::Map<>::local_ordinal_type,
          class GO = Tpetra::Map<>::global_ordinal_type,
          class Node = Tpetra::Map<>::node_type>
class Prior_Operators_Interface_MrHyDE
{

public:
  Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> M;
  Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> S;

  Prior_Operators_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int>> &comm, Teuchos::RCP<Teuchos::ParameterList> &settings, std::vector<string> &blockNames)
  {
    Teuchos::RCP<Teuchos::ParameterList> settings_prior = HDSA::makePtr<Teuchos::ParameterList>(*settings);
    settings_prior->remove("Physics");
    settings_prior->sublist("Physics").set("modules", "ellipticPrior");
    settings_prior->sublist("Solver").set("solver", "steady-state");
    settings_prior->sublist("Solver").set("matrix free", true);
    settings_prior->remove("Analysis");
    settings_prior->sublist("Analysis").set("Analysis type", "forward");
    settings_prior->remove("Functions");
    settings_prior->sublist("Functions").set("ellipticPrior diffusion", "0.0");
    settings_prior->sublist("Functions").set("ellipticPrior reaction", "1.0");
    settings_prior->sublist("Functions").set("specific heat", "0.0");
    settings_prior->remove("Postprocess");
    settings_prior->sublist("Postprocess").set("write solution", false);
    settings_prior->sublist("Postprocess").set("create optimization movie", false);

    M = Instantiate_Prior_Operators(comm, settings_prior, blockNames);
    settings_prior->sublist("Functions").set("ellipticPrior diffusion", "1.0");
    settings_prior->sublist("Functions").set("ellipticPrior reaction", "0.0");
    S = Instantiate_Prior_Operators(comm, settings_prior, blockNames);
  }

  virtual ~Prior_Operators_Interface_MrHyDE()
  {
  }

  Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int>> comm, Teuchos::RCP<Teuchos::ParameterList> &settings, std::vector<string> &blockNames)
  {

    Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(settings, comm));

    Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp(new MrHyDE::PhysicsInterface(settings, comm,
                                                                                               mesh->getBlockNames(),
                                                                                               mesh->getSideNames(),
                                                                                               mesh->getDimension()));

    mesh->finalize(physics->getVarList(), physics->getVarTypes(), physics->getDerivedList());

    Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp(new MrHyDE::DiscretizationInterface(settings, comm,
                                                                                                          mesh, physics));

    Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> params = Teuchos::rcp(new MrHyDE::ParameterManager<SolverNode>(comm, settings, mesh, physics, disc));

    Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode>> assembler = Teuchos::rcp(new MrHyDE::AssemblyManager<SolverNode>(comm, settings, mesh, disc, physics, params));

    assembler->setMeshData();

    Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp(new MrHyDE::MultiscaleManager(comm, mesh, settings,
                                                                                                            assembler->groups,
#ifndef MrHyDE_NO_AD
                                                                                                            assembler->function_managers_AD));
#else
                                                                                                            assembler->function_managers));
#endif

    Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>>
        postproc = Teuchos::rcp(new MrHyDE::PostprocessManager<SolverNode>(comm, settings, mesh,
                                                                           disc, physics,
                                                                           multiscale_manager,
                                                                           assembler, params));

    Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solve = Teuchos::rcp(new MrHyDE::SolverManager<SolverNode>(comm, settings, mesh, disc, physics, assembler, params));

    solve->multiscale_manager = multiscale_manager;
    assembler->multiscale_manager = multiscale_manager;
    solve->postproc = postproc;

    mesh->allocateMeshDataStructures();
    assembler->allocateGroupStorage();
    solve->completeSetup();
    postproc->linalg = solve->linalg;
    solve->setupExplicitMass();
    assembler->finalizeFunctions();
    solve->finalizeMultiscale();

    Kokkos::fence();
    comm->barrier();

    std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT, LO, GO, Node>>> veci_out;
    veci_out.resize(1);
    veci_out[0] = solve->linalg->getNewVector(0);
    Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> A_over = solve->linalg->getNewOverlappedMatrix(0);

    assembler->assembleJacRes(0, 0, veci_out, veci_out, veci_out, veci_out, veci_out, veci_out, true, false, false,
                              veci_out[0], A_over, false, 0, false, false, veci_out[0]->getGlobalLength(), veci_out[0], veci_out[0], false, 0.0);
    solve->linalg->fillComplete(A_over);

    Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> A = solve->linalg->getNewMatrix(0);
    A->resumeFill();
    solve->linalg->exportMatrixFromOverlapped(0, A, A_over);
    solve->linalg->fillComplete(A);
    return A;
  }
};
#endif
