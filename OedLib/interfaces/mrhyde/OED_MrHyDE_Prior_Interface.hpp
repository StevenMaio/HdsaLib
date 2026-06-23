#pragma once

#include "OED_Prior_Interface.hpp"
#include "OED_Vector.hpp"
#include "OED_Trilinos_Sparse_Matrix.hpp"
#include "OED_Sparse_Matrix_Solver.hpp"

#include "Tpetra_Map_decl.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT,
              class LO = Tpetra::Map<>::local_ordinal_type,
              class GO = Tpetra::Map<>::global_ordinal_type,
              class Node = Tpetra::Map<>::node_type>
    class MrHyDE_Prior_Interface : public OED::Prior_Interface<RealT>
    {
    private:
        // TODO: do I need to have these around?
        // Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> M_;
        // Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> S_;

        ScalarT gamma_{0};
        ScalarT delta_{0};

        // TODO: need to create solvers
        OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>> M_;
        OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>> S_;
        OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>> L_;

        OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix_Solver<RealT>> M_solver_;
        OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix_Solver<RealT>> L_solver_;

        Teuchos::RCP<Teuchos::MpiComm<int>> comm_;
        Teuchos::ParameterList &settings_;
        std::vector<string> &blockNames_;

    public:
        MrHyDE_Prior_Interface(Teuchos::RCP<Teuchos::MpiComm<int>> &comm,
            Teuchos::ParameterList &settings,
            std::vector<string> &blockNames)
            : comm_(comm), settings_(settings), blockNames_(blockNames)
        {
            // load OED settings
            Teuchos::ParameterList &oedSettings = settings_.sublist("Analysis").sublist("OED");

            if (!oedSettings.isSublist("prior"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified prior model for OED! Abort!");
            }
            // TODO: do error checking
            this->gamma_ = oedSettings.get<RealT>("gamma", 1e-3);
            this->delta_ = oedSettings.get<RealT>("delta", 1.0);

            // TODO: consider a mean as well

            // Create bilinar forms
            Teuchos::RCP<Teuchos::ParameterList> priorSettings = OED::makePtr<Teuchos::ParameterList>(settings);
            priorSettings->remove("Physics");
            priorSettings->sublist("Physics").set("modules", "ellipticPrior");
            priorSettings->sublist("Solver").set("solver", "steady-state");
            priorSettings->sublist("Solver").set("matrix free", true);
            priorSettings->remove("Analysis");
            priorSettings->sublist("Analysis").set("Analysis type", "forward");
            priorSettings->remove("Functions");
            priorSettings->sublist("Functions").set("ellipticPrior diffusion", "0.0");
            priorSettings->sublist("Functions").set("ellipticPrior reaction", "1.0");
            priorSettings->sublist("Functions").set("specific heat", "0.0");
            priorSettings->remove("Postprocess");
            priorSettings->sublist("Postprocess").set("write solution", false);
            priorSettings->sublist("Postprocess").set("create optimization movie", false);

            Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> M;
            Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> S;

            // TODO: should maybe convert these to other sparse matrices
            M = Instantiate_Prior_Operators(comm, priorSettings, blockNames);
            priorSettings->sublist("Functions").set("ellipticPrior diffusion", "1.0");
            priorSettings->sublist("Functions").set("ellipticPrior reaction", "0.0");
            S = Instantiate_Prior_Operators(comm, priorSettings, blockNames);

            this->M_ = OED::makePtr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>>(M, true);
            this->S_ = OED::makePtr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>>(S, true);

            this->L_ = OED::makePtr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>>(S, true);
            this->L_->Scale(this->gamma_);
            this->L_->Scaled_Plus(this->delta_, *this->M_);

            // TODO: create solvers for M and L
            this->M_solver_ = OED::makePtr<OED::Trilinos_Adapter::Sparse_Matrix_Solver<RealT>>(this->M_);
            this->L_solver_ = OED::makePtr<OED::Trilinos_Adapter::Sparse_Matrix_Solver<RealT>>(this->L_);
        }

        void Prior_Precision_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        void Prior_Covariance_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {
        }

        void Get_Prior_Mean(OED::Vector<RealT> &m_out) override
        {
            // TODO: support user provided means later on
            m_out.Zeros();
        }

        void Prior_Covariance_Factor_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        int Param_Dimension() override
        {
            // TODO: actually implement this
            return 0;
        }

        // TODO: move this at some point to a potentially new class
        void Mass_Matrix_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {
            this->M_->Apply(m_out, m_in);
        }

        void Mass_Matrix_Inverse_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {
            this->M_solver_->Apply_A_Inverse(m_out, m_in);
        }

        Ptr<OED::Vector<RealT>> Sample_Vector() override
        {
            // TODO: do something like a multi-vector I suppose
            return OED::nullPtr;
        }

    private:
        Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> Instantiate_Prior_Operators(Teuchos::RCP<Teuchos::MpiComm<int>> comm,
            Teuchos::RCP<Teuchos::ParameterList> &Settings,
            std::vector<string> &blockNames)
        {

            Teuchos::RCP<MrHyDE::MeshInterface> mesh = Teuchos::rcp(new MrHyDE::MeshInterface(Settings, comm));

            Teuchos::RCP<MrHyDE::PhysicsInterface> physics = Teuchos::rcp(new MrHyDE::PhysicsInterface(Settings, comm,
                                                                                                    mesh->getBlockNames(),
                                                                                                    mesh->getSideNames(),
                                                                                                    mesh->getDimension()));

            mesh->finalize(physics->getVarList(), physics->getVarTypes(), physics->getDerivedList());

            Teuchos::RCP<MrHyDE::DiscretizationInterface> disc = Teuchos::rcp(new MrHyDE::DiscretizationInterface(Settings, comm,
                                                                                                                mesh, physics));

            Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> params = Teuchos::rcp(new MrHyDE::ParameterManager<SolverNode>(comm, Settings, mesh, physics, disc));

            Teuchos::RCP<MrHyDE::AssemblyManager<SolverNode>> assembler = Teuchos::rcp(new MrHyDE::AssemblyManager<SolverNode>(comm, Settings, mesh, disc, physics, params));

            assembler->setMeshData();

            Teuchos::RCP<MrHyDE::MultiscaleManager> multiscale_manager = Teuchos::rcp(new MrHyDE::MultiscaleManager(comm, mesh, Settings,
                                                                                                                    assembler->groups,
#ifndef MrHyDE_NO_AD
                                                                                                                    assembler->function_managers_AD));
#else
                                                                                                                    assembler->function_managers));
#endif

            Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>>
                postproc = Teuchos::rcp(new MrHyDE::PostprocessManager<SolverNode>(comm, Settings, mesh,
                                                                                disc, physics,
                                                                                multiscale_manager,
                                                                                assembler, params));

            Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solve = Teuchos::rcp(new MrHyDE::SolverManager<SolverNode>(comm, Settings, mesh, disc, physics, assembler, params));

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

            int set = 0;
            int stage = 0;

            OED::Ptr<Tpetra::MultiVector<RealT, LO, GO, Node>> current_res, current_res_over;
            current_res = solve->linalg->getNewVector(set);
            current_res_over = solve->linalg->getNewOverlappedVector(set);

            Teuchos::RCP<Tpetra::CrsMatrix<ScalarT, LO, GO, Node>> J, J_over;
            J = solve->linalg->getNewMatrix(set);
            J_over = solve->linalg->getNewOverlappedMatrix(set);
            solve->linalg->fillComplete(J_over);
            J_over->resumeFill();
            J_over->setAllToScalar(0.0);

            auto paramvec = params->getDiscretizedParamsOver();
            auto paramdot = params->getDiscretizedParamsDotOver();

            std::vector<OED::Ptr<Tpetra::MultiVector<RealT, LO, GO, Node>>> zero_soln;
            assembler->assembleJacRes(set, stage, zero_soln, zero_soln, zero_soln, zero_soln, zero_soln, zero_soln, true, false, false, false, 0,
                                    current_res_over, J_over, false, 0.0, false, false,
                                    params->num_active_params, paramvec, paramdot, false, 0.0);
            solve->linalg->fillComplete(J_over);
            J->resumeFill();
            solve->linalg->exportMatrixFromOverlapped(set, J, J_over);
            solve->linalg->fillComplete(J);

            return J;
        }

    };
}