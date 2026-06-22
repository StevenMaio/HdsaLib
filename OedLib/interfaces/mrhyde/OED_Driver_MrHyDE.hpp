#ifndef OEDLIB_DRIVER_MRHYDE_HPP
#define OEDLIB_DRIVER_MRHYDE_HPP

#include <iostream>

#include "Tpetra_Map_decl.hpp"
#include "Teuchos_RCPDecl.hpp"

namespace OED
{
    namespace MrHyDE_Interface {

        template <class RealT,
                class LO = Tpetra::Map<>::local_ordinal_type,
                class GO = Tpetra::Map<>::global_ordinal_type,
                class Node = Tpetra::Map<>::node_type>
        class Driver_MrHyDE
        {
        private:
            Teuchos::RCP<MpiComm> comm_;
            Teuchos::RCP<Teuchos::ParameterList> settings_;
            Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solver_;
            Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> postproc_;
            Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> params_;
            int oed_verbosity_{0};

        public:
            Driver_MrHyDE(Teuchos::RCP<MpiComm> &comm, Teuchos::RCP<Teuchos::ParameterList> &settings, Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solver,
                        Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> &postproc, Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> &params)
                : comm_(comm), settings_(settings), solver_(solver), postproc_(postproc), params_(params)
            {
            }

            void OED_Solve()
            {
                // TODO: insert OED here
                std::cout << "Hello OED!" << std::endl;

                // TODO: determine general configurations

                this->Create_Model_Interface();
                this->Create_Observation_Operator_Interface();
                this->Create_Error_Model_Interface();
                this->Create_Prior_Interface();

                // TODO: do some OED -- and then do somehing with the result
            }

        private:
            inline void Create_Model_Interface()
            {
            }

            inline void Create_Observation_Operator_Interface()
            {
            }

            inline void Create_Error_Model_Interface()
            {
            }

            inline void Create_Prior_Interface()
            {
            }
        };

    }

}

#endif // OEDLIB_DRIVER_MRHYDE_HPP
