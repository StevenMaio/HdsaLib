#ifndef OEDLIB_DRIVER_MRHYDE_HPP
#define OEDLIB_DRIVER_MRHYDE_HPP

#include <iostream>

#include "Tpetra_Map_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "OED_MrHyDE_Model.hpp"
#include "OED_MrHyDE_Observation_Operator.hpp"
#include "OED_MrHyDE_Prior_Interface.hpp"
#include "OED_Gaussian_Error.hpp"

namespace OED::MrHyDE_Interface
{

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

        // TODO: may need to add other templates
        OED::Ptr<MrHyDE_Model_Interface<RealT>> model_;
        OED::Ptr<MrHyDE_Observation_Operator_Interface<RealT>> obervation_operator_;
        OED::Ptr<MrHyDE_Prior_Interface<RealT>> prior_;
        OED::Ptr<OED::Error_Model_Interface<RealT>> error_model_;

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

            // Exit if not OED settings are available
            if (!settings_->sublist("Analysis").isSublist("OED"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE could not find the OED sublist in the input file!  Abort!");
            }
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            // TODO: Load RNG settings (maybe for another time)

            this->Create_Model_Interface();
            this->Create_Observation_Operator_Interface();
            this->Create_Error_Model_Interface();
            this->Create_Prior_Interface();

            // TODO: do some OED -- and then do somehing with the result
        }

    private:
        inline void Create_Model_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Model_Interface Hello!" << std::endl;
        }

        inline void Create_Observation_Operator_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Observation_Operator_Interface Hello!" << std::endl;
        }

        inline void Create_Error_Model_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            // TODO: eventually do more here
            if (!oedSettings.isSublist("error model"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified error model for OED! Abort!");
            }
            Teuchos::ParameterList &error_parameters = oedSettings.sublist("error model");
            RealT noise_std = error_parameters.get<RealT>("std", 0.0);
            int dimension = error_parameters.get<int>("dimension", 0);

            this->error_model_ = OED::makePtr<OED::Gaussian_Error<RealT>>(dimension, noise_std);
            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Error_Model_Interface Hello!" << std::endl;
        }

        inline void Create_Prior_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            std::vector<std::string> blockNames = this->solver_->mesh->getBlockNames();
            this->prior_ = OED::makePtr<MrHyDE_Prior_Interface<RealT>>(this->comm_, *this->settings_, blockNames);
            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Prior_Interface Hello!" << std::endl;
        }
    };

}

#endif // OEDLIB_DRIVER_MRHYDE_HPP
