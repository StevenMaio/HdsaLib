#ifndef OEDLIB_DRIVER_MRHYDE_HPP
#define OEDLIB_DRIVER_MRHYDE_HPP

#include <iostream>
#include <format>

#include "Tpetra_Map_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "OED_Bayesian_Inversion_Interface.hpp"
#include "OED_Lazy_Greedy.hpp"
#include "OED_Linear_OED_D_Opt.hpp"
#include "OED_Active_Sensors.hpp"

#include "OED_MrHyDE_Model_Interface.hpp"
#include "OED_MrHyDE_Observation_Operator.hpp"
#include "OED_MrHyDE_Prior_Interface.hpp"
#include "OED_Gaussian_Error.hpp"

#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Component_Observation_Operator.hpp"

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
        OED::Ptr<OED::Trilinos_Adapter::Random_Number_Generator<ScalarT>> rng_;
        int oed_verbosity_{0};

        // TODO: may need to add other templates
        OED::Ptr<MrHyDE_Model_Interface<RealT>> model_;
        OED::Ptr<OED::Observation_Operator_Interface<RealT>> observation_operator_;
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
            postproc_->write_solution = false;
            postproc_->write_optimization_solution = false;

            // TODO: insert OED here
            std::cout << "Hello OED!" << std::endl;

            // Exit if not OED settings are available
            if (!settings_->sublist("Analysis").isSublist("OED"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE could not find the OED sublist in the input file!  Abort!");
            }
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            // TODO: Load RNG settings (maybe for another time)
            OED::Ptr<OED::Trilinos_Adapter::Comm<int>> oed_comm = OED::makePtr<OED::Trilinos_Adapter::Comm<int>>(comm_);
            this->rng_ = OED::makePtr<OED::Trilinos_Adapter::Random_Number_Generator<ScalarT>>(oed_comm);

            this->Create_Model_Interface();
            this->Create_Observation_Operator_Interface();
            this->Create_Error_Model_Interface();
            this->Create_Prior_Interface();

            std::string dir = "/home/smaio/Documents/Projects/AIVIS/HdsaLib/OedLib/examples/mrhyde/thermal_1D/test";
            // TODO: Test mass matrix and creating empty vectors
            {

                // should wrap this in one of our wrappers -- then we can at least confine things there
                auto v_in = this->model_->Get_Empty_Parameter_Vector();
                auto v_out = this->model_->Get_Empty_Parameter_Vector();

                int dim = this->model_->Param_Dimension();
                std::cout << "dimension: " << dim << std::endl;

                for (int i = 0; i < dim; i++)
                {
                    v_in->Zeros();
                    v_out->Zeros();
                    v_in->Set_Entry(i, 1.0);

                    this->prior_->Mass_Matrix_Apply(*v_out, *v_in);
                    std::string output_file = std::format("{}/Me_{}.txt", dir, i);

                    v_out->Write_To_File(output_file);
                }
            }

            // TODO: test the forward solve
            {
                std::string output_file = std::format("{}/u.txt", dir);
                auto m = this->model_->Get_Empty_Parameter_Vector();
                auto u = this->model_->Get_Empty_State_Vector();

                m->Set_Scalar(2.0);
                this->model_->State_Solve(*u, *m);
                u->Write_To_File(output_file);
            }

            // TODO: do some OED -- and then do somehing with the result
            auto inversion_problem = OED::makePtr<OED::Bayesian_Inversion_Interface<RealT>>(this->model_, this->observation_operator_, this->prior_, this->error_model_);

            auto oed_problem = OED::makePtr<OED::Linear_OED_D_Opt<RealT>>(inversion_problem);
            int budget = 5;
            int data_dim = this->error_model_->Data_Dimension();
            OED::Active_Sensors design = OED::Lazy_Greedy_Solve(*oed_problem, data_dim, budget);
            std::cout << "EIG: " << oed_problem->Evaluate(design) << std::endl;
            design.Print_Sensors();
        }

    private:
        inline void Create_Model_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");
            this->model_ = OED::makePtr<MrHyDE_Model_Interface<RealT>>(this->comm_, this->settings_, this->solver_, this->postproc_, this->params_, this->rng_);
            
            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Model_Interface Hello!" << std::endl;
        }

        inline void Create_Observation_Operator_Interface()
        {
            Teuchos::ParameterList &oedSettings = settings_->sublist("Analysis").sublist("OED");

            if (!oedSettings.isSublist("observation operator"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified observation operator for OED! Abort!");
            }
            Teuchos::ParameterList obs_settings = oedSettings.sublist("observation operator");
            std::string &type = obs_settings.get<std::string>("type", "ERROR");
            if (type == "components")
            {
                std::string indices_file = obs_settings.get<std::string>("indices file", "ERROR");
                // TODO: do error checking
                this->observation_operator_ = OED::makePtr<OED::Component_Observation_Operator<RealT>>(this->model_->State_Dimension(), indices_file);
            }
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
            int dimension = this->observation_operator_->Data_Dimension();

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
