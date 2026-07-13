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
#include "OED_Transient_Observation_Operator.hpp"

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
            std::cout << "Hello OED!" << std::endl;

            // Exit if not OED settings are available
            if (!settings_->sublist("Analysis").isSublist("OED"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE could not find the OED sublist in the input file!  Abort!");
            }
            Teuchos::ParameterList &oed_settings = settings_->sublist("Analysis").sublist("OED");

            // TODO: Load custom RNG settings
            OED::Ptr<OED::Trilinos_Adapter::Comm<int>> oed_comm = OED::makePtr<OED::Trilinos_Adapter::Comm<int>>(comm_);
            this->rng_ = OED::makePtr<OED::Trilinos_Adapter::Random_Number_Generator<ScalarT>>(oed_comm);

            this->Create_Model_Interface();
            this->Create_Observation_Operator_Interface();
            this->Create_Error_Model_Interface();
            this->Create_Prior_Interface();

            auto inversion_problem = OED::makePtr<OED::Bayesian_Inversion_Interface<RealT>>(this->model_, this->observation_operator_, this->prior_, this->error_model_);

            // TODO: look at design criterion settings
            if (!oed_settings.isSublist("design criterion"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified design criterion for OED! Abort!");
            }
            Teuchos::ParameterList &obj_settings = oed_settings.sublist("design criterion");
            std::string &type = obj_settings.get<std::string>("type", "ERROR");
            if (type == "D-Optimality")
            {
                std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::OED_Solve criterion=D-Optimality" << std::endl;
                // TODO: handle other types of constraints
                auto oed_problem = OED::makePtr<OED::Linear_OED_D_Opt<RealT>>(inversion_problem);
                int budget = obj_settings.get<int>("budget", 0);;
                int data_dim = this->error_model_->Data_Dimension();
                OED::Active_Sensors design = OED::Lazy_Greedy_Solve(*oed_problem, data_dim, budget);
                std::cout << "EIG: " << oed_problem->Evaluate(design) << std::endl;
                design.Print_Sensors();
            }
            else
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: Specified design criterion not supported!  Abort!");
            }
        }

    private:
        inline void Create_Model_Interface()
        {
            Teuchos::ParameterList &oed_settings = settings_->sublist("Analysis").sublist("OED");
            this->model_ = OED::makePtr<MrHyDE_Model_Interface<RealT>>(this->comm_, this->settings_, this->solver_, this->postproc_, this->params_, this->rng_);
            
            std::cout << "Driver_MrHyDE::Create_Model_Interface Hello!" << std::endl;
        }

        inline void Create_Observation_Operator_Interface()
        {
            Teuchos::ParameterList &oed_settings = this->settings_->sublist("Analysis").sublist("OED");

            if (!oed_settings.isSublist("observation operator"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified observation operator for OED! Abort!");
            }
            Teuchos::ParameterList &obs_settings = oed_settings.sublist("observation operator");
            std::string &type = obs_settings.get<std::string>("type", "ERROR");
            if (type == "components")
            {
                std::string indices_file = obs_settings.get<std::string>("indices file", "ERROR");
                // TODO: do error checking
                this->observation_operator_ = OED::makePtr<OED::Component_Observation_Operator<RealT>>(this->model_->State_Dimension(), indices_file);
            }
            else if (type == "sensors")
            {
                // TODO: do error checking -- make sure that an objective is already defined
                this->observation_operator_ = OED::makePtr<MrHyDE_Observation_Operator_Interface<RealT>>(this->model_, this->solver_, this->postproc_);
            }

            // determine if the problem is transient
            if (this->solver_->isTransient)
            {
                this->Configure_Transient_Settings();
            }
            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Observation_Operator_Interface Hello!" << std::endl;
        }

        inline void Create_Error_Model_Interface()
        {
            Teuchos::ParameterList &oed_settings = settings_->sublist("Analysis").sublist("OED");

            // TODO: eventually do more here
            if (!oed_settings.isSublist("error model"))
            {
                TEUCHOS_TEST_FOR_EXCEPTION(true, std::runtime_error, "Error: MrHyDE No specified error model for OED! Abort!");
            }
            Teuchos::ParameterList &error_settings = oed_settings.sublist("error model");
            RealT noise_std = error_settings.get<RealT>("std", 0.0);
            int dimension = this->observation_operator_->Data_Dimension();

            this->error_model_ = OED::makePtr<OED::Gaussian_Error<RealT>>(dimension, noise_std);

            std::cout << "OED::MrHyDE_Interface::Driver_MrHyDE::Create_Error_Model_Interface Hello!" << std::endl;
        }

        inline void Create_Prior_Interface()
        {
            Teuchos::ParameterList &oed_settings = settings_->sublist("Analysis").sublist("OED");

            std::vector<std::string> block_names = this->solver_->mesh->getBlockNames();
            this->prior_ = OED::makePtr<MrHyDE_Prior_Interface<RealT>>(this->comm_, *this->settings_, block_names);

            std::cout << "Driver_MrHyDE::Create_Prior_Interface Hello!" << std::endl;
        }

    private:
        inline void Configure_Transient_Settings()
        {
            std::cout << "Driver_MrHyDE::Configure_Transient_Settings begin configuring transient problem" << std::endl;
            std::vector<RealT> meas_times;

            Teuchos::ParameterList &oed_settings = this->settings_->sublist("Analysis").sublist("OED");
            Teuchos::ParameterList &obs_settings = oed_settings.sublist("observation operator");

            // load number of time steps
            int num_steps = this->settings_->sublist("Solver").get<int>("number of steps", 0);
            double initial_time = this->settings_->sublist("Solver").get<double>("initial time", 0.0);
            double final_time = this->settings_->sublist("Solver").get<double>("final time", 1.0);

            // load measurement times
            std::string measurement_times_file = obs_settings.get<std::string>("measurement times", "ERROR");
            this->Load_Measurement_Times(measurement_times_file, meas_times);
            // Start: DELETE THIS LATER
            for (int i = 0; i < meas_times.size(); i++)
            {
                std::cout << "Driver_MrHyDE::Configure_Transient_Settings meas_time=" << meas_times[i] << std::endl;
            }
            // End:   DELETE THIS LATER

            // TODO: update observation operator to the transient operator
            this->observation_operator_ = OED::makePtr<OED::Transient_Observation_Operator<RealT>>(
                this->observation_operator_,
                initial_time,
                final_time,
                num_steps,
                meas_times
            );

            // update dimension numSensors * measTimes
            std::cout << "Driver_MrHyDE::Configure_Transient_Settings finish configuring transient problem data_dim=" << this->observation_operator_->Data_Dimension() << std::endl;
        }

        inline void Load_Measurement_Times(std::string &meas_times_filename, std::vector<double> &meas_times)
        {
            // TODO: load data file, etc.
            double t;
            std::string line;
            std::fstream input_file;
            input_file.open(meas_times_filename, std::ios::in);
            // TODO: add error checking
            if (input_file.is_open())
            {
                while (getline(input_file, line))
                {
                    t = std::stod(line);
                    meas_times.push_back(t);
                }
            }
        }
    };

}

#endif // OEDLIB_DRIVER_MRHYDE_HPP
