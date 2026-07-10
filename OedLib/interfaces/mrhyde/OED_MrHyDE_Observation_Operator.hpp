#pragma once

#include <vector>
#include <algorithm>

#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Vector.hpp"
#include "OED_Tpetra_Vector.hpp"
#include "OED_Sparse_Matrix.hpp"

#include "Tpetra_Map_decl.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT,
              class LO = Tpetra::Map<>::local_ordinal_type,
              class GO = Tpetra::Map<>::global_ordinal_type,
              class Node = Tpetra::Map<>::node_type>
    class MrHyDE_Observation_Operator_Interface
        : public OED::Observation_Operator_Interface<RealT>
    {
    private:
        int state_dim_;
        int data_dim_;
        // OED::Ptr<OED::Trilinos_Adapter::Sparse_Matrix<RealT, LO, GO, Node>> B_;  // Maybe use this later on
        OED::Ptr<OED::Sparse_Matrix<RealT>> B_;

        // Transient Problem Properties
        bool is_transient_{false};
        int num_t_{0};
        double initial_time_{0};
        double final_time_{0};
        OED::Ptr<std::vector<RealT>> meas_times_;
        OED::Ptr<std::vector<int>> meas_indices_;

    public:
        MrHyDE_Observation_Operator_Interface(
            OED::Ptr<MrHyDE_Model_Interface<RealT>> &model,
            Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solver,
            Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> &postproc,
            Teuchos::RCP<Teuchos::ParameterList> &settings
        )
        {
            this->B_ = OED::makePtr<OED::Sparse_Matrix<RealT>>();
            auto u_vec = model->Get_Empty_State_Vector();
            this->state_dim_ = u_vec->Dimension();
            this->data_dim_ = postproc->objectives[0].numSensors;
            this->is_transient_ = solver->isTransient;

            OED::Trilinos_Adapter::Tpetra_Vector<RealT> *u;
            if (!this->is_transient_)
            {
                u = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> *>(&(*u_vec));
            }
            else
            {
                this->Configure_Transient_Settings(settings);
                auto &u_trans = dynamic_cast<OED::Transient_Vector<RealT> &>(*u_vec);
                u = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> *>(&(*u_trans.Get_Vector_Const(0)));
            }

            std::vector<Teuchos::RCP<Tpetra::MultiVector<RealT, LO, GO, Node>>> soln;
            soln.push_back(u->getVector());
            // TODO: need to check later on that this is correct
            auto &response_data = postproc->objectives[0].response_data;

            for (int i = 0; i < u->Dimension(); i++)
            {
                response_data.clear();
                u->Zeros();
                u->Set_Entry(i, 1.0); 
                postproc->computeObjective(soln, 0);
                for (int j = 0; j < response_data.size(); j++)
                {
                    for (int k = 0; k < response_data[j].extent(0); k++)
                    {
                        RealT val = response_data[j](k);
                        this->B_->Add_Entry(k, i, val);
                    }
                }
            }

        }

        void Observation_Operator_Apply(OED::Vector<RealT> &d_out,
            OED::Vector<RealT> &u_in) override
        {
            // TODO: later on, all this should do is cast u_in and d_out appropriately and then move on
            this->B_->Apply(d_out, u_in);
        }

        void Observation_Operator_Transpose_Apply(OED::Vector<RealT> &u_out,
            OED::Vector<RealT> &d_in) override
        {
            this->B_->Transpose_Apply(u_out, d_in);
        }

        int Data_Dimension() override
        {
            return this->data_dim_;
        }

        int State_Dimension() override
        {
            return this->state_dim_;
        }
    
    private:

        inline void Configure_Transient_Settings(Teuchos::RCP<Teuchos::ParameterList> &settings)
        {
            std::cout << "MrHyDE_Observation_Operator_Interface::Configure_Transient_Settings begin configuring transient problem" << std::endl;
            this->meas_times_ = OED::makePtr<std::vector<RealT>>();
            this->meas_indices_ = OED::makePtr<std::vector<int>>();

            Teuchos::ParameterList &oed_settings = settings->sublist("Analysis").sublist("OED");
            Teuchos::ParameterList &obs_settings = oed_settings.sublist("observation operator");

            // load number of time steps
            this->num_t_ = settings->sublist("Solver").get<int>("number of steps", 0) + 1;
            this->initial_time_ = settings->sublist("Solver").get<double>("initial time", 0.0);
            this->final_time_ = settings->sublist("Solver").get<double>("final time", 1.0);

            // load measurement times
            std::string measurement_times_file = obs_settings.get<std::string>("measurement times", "ERROR");
            this->Load_Measurement_Times(measurement_times_file);
            // Start: DELETE THIS LATER
            for (int i = 0; i < this->meas_times_->size(); i++)
            {
                std::cout << "MrHyDE_Observation_Operator_Interface::Configure_Transient_Settings meas_time=" << this->meas_times_->at(i) << std::endl;
            }
            // End:   DELETE THIS LATER

            // TODO: determine measurement indices
            this->Determine_Measurement_Indices();

            // update dimension numSensors * measTimes
            this->data_dim_ *= this->meas_times_->size();
            std::cout << "MrHyDE_Observation_Operator_Interface::Configure_Transient_Settings finish configuring transient problem data_dim=" << this->data_dim_ << std::endl;
        }

        inline void Load_Measurement_Times(std::string &meas_times_filename)
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
                    this->meas_times_->push_back(t);
                }
            }
        }

        inline void Determine_Measurement_Indices()
        {
            double t = this->initial_time_;
            int idx = 0;
            // TODO: iterate over measurement times?
        }
    };

}
