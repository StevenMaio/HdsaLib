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

    public:
        MrHyDE_Observation_Operator_Interface(
            OED::Ptr<MrHyDE_Model_Interface<RealT>> &model,
            Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solver,
            Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> &postproc
        )
        {
            this->B_ = OED::makePtr<OED::Sparse_Matrix<RealT>>();
            auto u_vec = model->Get_Empty_State_Vector();
            this->state_dim_ = u_vec->Dimension();
            this->data_dim_ = postproc->objectives[0].numSensors;

            OED::Trilinos_Adapter::Tpetra_Vector<RealT> *u;
            if (!solver->isTransient)
            {
                u = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> *>(&(*u_vec));
            }
            else
            {
                auto &u_trans = dynamic_cast<OED::Transient_Vector<RealT> &>(*u_vec);
                u = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> *>(&(*u_trans.Get_Vector_Const(0)));
            }

            std::vector<Teuchos::RCP<Tpetra::MultiVector<RealT, LO, GO, Node>>> soln;
            soln.push_back(u->getVector());
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
    };

}
