#pragma once

#include "OED_Model_Interface.hpp"
#include "OED_Vector.hpp"
#include "OED_Ptr.hpp"
#include "OED_Random_Number_Generator.hpp"
#include "OED_Tpetra_Vector.hpp"

#include "Tpetra_Map_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT,
              class LO = Tpetra::Map<>::local_ordinal_type,
              class GO = Tpetra::Map<>::global_ordinal_type,
              class Node = Tpetra::Map<>::node_type>
    class MrHyDE_Model_Interface : public OED::Model_Interface<RealT>
    {
    private:
        Teuchos::RCP<MpiComm> comm_;
        Teuchos::RCP<Teuchos::ParameterList> settings_;
        Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solver_;
        Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> postproc_;
        Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> params_;
        OED::Ptr<OED::Trilinos_Adapter::Random_Number_Generator<ScalarT>> rng_;
        int oed_verbosity_{0};

        OED::Ptr<OED::Trilinos_Adapter::Tpetra_Vector<RealT>> parameter_vec_;
        OED::Ptr<OED::Trilinos_Adapter::Tpetra_Vector<RealT>> state_vec_;

        int dim_;

    public:
        MrHyDE_Model_Interface(Teuchos::RCP<MpiComm> &comm,
                               Teuchos::RCP<Teuchos::ParameterList> &settings,
                               Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solver,
                               Teuchos::RCP<MrHyDE::PostprocessManager<SolverNode>> &postproc,
                               Teuchos::RCP<MrHyDE::ParameterManager<SolverNode>> &params,
                               OED::Ptr<OED::Trilinos_Adapter::Random_Number_Generator<ScalarT>> &rng)
            : comm_(comm), settings_(settings), solver_(solver), postproc_(postproc), params_(params), rng_(rng)
        {
            auto tpetra_vec = solver_->linalg->getNewVector(0);
            this->parameter_vec_ = OED::makePtr<OED::Trilinos_Adapter::Tpetra_Vector<RealT>>(tpetra_vec, this->rng_);
            this->state_vec_ = this->parameter_vec_;
            this->dim_ = tpetra_vec->getGlobalLength();
        }

        void State_Solve(OED::Vector<RealT> &u_out, OED::Vector<RealT> &z) override
        {

        }

        void c_u_Transpose_Inverse_Apply(OED::Vector<RealT> &u_out, OED::Vector<RealT> &u_in,
            OED::Vector<RealT> &u, OED::Vector<RealT> &z) override
        {

        }

        void c_z_Transpose_Apply(OED::Vector<RealT> &z_out, OED::Vector<RealT> &u_in,
            OED::Vector<RealT> &u, OED::Vector<RealT> &z) override
        {

        }

        void c_u_Inverse_Apply(OED::Vector<RealT> &u_out, OED::Vector<RealT> &u_in,
            OED::Vector<RealT> &u, OED::Vector<RealT> &z) override
        {

        }

        void c_z_Apply(OED::Vector<RealT> &u_out, OED::Vector<RealT> &z_in,
            OED::Vector<RealT> &u, OED::Vector<RealT> &z) override
        {

        }

        int Param_Dimension() override
        {
            return this->dim_;
        }

        int State_Dimension() override
        {
            return this->dim_;
        }

        Ptr<Vector<RealT>> Get_Empty_Parameter_Vector() override
        {
            return this->parameter_vec_->Clone();
        }

        Ptr<Vector<RealT>> Get_Empty_State_Vector() override
        {
            return this->state_vec_->Clone();
        }

    private:
        OED::Ptr<MrHyDE_OptVector> Map_OED_Vector_to_MrHyDE_OptVector(const OED::Vector<RealT> &z) const
        {
            ROL::Ptr<MrHyDE_OptVector> z_rol;
            if (const OED::Tpetra_Vector<RealT> *ez = dynamic_cast<const OED::Tpetra_Vector<RealT> *>(&z))
            {
            OED::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> fvec = ez->getVector();
            ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<RealT>>(0);
            z_rol = ROL::makePtr<MrHyDE_OptVector>(fvec, svec);
            }
            else if (const OED::Transient_Vector<RealT> *ez = dynamic_cast<const OED::Transient_Vector<RealT> *>(&z))
            {
            std::vector<ROL::Ptr<Tpetra::MultiVector<RealT, LO, GO, SolverNode>>> f_vec;
            std::vector<ROL::Ptr<std::vector<RealT>>> s_vec;
            int n_t = ez->Get_n_t();
            s_vec.resize(n_t);
            for (int k = 0; k < n_t; k++)
            {
                OED::Ptr<OED::Vector<RealT>> z_k = (*ez)[k];
                const OED::Ptr<OED::Std_Vector<RealT>> ez_k = OED::dynamicPtrCast<OED::Std_Vector<RealT>>(z_k);
                s_vec[k] = ez_k->get_std_vec();
            }
            RealT dt = solve_->deltat; // Assumes that z is discretized on the same time nodes as the state
            z_rol = ROL::makePtr<MrHyDE_OptVector>(f_vec, s_vec, dt);
            }
            else if (const OED::Std_Vector<RealT> *ez = dynamic_cast<const OED::Std_Vector<RealT> *>(&z))
            {
            ROL::Ptr<std::vector<ScalarT>> svec = ez->get_std_vec();
            z_rol = ROL::makePtr<MrHyDE_OptVector>(svec);
            }
            return z_rol;
        }
    };

}