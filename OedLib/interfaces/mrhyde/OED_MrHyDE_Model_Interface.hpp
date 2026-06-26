#pragma once

#include "OED_Model_Interface.hpp"
#include "OED_Vector.hpp"
#include "OED_Std_Vector.hpp"
#include "OED_Ptr.hpp"
#include "OED_Random_Number_Generator.hpp"
#include "OED_Tpetra_Vector.hpp"

#include "Tpetra_Map_decl.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Teuchos_ParameterList.hpp"

#include "solverManager.hpp"

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

            // TODO: this seems to be necessary for getting the gradient
            postproc_->hdsa_solop_data.resize(solver_->setnames.size());
            for (int set = 0; set < solver_->setnames.size(); set++)
            {
                postproc_->hdsa_solop_data[set] = OED::makePtr<MrHyDE::SolutionStorage<SolverNode>>(solver_->settings);
            }
        }

        void State_Solve(OED::Vector<RealT> &u_out, OED::Vector<RealT> &z) override
        {
            ROL::Ptr<MrHyDE_OptVector> z_rol = Map_OED_Vector_to_MrHyDE_OptVector(z);
            params_->updateParams(*z_rol);
            ScalarT val = 0.0;
            solver_->forwardModel(val);
            // TODO: handle transient problems
            OED::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_vec;
            solver_->postproc->soln[0]->extract(u_vec, 0);
            OED::Trilinos_Adapter::Tpetra_Vector<RealT> &eu = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> &>(u_out);
            eu.getVector()->update(1.0, *u_vec, 0.0);
        }

        void State_Transpose_Apply(OED::Vector<RealT> &z_out, OED::Vector<RealT> &u_in, OED::Vector<RealT> &u, OED::Vector<RealT> &z) override
        {
            this->Apply_Solution_Operator_z_Jacobian_Transpose(z_out, u_in, z);
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
        OED::Ptr<MrHyDE_OptVector> Map_OED_Vector_to_MrHyDE_OptVector(OED::Vector<RealT> &z) const
        {
            // TODO: implement transient features
            ROL::Ptr<MrHyDE_OptVector> z_rol;
            if (OED::Trilinos_Adapter::Tpetra_Vector<RealT> *ez = dynamic_cast<OED::Trilinos_Adapter::Tpetra_Vector<RealT> *>(&z))
            {
                OED::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> fvec = ez->getVector();
                ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<RealT>>(0);
                z_rol = ROL::makePtr<MrHyDE_OptVector>(fvec, svec);
            }
            // else if (const OED::Transient_Vector<RealT> *ez = dynamic_cast<const OED::Transient_Vector<RealT> *>(&z))
            // {
            // std::vector<ROL::Ptr<Tpetra::MultiVector<RealT, LO, GO, SolverNode>>> f_vec;
            // std::vector<ROL::Ptr<std::vector<RealT>>> s_vec;
            // int n_t = ez->Get_n_t();
            // s_vec.resize(n_t);
            // for (int k = 0; k < n_t; k++)
            // {
            //     OED::Ptr<OED::Vector<RealT>> z_k = (*ez)[k];
            //     const OED::Ptr<OED::Std_Vector<RealT>> ez_k = OED::dynamicPtrCast<OED::Std_Vector<RealT>>(z_k);
            //     s_vec[k] = ez_k->get_std_vec();
            // }
            // RealT dt = solve_->deltat; // Assumes that z is discretized on the same time nodes as the state
            // z_rol = ROL::makePtr<MrHyDE_OptVector>(f_vec, s_vec, dt);
            // }
            else if (OED::Std_Vector<RealT> *ez = dynamic_cast<OED::Std_Vector<RealT> *>(&z))
            {
                // TODO: need to create a Tpetra multivector and then pass that to the MrHyDE_OptVector constructor
                // ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<RealT>>(0);

                // TODO: this needs to be fixed
                auto &ez_vec = ez->Vec(); 
                ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<ScalarT>>(ez_vec.begin(), ez_vec.end());
                z_rol = ROL::makePtr<MrHyDE_OptVector>(svec);
            }
            return z_rol;
        }

        void Apply_Solution_Operator_z_Jacobian_Transpose(OED::Vector<RealT> &z_out, OED::Vector<RealT> &u_in, OED::Vector<RealT> &z) const
        {
            Write_Data_Solution_Operator(u_in);
            Do_Solution_Operator(true);
            RS_Gradient(z_out, z);
            z_out.Scale(-1.0);
        }

        void Do_Solution_Operator(bool solop_flag) const
        {
            // TODO: what's this doing? I'll leave it on for now
            postproc_->hdsa_solop = solop_flag;
        }

        void Write_Data_Solution_Operator(const OED::Vector<RealT> &u) const
        {
            // if (solver_->isTransient)
            // {
            // const OED::Transient_Vector<RealT> &eu = dynamic_cast<const OED::Transient_Vector<RealT> &>(u);
            // int n_t = solver_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
            // for (int i = 0; i < n_t; i++)
            // {
            //     const OED::Tpetra_Vector<RealT> &eu_i = dynamic_cast<const OED::Tpetra_Vector<RealT> &>(*eu[i]);
            //     RealT currenttime = solver_->initial_time + (double)i * solver_->deltat;
            //     postproc_->hdsa_solop_data[0]->store(eu_i.getVector(), currenttime, 0);
            // }
            // }
            // else
            // {
            const OED::Trilinos_Adapter::Tpetra_Vector<RealT> &eu = dynamic_cast<const OED::Trilinos_Adapter::Tpetra_Vector<RealT> &>(u);
            postproc_->hdsa_solop_data[0]->store(eu.getVector(), 0.0, 0);
            // }
        }


        void RS_Gradient(OED::Vector<RealT> &grad_z, OED::Vector<RealT> &z) const
        {
            bool new_z = Check_New_Params(z);
            if (new_z)
            {
                OED::Ptr<MrHyDE_OptVector> z_rol = this->Map_OED_Vector_to_MrHyDE_OptVector(z);
                MrHyDE_OptVector curr_z = params_->getCurrentVector();
                ROL::Ptr<ROL::Vector<RealT>> z_tmp = curr_z.clone();
                MrHyDE_OptVector ez_tmp = Teuchos::dyn_cast<MrHyDE_OptVector>(dynamic_cast<ROL::Vector<RealT> &>(*z_tmp));
                ez_tmp.set(*z_rol);

                params_->updateParams(ez_tmp);
                ScalarT val = 0.0;
                solver_->forwardModel(val);
            }

            OED::Ptr<MrHyDE_OptVector> grad_z_rol = this->Map_OED_Vector_to_MrHyDE_OptVector(grad_z);
            grad_z_rol->zero();

            solver_->adjointModel(*grad_z_rol);
        }

        bool Check_New_Params(OED::Vector<RealT> &z) const
        {
            OED::Ptr<MrHyDE_OptVector> z_rol = this->Map_OED_Vector_to_MrHyDE_OptVector(z);
            MrHyDE_OptVector curr_z = params_->getCurrentVector();

            ROL::Ptr<ROL::Vector<RealT>> diff = curr_z.clone();
            diff->zero();
            diff->set(curr_z);
            diff->axpy(-1.0, *z_rol);

            ScalarT dnorm = diff->norm();
            ScalarT refnorm = curr_z.norm();
            dnorm = dnorm / refnorm;
            ScalarT reltol = 1.0e-12;

            bool new_z = false;
            if (dnorm > reltol)
            {
                new_z = true;
            }
            return new_z;
        }
    };

}