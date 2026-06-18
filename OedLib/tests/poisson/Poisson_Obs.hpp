#ifndef OEDLIB_TEST_POISSON_OBS_HPP
#define OEDLIB_TEST_POISSON_OBS_HPP

#include "OED_Observation_Operator_Interface.hpp"

namespace OED_TEST
{
    template <class RealT>
    class Poisson_Observation_Operator : public OED::Observation_Operator_Interface<RealT>
    {
    private:
        std::vector<int> observation_vec_;
        int state_dim_;
        int data_dim_;

    public:
        Poisson_Observation_Operator(int state_dim, std::vector<int> &observation_vec)
            : observation_vec_(observation_vec), state_dim_(state_dim), data_dim_(observation_vec.size())
            {}

        int State_Dimension() { return this->state_dim_; }

        int Data_Dimension() { return this->data_dim_; }

        void Observation_Operator_Apply(Vector<RealT> &d_out, Vector<RealT> &u_in) override
        {
        auto &u_in_impl = dynamic_cast<Std_Vector<RealT> &>(u_in);
        auto &d_out_impl = dynamic_cast<Std_Vector<RealT> &>(d_out);
        for (int i = 0; i < this->observation_vec_.size(); ++i)
        {
            int j = this->observation_vec_[i];
            d_out_impl.Vec()[i] = u_in_impl.Vec()[j];
        }
        }

        void Observation_Operator_Transpose_Apply(Vector<RealT> &u_out, Vector<RealT> &d_in) override
        {
        auto &u_out_impl = dynamic_cast<Std_Vector<RealT> &>(u_out);
        auto &d_in_impl = dynamic_cast<Std_Vector<RealT> &>(d_in);
        for (int i = 0; i < this->observation_vec_.size(); ++i)
        {
            int j = this->observation_vec_[i];
            u_out_impl.Vec()[j] = d_in_impl.Vec()[i];
        }
        }
    };
}

#endif // OEDLIB_TEST_POISSON_OBS_HPP