#ifndef OED_COMPONENT_OBSERVATION_OPERATOR_HPP
#define OED_COMPONENT_OBSERVATION_OPERATOR_HPP

#include "OED_Observation_Operator_Interface.hpp"

namespace OED
{
    template <class RealT>
    class Component_Observation_Operator : public Observation_Operator_Interface<RealT>
    {
    private:
        std::vector<int> observation_vec_;
        int state_dim_;
        int data_dim_;

    public:
        Component_Observation_Operator(int state_dim, std::vector<int> &observation_vec)
            : observation_vec_(observation_vec), state_dim_(state_dim), data_dim_(observation_vec.size())
            {}

        int State_Dimension() { return this->state_dim_; }

        int Data_Dimension() { return this->data_dim_; }

        void Observation_Operator_Apply(Vector<RealT> &d_out, Vector<RealT> &u_in) override
        {
            for (int i = 0; i < this->observation_vec_.size(); ++i)
            {
                int j = this->observation_vec_[i];
                d_out.Set_Entry(i, u_in.Get_Entry(j));
            }
        }

        void Observation_Operator_Transpose_Apply(Vector<RealT> &u_out, Vector<RealT> &d_in) override
        {
            for (int i = 0; i < this->observation_vec_.size(); ++i)
            {
                int j = this->observation_vec_[i];
                u_out.Set_Entry(j, d_in.Get_Entry(i));
            }
        }
    };
}

#endif // OED_COMPONENT_OBSERVATION_OPERATOR_HPP