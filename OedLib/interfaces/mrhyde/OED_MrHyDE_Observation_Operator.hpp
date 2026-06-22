#pragma once

#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Vector.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT>
    class MrHyDE_Observation_Operator_Interface
        : public OED::Observation_Operator_Interface<RealT>
    {
    public:
        void Observation_Operator_Apply(OED::Vector<RealT> &d_out,
            OED::Vector<RealT> &u_in) override
        {

        }

        void Observation_Operator_Transpose_Apply(OED::Vector<RealT> &u_out,
            OED::Vector<RealT> &d_in) override
        {

        }
    };

}