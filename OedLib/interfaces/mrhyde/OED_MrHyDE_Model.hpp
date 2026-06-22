#pragma once

#include "OED_Model_Interface.hpp"
#include "OED_Vector.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT>
    class MrHyDE_Model_Interface : public OED::Model_Interface<RealT>
    {
    private:
    public:

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

        }

        int State_Dimension() override
        {

        }

    };

}