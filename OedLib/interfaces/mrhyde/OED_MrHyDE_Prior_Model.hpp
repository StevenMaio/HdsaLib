#pragma once

#include "OED_Prior_Interface.hpp"
#include "OED_Vector.hpp"

namespace OED::MrHyDE_Interface
{

    template <class RealT>
    class MrHyDE_Prior_Interface : public OED::Prior_Interface<RealT>
    {
    public:

        void Prior_Precision_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        void Prior_Covariance_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        void Get_Prior_Mean(OED::Vector<RealT> &m_out) override
        {

        }

        void Prior_Covariance_Factor_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        int Param_Dimension() override
        {
            // TODO: actually implement this
            return 0;
        }

        // TODO: move this at some point to a potentially new class
        void Mass_Matrix_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        void Mass_Matrix_Inverse_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &m_in) override
        {

        }

        Ptr<OED::Vector<RealT>> Sample_Vector() override
        {
            // TODO: do something like a multi-vector I suppose
        }

    }
}