#pragma once

#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Vector.hpp"

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