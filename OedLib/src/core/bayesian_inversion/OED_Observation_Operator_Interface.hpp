#ifndef OEDLIB_OBSERVATION_OPERATOR_INTERFACE_HPP
#define OEDLIB_OBSERVATION_OPERATOR_INTERFACE_HPP

namespace OED
{

    template <class RealT>
    class Observation_Operator_Interface
    {
    public:
        virtual void Observation_Operator_Apply(Vector<RealT> &d_out, Vector<RealT> &u_in) = 0;

        virtual void Observation_Operator_Transpose_Apply(Vector<RealT> &u_out, Vector<RealT> &d_in) = 0;
    };

}

#endif // OEDLIB_OBSERVATION_OPERATOR_INTERFACE_HPP