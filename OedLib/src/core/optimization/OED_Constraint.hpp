//
// Created by Steven Maio on 5/22/26.
//
#ifndef OEDLIB_CONSTRAINT_HPP
#define OEDLIB_CONSTRAINT_HPP
#include "../base/vectors/OED_Vector.hpp"


namespace OED
{
  template<class RealT>
  class Constraint
  {
  public:
    virtual void State_Solve(Vector<RealT> &u_out, Vector<RealT> &z) = 0;

    virtual void c_u_Transpose_Inverse_Apply(Vector<RealT> &u_out, Vector<RealT> &u_in, Vector<RealT> &u,
                                             Vector<RealT> &z) = 0;

    virtual void c_z_Transpose_Apply(Vector<RealT> &z_out, Vector<RealT> &u_in, Vector<RealT> &u, Vector<RealT> &z) = 0;

    virtual void c_u_Inverse_Apply(Vector<RealT> &u_out, Vector<RealT> &u_in, Vector<RealT> &u, Vector<RealT> &z) = 0;

    virtual void c_z_Apply(Vector<RealT> &u_out, Vector<RealT> &z_in, Vector<RealT> &u, Vector<RealT> &z) = 0;

    virtual int Param_Dimension() = 0;
    virtual int State_Dimension() = 0;

    virtual ~Constraint()
    {
    }

    // TODO: add semi-abstract methods from Sola
    // virtual void Adjoint_Solve(Vector<RealT> &u_out, Vector<RealT> &u, Vector<RealT> &z) = 0;
  };
}

#endif //OEDLIB_CONSTRAINT_HPP
