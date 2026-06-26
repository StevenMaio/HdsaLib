//
// Created by Steven Maio on 5/22/26.
//
#ifndef OEDLIB_CONSTRAINT_HPP
#define OEDLIB_CONSTRAINT_HPP

#include "OED_Vector.hpp"
#include "OED_Ptr.hpp"


namespace OED
{
  template<class RealT>
  class Model_Interface
  {
  public:
    virtual void State_Solve(Vector<RealT> &u_out, Vector<RealT> &z) = 0;

    virtual void State_Transpose_Apply(Vector<RealT> &z_out, Vector<RealT> &u_in, Vector<RealT> &u, Vector<RealT> &z)
    {
      Ptr<Vector<RealT>> u_temp = u.Clone();
      this->c_u_Transpose_Inverse_Apply(*u_temp, u_in, u, z);
      this->c_z_Transpose_Apply(z_out, *u_temp, u, z);
      z_out.Scale(-1);
    }

    virtual void c_u_Transpose_Inverse_Apply(Vector<RealT> &u_out, Vector<RealT> &u_in, Vector<RealT> &u,
                                             Vector<RealT> &z) {}

    virtual void c_z_Transpose_Apply(Vector<RealT> &z_out, Vector<RealT> &u_in, Vector<RealT> &u, Vector<RealT> &z) {}

    virtual void c_u_Inverse_Apply(Vector<RealT> &u_out, Vector<RealT> &u_in, Vector<RealT> &u, Vector<RealT> &z) {}

    virtual void c_z_Apply(Vector<RealT> &u_out, Vector<RealT> &z_in, Vector<RealT> &u, Vector<RealT> &z) {}

    virtual int Param_Dimension() = 0;
    virtual int State_Dimension() = 0;

    virtual ~Model_Interface()
    {
    }

    virtual Ptr<Vector<RealT>> Get_Empty_Parameter_Vector() = 0;

    virtual Ptr<Vector<RealT>> Get_Empty_State_Vector() = 0;

    // TODO: add semi-abstract methods from Sola
    // virtual void Adjoint_Solve(Vector<RealT> &u_out, Vector<RealT> &u, Vector<RealT> &z) = 0;
  };
}

#endif //OEDLIB_CONSTRAINT_HPP
