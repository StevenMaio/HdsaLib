/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_BF_SOL_OP_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_BF_SOL_OP_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_BF_Sol_Op_Interface.hpp"

template <class RealT>
class BF_Sol_Op_Interface_synthetic_test : public HDSA::BF_Sol_Op_Interface<RealT>
{

private:
  int m_;                                  // Mesh resolution

public:
  BF_Sol_Op_Interface_synthetic_test()
  {
    m_ = 51;
  }

  virtual ~BF_Sol_Op_Interface_synthetic_test()
  {
  }

  // Assume a constraint u = z^3 nodewise on the mesh defined by nodes in x_
  // Assume an objective (1/2)*(u-T)^t*M*(u-T) where T = (x_+1.0)^3 so that the optimal solution is u_opt=(x_+1.0)^3 and z_opt=x_+1.0
  // Assume a high-fidelity model u = 1.2 * z^3

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {
    HDSA::Std_Vector<RealT> u_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u);
    const HDSA::Std_Vector<RealT> z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    for (int k = 0; k < m_; k++)
    {
      u_std.Set_Entry(k, 1.2 * std::pow(z_std(k), 3.0));
    }
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &z) const
  {
    const HDSA::Std_Vector<RealT> u_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u_in);
    const HDSA::Std_Vector<RealT> z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    HDSA::Std_Vector<RealT> z_out_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z_out);
    for (int k = 0; k < m_; k++)
    {
      z_out_std.Set_Entry(k, 1.2 * 3.0 * std::pow(z_std(k), 2.0) * u_in_std(k));
    }
  }
};

#endif
