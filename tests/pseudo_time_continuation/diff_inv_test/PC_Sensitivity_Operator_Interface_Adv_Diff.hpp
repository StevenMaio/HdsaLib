/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef PC_SENSITIVITY_OPERATOR_INTERFACE_ADV_DIFF_HPP
#define PC_SENSITIVITY_OPERATOR_INTERFACE_ADV_DIFF_HPP

#include "HDSA_PC_Euclidean_Sensitivity_Operator_Interface.hpp"

template <class RealT>
class PC_Sensitivity_Operator_Interface_Adv_Diff : public HDSA::PC_Euclidean_Sensitivity_Operator_Interface<RealT>
{

private:
  HDSA::Ptr<Reduced_Space_Objective<RealT>> obj_;

public:
  PC_Sensitivity_Operator_Interface_Adv_Diff(HDSA::Ptr<Reduced_Space_Objective<RealT>> &obj) : obj_(obj)
  {
  }

  virtual ~PC_Sensitivity_Operator_Interface_Adv_Diff()
  {
  }

  void Euclidean_Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    obj_->Gradient(grad, z, theta);
  }

  void Euclidean_Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    obj_->Apply_Hessian(z_out, z_in, z, theta);
  }

  void Euclidean_Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    obj_->Apply_B(z_out, theta_in, z, theta);
  }
};

#endif
