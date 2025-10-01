#ifndef HDSA_PC_LIS_INTERFACE_ADV_DIFF_HPP
#define HDSA_PC_LIS_INTERFACE_ADV_DIFF_HPP

#include "HDSA_PC_LIS_Interface.hpp"

template <class RealT>
class PC_LIS_Interface_Adv_Diff : public HDSA::PC_LIS_Interface<RealT>
{

private:
  HDSA::Ptr<Reduced_Space_Objective<RealT>> obj_;

public:
  PC_LIS_Interface_Adv_Diff(HDSA::Ptr<Reduced_Space_Objective<RealT>> &obj) : obj_(obj)
  {
  }

  virtual ~PC_LIS_Interface_Adv_Diff()
  {
  }

  void Apply_Misfit_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    obj_->Apply_Misfit_Hessian(z_out, z_in, z, theta);
  }

  void Apply_Prior_Precision(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    obj_->prior_and_like_->Apply_Prior_Precision(z_out, z_in);
  }

  void Apply_Prior_Covariance(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
  {
    obj_->prior_and_like_->Apply_Prior_Covariance(z_out, z_in);
  }

  void Generate_Prior_Samples(HDSA::MultiVector<RealT> &samples) const
  {
    obj_->prior_and_like_->Generate_Prior_Samples(samples);
  }
};

#endif
