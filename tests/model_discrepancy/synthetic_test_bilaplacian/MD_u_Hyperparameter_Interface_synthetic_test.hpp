/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_u_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_u_Hyperparameter_Interface_synthetic_test : public HDSA::MD_u_Hyperparameter_Interface<RealT>
{

public:
  std::vector<std::vector<RealT>> Spatial_Domain_Bounds(void) const
  {
    std::vector<std::vector<RealT>> vec; // vec.size() = spatial Dimension, e.g. 1,2, or 3, [ vec[i][0],vec[i][1] ] is an interval bounding the ith spatial coordinate
    vec.resize(1);
    vec[0].resize(2);
    vec[0][0] = 0.0;
    vec[0][1] = 1.0;
    return vec;
  }

  MD_u_Hyperparameter_Interface_synthetic_test() : HDSA::MD_u_Hyperparameter_Interface<RealT>(false, true)
  {
  }

  virtual ~MD_u_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif
