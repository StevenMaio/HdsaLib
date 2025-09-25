#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_z_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_z_Hyperparameter_Interface_synthetic_test : public HDSA::MD_z_Hyperparameter_Interface<RealT>
{

public:

  MD_z_Hyperparameter_Interface_synthetic_test(const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, "spatial field")
  {
  }

  virtual ~MD_z_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif