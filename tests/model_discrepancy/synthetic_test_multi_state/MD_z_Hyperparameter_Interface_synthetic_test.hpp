#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

template <class RealT>
class MD_z_Hyperparameter_Interface_synthetic_test : public HDSA::MD_z_Hyperparameter_Interface<RealT>
{

public:
  std::vector<std::vector<RealT>> Spatial_Domain_Bounds(void) const
  {
    std::vector<std::vector<RealT>> vec; // vec.size() = spatial dimension, e.g. 1,2, or 3, [ vec[i][0],vec[i][1] ] is an interval bounding the ith spatial coordinate
    vec.resize(1);
    vec[0].resize(2);
    vec[0][0] = 0.0;
    vec[0][1] = 1.0;
    return vec;
  }

  MD_z_Hyperparameter_Interface_synthetic_test(const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, "spatial field")
  {
  }

  virtual ~MD_z_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif