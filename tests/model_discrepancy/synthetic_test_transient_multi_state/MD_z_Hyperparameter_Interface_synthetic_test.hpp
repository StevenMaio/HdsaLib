#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_z_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_z_Hyperparameter_Interface_synthetic_test : public HDSA::MD_z_Hyperparameter_Interface<RealT>
{

private:
  int n_y_;
  int n_t_;
  RealT c_low_;

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

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {
    Transient_Vector<RealT> u_trans = dynamic_cast<Transient_Vector<RealT> &>(u);
    const HDSA_Tpetra_Vector<RealT> z_tpetra = dynamic_cast<const HDSA_Tpetra_Vector<RealT> &>(z);
    Teuchos::ArrayRCP<const RealT> z_view = z_tpetra.getVector()->get1dView();

    RealT coeff_low = 1.0;
    for (int j = 0; j < n_t_; j++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> u_j = u_trans[j];
      HDSA_Tpetra_Vector<RealT> u_tpetra = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(*u_j);
      for (int k = 0; k < n_y_; k++)
      {
        RealT zk = std::pow(z_view[k],3.0);
        u_tpetra.getVector()->replaceGlobalValue(k, 0, coeff_low* zk);
        u_tpetra.getVector()->replaceGlobalValue(n_y_ + k, 0, c_low_ * coeff_low * zk);
      }
      coeff_low = c_low_ * coeff_low;
    }
  }

  MD_z_Hyperparameter_Interface_synthetic_test(const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator, int num_state_solves, int n_y, int n_t, RealT c_low) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, "spatial field",num_state_solves), n_y_(n_y), n_t_(n_t), c_low_(c_low)
  {
  }

  virtual ~MD_z_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif