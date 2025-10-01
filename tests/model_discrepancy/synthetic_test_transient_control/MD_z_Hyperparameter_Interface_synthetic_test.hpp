#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_SYNTHETIC_TEST_HPP

#include "HDSA_MD_z_Hyperparameter_Interface.hpp"

template <class RealT>
class MD_z_Hyperparameter_Interface_synthetic_test : public HDSA::MD_z_Hyperparameter_Interface<RealT>
{

private:
  int n_y_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x_;

public:
  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const override
  {
    HDSA::Transient_Vector<RealT> u_trans = dynamic_cast<HDSA::Transient_Vector<RealT> &>(u);
    const HDSA::Transient_Vector<RealT> z_trans = dynamic_cast<const HDSA::Transient_Vector<RealT> &>(z);
    int n_t = u_trans.Get_n_t();
    for (int j = 0; j < n_t; j++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> uj = u_trans[j];
      HDSA::Ptr<HDSA::Vector<RealT>> zj = z_trans[j];
      HDSA::Tpetra_Vector<RealT> uj_tpetra = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(*uj);
      HDSA::Std_Vector<RealT> zj_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*zj);
      for (int k = 0; k < n_y_; k++)
      {
        RealT val = zj_std.Get_Entry(0) * (1.0 - (*x_)(k, 0)) + zj_std.Get_Entry(1) * (*x_)(k, 0);
        uj_tpetra.getVector()->replaceGlobalValue(k, 0, val);
      }
    }
  }

  MD_z_Hyperparameter_Interface_synthetic_test(const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, "transient vector", 100)
  {
    n_y_ = 50;
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(n_y_, 1);
    for (int k = 0; k < n_y_; k++)
    {
      x_->Set_Entry(k, 0, static_cast<RealT>(k) / static_cast<RealT>(n_y_ - 1));
    }
  }

  virtual ~MD_z_Hyperparameter_Interface_synthetic_test()
  {
  }
};

#endif