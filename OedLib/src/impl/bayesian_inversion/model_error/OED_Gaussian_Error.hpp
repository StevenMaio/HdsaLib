//
// Centered and isotropic Gaussian error model
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_POISSON_LIKELIHOOD_HPP
#define OEDLIB_POISSON_LIKELIHOOD_HPP

#include <iterator>

#include "OED_Error_Model_Interface.hpp"

#include <vector>

namespace OED
{
  template <class RealT>
  class Gaussian_Error : public Error_Model_Interface<RealT>
  {
  private:
    // TODO: add some more options
    double noise_std_;
    int data_dim_;

  public:
    Gaussian_Error(int data_dim, double noise_std) :
        data_dim_(data_dim), noise_std_{noise_std} {}

    int Data_Dimension() override
    {
      return this->data_dim_;
    }

    void Noise_Precision_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) override
    {
      d_out.Scaled_Plus(1 / (this->noise_std_ * this->noise_std_), d_in);
    }

    void Noise_Covariance_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) override
    {
      d_out.Scaled_Plus(this->noise_std_ * this->noise_std_, d_in);
    }

  };

}

#endif //OEDLIB_POISSON_LIKELIHOOD_HPP
