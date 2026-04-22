/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_ROSENBROCK_TEST_HPP
#define HDSA_ROSENBROCK_TEST_HPP

template <class RealT>
class Rosenbrock
{

private:
  int d_;
  RealT coeff_;

public:
  Rosenbrock(int d) : d_(d)
  {
    coeff_ = 1.0;
  }

  void J(RealT &val, std::vector<RealT> &grad, std::vector<std::vector<RealT>> &hess,
         const std::vector<RealT> &z, const std::vector<RealT> &theta)
  {
    val = 0.0;
    for (int i = 0; i < d_ - 1; i++)
    {
      val += std::pow(z[i] - theta[i], 2.0) + coeff_ * std::pow(z[i + 1] - std::pow(z[i], 2.0), 2.0);
      RealT tmp = 2.0 * (z[i] - theta[i]) - 4.0 * coeff_ * (z[i + 1] - std::pow(z[i], 2.0)) * z[i];
      grad[i] = tmp;
      tmp = 2.0 - 4.0 * coeff_ * (z[i + 1] - 3.0 * std::pow(z[i], 2.0));
      hess[i][i] = tmp;
      tmp = -4.0 * coeff_ * z[i];
      hess[i][i + 1] = tmp;
      if (i > 0)
      {
        tmp = grad[i] + 2 * coeff_ * (z[i] - std::pow(z[i - 1], 2.0));
        grad[i] = tmp;
        tmp = hess[i][i] + 2 * coeff_;
        hess[i][i] = tmp;
        tmp = -4.0 * coeff_ * z[i - 1];
        hess[i][i - 1] = tmp;
      }
    }
    RealT tmp = 2.0 * coeff_ * (z[d_ - 1] - std::pow(z[d_ - 2], 2.0));
    grad[d_ - 1] = tmp;
    tmp = 2.0 * coeff_;
    hess[d_ - 1][d_ - 1] = tmp;
    tmp = -4.0 * coeff_ * z[d_ - 2];
    hess[d_ - 1][d_ - 2] = tmp;
  }

  std::vector<std::vector<RealT>> Compute_B(void)
  {
    std::vector<std::vector<RealT>> B;
    B.resize(d_);
    for (int i = 0; i < d_; i++)
    {
      B[i].resize(d_ - 1);
      if (i < d_ - 1)
      {
        B[i][i] = -2.0;
      }
    }
    return B;
  }
};

#endif
