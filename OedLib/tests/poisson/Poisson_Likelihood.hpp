//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_POISSON_LIKELIHOOD_HPP
#define OEDLIB_POISSON_LIKELIHOOD_HPP

#include <iterator>

#include "../../src/core/bayesian_inversion/OED_Likelihood_Model.hpp"
#include "../impl/base/vectors/OED_Test_Vector.hpp"

#include <vector>
using OED::Likelihood_Model;

namespace OED_TEST
{
  class Poisson_Likelihood : public Likelihood_Model<double>
  {
  private:
    double noise_std_;
    std::vector<int> observation_vec_;
    int state_dim_;
    int data_dim_;
  public:
    Poisson_Likelihood(int dim, double noise_std, const std::vector<int> &observation_vec) :
        noise_std_{noise_std}, state_dim_{dim}, data_dim_(observation_vec.size())
    {
      for (int i : observation_vec)
      {
        this->observation_vec_.push_back(i);
      }
    }

    int State_Dimension() override
    {
      return this->state_dim_;
    }

    int Data_Dimension() override
    {
      return this->data_dim_;
    }

    void Noise_Precision_Apply(Vector<double> &d_out, Vector<double> &d_in) override
    {
      auto &d_out_impl = dynamic_cast<Test_Vector<double> &>(d_out);
      auto &d_in_impl = dynamic_cast<Test_Vector<double> &>(d_in);
      d_out_impl.Vec() = d_in_impl.Vec();
      d_out_impl.Scale(1 / (this->noise_std_ * this->noise_std_));
    }

    void Observation_Operator_Apply(Vector<double> &d_out, Vector<double> &u_in) override
    {
      auto &u_in_impl = dynamic_cast<Test_Vector<double> &>(u_in);
      auto &d_out_impl = dynamic_cast<Test_Vector<double> &>(d_out);
      for (int i = 0; i < this->observation_vec_.size(); ++i)
      {
        int j = this->observation_vec_[i];
        d_out_impl.Vec()[i] = u_in_impl.Vec()[j];
      }
    }

    void Observation_Operator_Transpose_Apply(Vector<double> &u_out, Vector<double> &d_in) override
    {
      auto &u_out_impl = dynamic_cast<Test_Vector<double> &>(u_out);
      auto &d_in_impl = dynamic_cast<Test_Vector<double> &>(d_in);
      for (int i = 0; i < this->observation_vec_.size(); ++i)
      {
        int j = this->observation_vec_[i];
        u_out_impl.Vec()[j] = d_in_impl.Vec()[i];
      }
    }

    void Get_Observed_Data(Vector<double> &d) override
    {
    }

    void Noise_Covariance_Apply(Vector<double> &d_out, Vector<double> &d_in) override
    {
      auto &d_out_impl = dynamic_cast<Test_Vector<double> &>(d_out);
      auto &d_in_impl = dynamic_cast<Test_Vector<double> &>(d_in);
      d_out_impl.Vec() = d_in_impl.Vec();
      d_out_impl.Scale(this->noise_std_ * this->noise_std_);
    };
  };

}

#endif //OEDLIB_POISSON_LIKELIHOOD_HPP
