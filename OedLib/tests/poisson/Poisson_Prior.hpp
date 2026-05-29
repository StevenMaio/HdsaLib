//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_POISSON_PRIOR_HPP
#define OEDLIB_POISSON_PRIOR_HPP

#include <Eigen/Dense>

#include "../../src/core/bayesian_inversion/OED_Prior_Model.hpp"
#include "../impl/base/vectors/OED_Test_Vector.hpp"
#include "Poisson_Constraint.hpp"

using Eigen::MatrixXd;
using Eigen::FullPivLU;

namespace OED_TEST
{

  class Poisson_Prior : public OED::Prior_Model<double>
  {
  private:
    int param_dim_;
    MatrixXd &M_;
    MatrixXd M_inv_;  // symmetric.. can use Cholesky?
    MatrixXd &S_;
    MatrixXd L_;
    MatrixXd L_inv_;
    double norm_scale_;
    double grad_scale_;
  public:
    Poisson_Prior(Poisson_Constraint &constraint, double norm_scale, double grad_scale)
      : param_dim_(constraint.Param_Dimension()), M_(constraint.M()), M_inv_(this->param_dim_, this->param_dim_),
        S_(constraint.S()), L_(param_dim_, param_dim_), L_inv_(param_dim_, param_dim_),
        norm_scale_(norm_scale), grad_scale_(grad_scale)

    {
      this->M_inv_ = this->M_.inverse(); // TODO: I know this is bad
      this->L_ = norm_scale * this->M_ + grad_scale * this->S_;
      this->L_inv_ = this->L_.inverse();  // TODO: I know this is also bad (just for now)
    }

    void Prior_Precision_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {

    };

    void Prior_Covariance_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      this->Prior_Factor_Apply(z_out, z_in);
      this->Prior_Factor_Apply(z_out, z_out);
    };

    void Get_Prior_Mean(Vector<double> &z_out) override
    {

    };

    void Prior_Factor_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      VectorXd v = this->M_ * z_in_impl.Vec();
      v = this->L_inv_ * v;
      z_out_impl.Vec() = v;
    };

    void Mass_Matrix_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      VectorXd v = this->M_ * z_in_impl.Vec();
      z_out_impl.Vec() = v;
    };

    void Mass_Matrix_Inverse_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      VectorXd v = this->M_inv_ * z_in_impl.Vec();
      z_out_impl.Vec() = v;
    }

    int Param_Dimension() override
    {
      return this->param_dim_;
    }
  };
}

#endif //OEDLIB_POISSON_PRIOR_HPP
