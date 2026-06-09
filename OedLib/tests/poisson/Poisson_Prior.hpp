//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_POISSON_PRIOR_HPP
#define OEDLIB_POISSON_PRIOR_HPP

#include "Eigen/Dense"

#include "OED_Prior_Interface.hpp"
#include "OED_Test_Vector.hpp"
#include "Poisson_Constraint.hpp"

using Eigen::MatrixXd;
using Eigen::SelfAdjointEigenSolver;
using Eigen::FullPivLU;

namespace OED_TEST
{

  class Poisson_Prior : public OED::Prior_Interface<double>
  {
  private:
    int param_dim_;
    MatrixXd &M_;
    SelfAdjointEigenSolver<MatrixXd> M_es_;
    MatrixXd &S_;
    MatrixXd L_;
    std::shared_ptr<FullPivLU<MatrixXd>> L_plu_;
    double norm_scale_;
    double grad_scale_;

  public:
    Poisson_Prior(std::shared_ptr<Poisson_Constraint> &constraint, double norm_scale, double grad_scale)
      : param_dim_(constraint->Param_Dimension()), M_(constraint->M()), M_es_(constraint->M()), 
        S_(constraint->S()), L_(param_dim_, param_dim_), 
        norm_scale_(norm_scale), grad_scale_(grad_scale)

    {
      this->L_ = norm_scale * this->M_ + grad_scale * this->S_;
      this->L_plu_ = std::make_shared<FullPivLU<MatrixXd>>(this->L_);
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
      VectorXd &v = z_out_impl.Vec();
      v = this->M_ * z_in_impl.Vec();
      v = this->L_plu_->solve(v);
    };

    void Mass_Matrix_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      auto &v = z_out_impl.Vec();
      v = this->M_ * z_in_impl.Vec();
    };

    void Mass_Matrix_Inverse_Apply(Vector<double> &z_out, Vector<double> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      auto &V = this->M_es_.eigenvectors();
      auto &w = this->M_es_.eigenvalues();
      // TODO: compute inverse
      auto &x = z_out_impl.Vec();
      x = V.transpose() * z_in_impl.Vec();
      x = x.cwiseProduct(w.cwiseInverse());
      x = V * x;
    }

    int Param_Dimension() override
    {
      return this->param_dim_;
    }

    std::shared_ptr<Vector<double>> Sample_Vector()
    {
      return nullptr;
    }
  };
}

#endif //OEDLIB_POISSON_PRIOR_HPP
