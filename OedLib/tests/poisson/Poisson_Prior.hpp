//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_POISSON_PRIOR_HPP
#define OEDLIB_POISSON_PRIOR_HPP

#include "Eigen/Dense"

#include "OED_Dense_Mass_Matrix.hpp"
#include "OED_Prior_Interface.hpp"
#include "OED_Std_Vector.hpp"
#include "Poisson_Model.hpp"

using Eigen::SelfAdjointEigenSolver;

namespace OED_TEST
{

  template <class RealT>
  class Poisson_Prior : public OED::Prior_Interface<RealT>
  {
  public:
    using Dense_Matrix = Eigen::Matrix<RealT, Eigen::Dynamic, Eigen::Dynamic>;
    using Dense_Vector = Eigen::Matrix<RealT, Eigen::Dynamic, 1>;

  private:
    int param_dim_;
    std::shared_ptr<Dense_Mass_Matrix<RealT>> mass_matrix_;
    Dense_Matrix &M_;
    Dense_Matrix &S_;
    Dense_Matrix L_;
    std::shared_ptr<Eigen::FullPivLU<Dense_Matrix>> L_plu_;
    double norm_scale_;
    double grad_scale_;

  public:
    Poisson_Prior(std::shared_ptr<Poisson_Model<RealT>> &constraint, double norm_scale, double grad_scale)
      : param_dim_(constraint->Param_Dimension()),
        mass_matrix_(constraint->Mass_Matrix()),
        M_(constraint->Mass_Matrix()->M()),
        S_(constraint->S()), L_(param_dim_, param_dim_), 
        norm_scale_(norm_scale), grad_scale_(grad_scale)

    {
      this->L_ = norm_scale * this->M_ + grad_scale * this->S_;
      this->L_plu_ = std::make_shared<Eigen::FullPivLU<Dense_Matrix>>(this->L_);
    }

    void Prior_Precision_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) override
    {

    };

    void Prior_Covariance_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) override
    {
      this->Prior_Covariance_Factor_Apply(z_out, z_in);
      this->Prior_Covariance_Factor_Apply(z_out, z_out);
    };

    void Get_Prior_Mean(Vector<RealT> &z_out) override
    {

    };

    void Prior_Covariance_Factor_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Std_Vector<RealT> &>(z_in);
      auto &z_out_impl = dynamic_cast<Std_Vector<RealT> &>(z_out);
      Dense_Vector &v = z_out_impl.Vec();
      v = this->M_ * z_in_impl.Vec();
      v = this->L_plu_->solve(v);
    };

    void Mass_Matrix_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) override
    {
      auto &z_in_impl = dynamic_cast<Std_Vector<RealT> &>(z_in);
      auto &z_out_impl = dynamic_cast<Std_Vector<RealT> &>(z_out);
      auto &v = z_out_impl.Vec();
      v = this->M_ * z_in_impl.Vec();
    };

    void Mass_Matrix_Inverse_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) override
    {
      this->mass_matrix_->Apply_Inverse(z_out, z_in);
      /*
      auto &z_in_impl = dynamic_cast<Std_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Std_Vector<double> &>(z_out);
      auto &V = this->M_es_.eigenvectors();
      auto &w = this->M_es_.eigenvalues();
      // TODO: compute inverse
      auto &x = z_out_impl.Vec();
      x = V.transpose() * z_in_impl.Vec();
      x = x.cwiseProduct(w.cwiseInverse());
      x = V * x;
      */
    }

    int Param_Dimension() override
    {
      return this->param_dim_;
    }

    std::shared_ptr<Vector<RealT>> Sample_Vector() override
    {
      return nullptr;
    }
  };
}

#endif //OEDLIB_POISSON_PRIOR_HPP
