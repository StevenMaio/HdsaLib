//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_POISSON_CONSTRAINT_HPP
#define OEDLIB_POISSON_CONSTRAINT_HPP
#include <Eigen/src/Core/Matrix.h>

#include "OED_Constraint_Interface.hpp"

using Eigen::MatrixXd;
using OED_TEST::Test_Vector;

namespace OED_TEST
{
  class Poisson_Constraint : public OED::Constraint_Interface<double>
  {
  public:

  private:
    int dim_;
    double h;
    MatrixXd M_;
    MatrixXd M0_;
    MatrixXd S_;
    MatrixXd A_;
    MatrixXd A_inv_;

  public:
    Poisson_Constraint(int dim) :
        dim_{dim}, M_{dim, dim}, M0_{dim, dim},
        S_{dim, dim}, A_{dim, dim}, A_inv_{dim, dim}
    {
      this->h = 1.0 / (dim - 1);
      this->M_(0, 0) = 1.0 / 3.0 * h;
      this->M_(0, 1) = 1.0 / 6.0 * h;
      for (int i = 1; i < dim - 1; i++)
      {
        this->M_(i, i) = (2.0 / 3.0) * h;
        this->M_(i, i - 1) = (1.0 / 6.0) * h;
        this->M_(i, i + 1) = (1.0 / 6.0) * h;
      }
      this->M_(dim - 1, dim - 2) = 1.0 / 6.0 * h;
      this->M_(dim - 1, dim - 1) = 1.0 / 3.0 * h;
      // Copy M and apply BCs
      this->M0_ = this->M_;
      this->M0_(0, 0) = 0;
      this->M0_(0, 1) = 0;
      this->M0_(dim - 1, dim - 2) = 0;
      this->M0_(dim - 1, dim - 1) = 0;

      this->S_(0, 0) = 1.0 / h;
      this->S_(0, 1) = -1.0 / h;
      for (int i = 1; i < dim - 1; i++)
      {
        this->S_(i, i) = 2.0 / h;
        this->S_(i, i - 1) = -1.0 / h;
        this->S_(i, i + 1) = -1.0 / h;
      }
      this->S_(dim - 1, dim - 2) = -1.0 / h;
      this->S_(dim - 1, dim - 1) = 1.0 / h;
      // std::cout << S << std::endl;

      // Apply BCs
      this->A_ = this->S_;
      this->A_(0, 0) = 1.0;
      this->A_(0, 1) = 0.0;
      this->A_(dim -1, dim - 2) = 0.0;
      this->A_(dim -1, dim - 1) = 1.0;
      this->A_inv_ = this->A_.inverse();    // TODO: not good, but this isn't a serious implementation
    }

    MatrixXd &M()
    {
      return this->M_;
    }

    MatrixXd &S()
    {
      return this->S_;
    }

    MatrixXd &A()
    {
      return this->A_;
    }

    int Param_Dimension() override
    {
      return this->dim_;
    }

    int State_Dimension() override
    {
      return this->dim_;
    };

    void State_Solve(OED::Vector<double> &u_out, OED::Vector<double> &z) override
    {
      auto &z_impl = dynamic_cast<Test_Vector<double> &>(z);
      auto &u_out_impl = dynamic_cast<Test_Vector<double> &>(u_out);
      VectorXd rhs = this->M0_ * z_impl.Vec();
      u_out_impl.Vec() = this->A_inv_ * rhs;
    };

    void c_u_Transpose_Inverse_Apply(OED::Vector<double> &u_out, OED::Vector<double> &u_in, OED::Vector<double> &u,
      OED::Vector<double> &z) override
    {
      auto &u_in_impl = dynamic_cast<Test_Vector<double> &>(u_in);
      auto &u_out_impl = dynamic_cast<Test_Vector<double> &>(u_out);
      VectorXd v = this->A_inv_.transpose() * u_in_impl.Vec();
      u_out_impl.Vec() = v;
    };

    void c_z_Transpose_Apply(OED::Vector<double> &z_out, OED::Vector<double> &z_in, OED::Vector<double> &u,
      OED::Vector<double> &z) override
    {
      auto &z_in_impl = dynamic_cast<Test_Vector<double> &>(z_in);
      auto &z_out_impl = dynamic_cast<Test_Vector<double> &>(z_out);
      VectorXd v = -this->M0_.transpose() * z_in_impl.Vec();
      z_out_impl.Vec() = v;
    };

    void c_u_Inverse_Apply(OED::Vector<double> &u_out, OED::Vector<double> &u_in, OED::Vector<double> &u,
      OED::Vector<double> &z) override
    {

    };

    void c_z_Apply(OED::Vector<double> &z_out, OED::Vector<double> &z_in, OED::Vector<double> &u,
      OED::Vector<double> &z) override
    {

    };
  };

}

#endif //OEDLIB_POISSON_CONSTRAINT_H
