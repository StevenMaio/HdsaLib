//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_POISSON_CONSTRAINT_HPP
#define OEDLIB_POISSON_CONSTRAINT_HPP

#include <memory>

#include "Eigen/Dense"

#include "OED_Std_Vector.hpp"
#include "OED_Model_Interface.hpp"
#include "OED_Dense_Mass_Matrix.hpp"

#include "OED_Ptr.hpp"


namespace OED_TEST
{
  template <class RealT>
  class Poisson_Model : public OED::Model_Interface<RealT>
  {
  public:
    using Dense_Vector = Eigen::VectorXd;
    using Dense_Matrix = Eigen::MatrixXd;

  private:
    int dim_;
    double h;
    Ptr<Dense_Mass_Matrix<RealT>> M_;
    Dense_Matrix S_;
    Dense_Matrix A_;
    Ptr<Eigen::FullPivLU<Dense_Matrix>> A_plu_;

  public:
    Poisson_Model(int dim) :
        dim_{dim}, S_{dim, dim}, A_{dim, dim}
    {
      this->h = 1.0 / (dim - 1);

      Dense_Matrix M(dim, dim);
      M(0, 0) = 1.0 / 3.0 * h;
      M(0, 1) = 1.0 / 6.0 * h;
      for (int i = 1; i < dim - 1; i++)
      {
        M(i, i) = (2.0 / 3.0) * h;
        M(i, i - 1) = (1.0 / 6.0) * h;
        M(i, i + 1) = (1.0 / 6.0) * h;
      }
      M(dim - 1, dim - 2) = 1.0 / 6.0 * h;
      M(dim - 1, dim - 1) = 1.0 / 3.0 * h;
      this->M_ = std::make_shared<Dense_Mass_Matrix<RealT>>(M);

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
      this->A_plu_ = std::make_shared<Eigen::FullPivLU<Dense_Matrix>>(this->A_);
    }

    Ptr<Dense_Mass_Matrix<RealT>> &Mass_Matrix()
    {
      return this->M_;
    }

    Dense_Matrix &S()
    {
      return this->S_;
    }

    Dense_Matrix &A()
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

    void State_Solve(OED::Vector<RealT> &u_out, OED::Vector<RealT> &z_in) override
    {
      auto &z = dynamic_cast<OED::Std_Vector<RealT> &>(z_in);
      auto &u = dynamic_cast<OED::Std_Vector<RealT> &>(u_out);
      Dense_Vector b = z.Vec();
      Dense_Matrix &M = this->M_->M();
      b = M * b;
      b(0) = 0;
      b(this->dim_ - 1) = 0;
      Dense_Vector v = this->A_plu_->solve(b);
      u.Set_Vec(v);
    };

    void State_Adjoint_Apply(OED::Vector<RealT> &m_out, OED::Vector<RealT> &u_in, OED::Vector<RealT> &m, OED::Vector<RealT> &u)
    {
      auto &u_in_impl = dynamic_cast<OED::Std_Vector<RealT> &>(u_in);
      auto &m_out_impl = dynamic_cast<OED::Std_Vector<RealT> &>(m_out);
      Dense_Vector v = this->A_plu_->transpose().solve(u_in_impl.Vec());
      m_out_impl.Set_Vec(v);
    };

    void c_u_Transpose_Inverse_Apply(OED::Vector<RealT> &u_out, OED::Vector<RealT> &u_in, OED::Vector<RealT> &u,
      OED::Vector<RealT> &z) override
    {
      auto &u_in_impl = dynamic_cast<OED::Std_Vector<RealT> &>(u_in);
      auto &u_out_impl = dynamic_cast<OED::Std_Vector<RealT> &>(u_out);
      Dense_Vector v = this->A_plu_->transpose().solve(u_in_impl.Vec());
      u_out_impl.Vec() = v;
    };

    void c_z_Transpose_Apply(OED::Vector<RealT> &z_out, OED::Vector<RealT> &z_in, OED::Vector<RealT> &u,
      OED::Vector<RealT> &z) override
    {
      auto &z_in_impl = dynamic_cast<OED::Std_Vector<RealT> &>(z_in);
      auto &z_out_impl = dynamic_cast<OED::Std_Vector<RealT> &>(z_out);
      // TODO: Replace this with transpose of mass matrix apply and settings BCs to zero
      Dense_Matrix &M = this->M_->M();
      Dense_Vector v = -z_in_impl.Vec();
      v(0) = 0;
      v(this->dim_ - 1) = 0;
      v = M * v;
      z_out_impl.Set_Vec(v);
      //z_out_impl.Set_Vec(v);
      //this->M_->Apply(z_out, z_out);
    };

    void c_u_Inverse_Apply(OED::Vector<RealT> &u_out, OED::Vector<RealT> &u_in, OED::Vector<RealT> &u,
      OED::Vector<RealT> &z) override
    {
      // TODO: implement this
    };

    void c_z_Apply(OED::Vector<RealT> &z_out, OED::Vector<RealT> &z_in, OED::Vector<RealT> &u,
      OED::Vector<RealT> &z) override
    {
      // TODO: implement this
    };

    OED::Ptr<OED::Vector<RealT>> Get_Empty_Parameter_Vector() override
    {
      OED::Ptr<OED::Vector<RealT>> param_vec = OED::makePtr<OED::Std_Vector<RealT>>(this->dim_);
      return param_vec;
    }

    OED::Ptr<OED::Vector<RealT>> Get_Empty_State_Vector() override
    {
      OED::Ptr<OED::Vector<RealT>> state_vec = OED::makePtr<OED::Std_Vector<RealT>>(this->dim_);
      return state_vec;
    }
  };

}

#endif //OEDLIB_POISSON_CONSTRAINT_H
