//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_OED_EIG_VECTOR_HPP
#define OEDLIB_OED_EIG_VECTOR_HPP

#include <Eigen/Dense>

#include "OED_Constraint_Interface.hpp"
#include "OED_Vector.hpp"

using Eigen::VectorXd;

namespace OED_TEST
{
  template<class RealT>
  class Test_Vector : public OED::Vector<RealT>
  {
  private:
    VectorXd vec_;
    int dim_;

  public:
    Test_Vector(int dim) : vec_(dim)
    {
      this->dim_ = dim;
    }

    RealT Dot(const OED::Vector<RealT> &x) const override
    {
      // TODO: check to make sure they are compatible
      const Test_Vector<RealT> x_eig = dynamic_cast<const Test_Vector<RealT> &>(x);
      RealT res = 0;
      for (int i = 0; i < this->dim_; i++)
      {
        res += this->vec_[i] * x_eig.vec_[i];
      }
      return res;
    };

    void Scaled_Plus(const RealT alpha, const OED::Vector<RealT> &x) override
    {
      const Test_Vector<RealT> x_eig = dynamic_cast<const Test_Vector<RealT> &>(x);
      for (int i = 0; i < this->dim_; i++)
      {
        this->vec_[i] += alpha * x_eig.vec_[i];
      }
    }

    int Dimension() const override
    {
      return this->dim_;
    }

    void Set_Scalar(const RealT val) override
    {
      for (int i = 0; i < this->dim_; i++)
      {
        this->vec_[i] = val;
      }
    }

    ~Test_Vector() override
    {
    }

    void Set_Entry(int k, RealT val) override
    {
      this->vec_[k] = val;
    }

    RealT Get_Entry(int k) const override
    {
      return this->vec_[k];
    }

    VectorXd &Vec()
    {
      return this->vec_;
    }

    void Set_Vec(const VectorXd &new_vec)
    {
      this->vec_ = new_vec;
    }
  };
}

#endif //OEDLIB_OED_EIG_VECTOR_HPP
