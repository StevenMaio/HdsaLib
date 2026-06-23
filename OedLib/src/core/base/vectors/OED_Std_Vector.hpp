//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_OED_EIG_VECTOR_HPP
#define OEDLIB_OED_EIG_VECTOR_HPP

#include "Eigen/Dense"

#include "OED_Vector.hpp"

namespace OED
{
  template<class RealT>
  class Std_Vector : public OED::Vector<RealT>
  {
  private:
    using Dense_Vector = Eigen::Matrix<RealT, Eigen::Dynamic, 1>;
    Dense_Vector vec_;
    int dim_;

  public:
    Std_Vector(int dim) : vec_(dim)
    {
      this->dim_ = dim;
    }

    RealT Dot(const OED::Vector<RealT> &x) const override
    {
      // TODO: check to make sure they are compatible
      const Std_Vector<RealT> x_eig = dynamic_cast<const Std_Vector<RealT> &>(x);
      RealT res = 0;
      for (int i = 0; i < this->dim_; i++)
      {
        res += this->vec_[i] * x_eig.vec_[i];
      }
      return res;
    };

    void Scaled_Plus(const RealT alpha, const OED::Vector<RealT> &x) override
    {
      const Std_Vector<RealT> x_eig = dynamic_cast<const Std_Vector<RealT> &>(x);
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

    ~Std_Vector() override
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

    Dense_Vector  &Vec()
    {
      return this->vec_;
    }

    void Set_Vec(const Dense_Vector &new_vec)
    {
      this->vec_ = new_vec;
    }

    void Randomize_Standard_Norm() override {
      OED_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in OED::Vector: Randomize_Standard_Norm has not been implemented for this vector type" << std::endl);
    }

    Ptr<Std_Vector<RealT>> Clone() const override {
      Ptr<Std_Vector<RealT>> clone = OED<Std_Vector<RealT>>(this->dim_);
      clone->Set_Vec(this->vec_);
    }
  };
}

#endif //OEDLIB_OED_EIG_VECTOR_HPP
