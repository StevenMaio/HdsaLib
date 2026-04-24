/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef MATRIX_SQRT_TEST_HPP
#define MATRIX_SQRT_TEST_HPP

#include "HDSA_Matrix_Sqrt.hpp"
#include "HDSA_Std_Vector.hpp"

template <class RealT>
class Matrix_Sqrt_test : public HDSA::Matrix_Sqrt<RealT>
{

private:
  int m_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;

public:
  Matrix_Sqrt_test(int m) : m_(m)
  {
    RealT h = 1.0 / static_cast<RealT>(m_ - 1);

    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_->Set_Entry(0, 0, (1.0 / 3.0) * h);
    M_->Set_Entry(0, 1, (1.0 / 6.0) * h);
    for (int i = 1; i < m_ - 1; i++)
    {
      M_->Set_Entry(i, i, (2.0 / 3.0) * h);
      M_->Set_Entry(i, i - 1, (1.0 / 6.0) * h);
      M_->Set_Entry(i, i + 1, (1.0 / 6.0) * h);
    }
    M_->Set_Entry(m_ - 1, m_ - 2, (1.0 / 6.0) * h);
    M_->Set_Entry(m_ - 1, m_ - 1, (1.0 / 3.0) * h);
  }

  void Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const HDSA::Std_Vector<RealT> &vec_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(vec_in);
    HDSA::Std_Vector<RealT> &vec_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(vec_out);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, vec_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    M_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      vec_out_std.Set_Entry(k, (*x)(k, 0));
    }
  }
};

#endif
