#ifndef OPERATOR_SQRT_TEST_HPP
#define OPERATOR_SQRT_TEST_HPP

#include "HDSA_Operator_Sqrt.hpp"
#include "HDSA_Std_Vector.hpp"

template <class RealT>
class Operator_Sqrt_test : public HDSA::Operator_Sqrt<RealT>
{

private:
  int m_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;

public:
  Operator_Sqrt_test(int m) : m_(m)
  {
    RealT h = 1.0 / static_cast<RealT>(m_ - 1);

    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_->Replace_Element(0, 0, (1.0 / 3.0) * h);
    M_->Replace_Element(0, 1, (1.0 / 6.0) * h);
    for (int i = 1; i < m_ - 1; i++)
    {
      M_->Replace_Element(i, i, (2.0 / 3.0) * h);
      M_->Replace_Element(i, i - 1, (1.0 / 6.0) * h);
      M_->Replace_Element(i, i + 1, (1.0 / 6.0) * h);
    }
    M_->Replace_Element(m_ - 1, m_ - 2, (1.0 / 6.0) * h);
    M_->Replace_Element(m_ - 1, m_ - 1, (1.0 / 3.0) * h);
  }

  void Apply_Operator(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const Std_Vector<RealT> &vec_in_std = dynamic_cast<const Std_Vector<RealT> &>(vec_in);
    Std_Vector<RealT> &vec_out_std = dynamic_cast<Std_Vector<RealT> &>(vec_out);
    for (int k = 0; k < m_; k++)
    {
      b->Replace_Element(k, 0, vec_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    M_->Multiply(*x, *b);
    for (int k = 0; k < m_; k++)
    {
      vec_out_std.Replace_Element(k, (*x)(k, 0));
    }
  }
};

#endif
