/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef RANDOMIZED_GEVP_TEST_HPP
#define RANDOMIZED_GEVP_TEST_HPP

#include "HDSA_Randomized_GEVP.hpp"

template <class RealT>
class Randomized_GEVP_test : public HDSA::Randomized_GEVP<RealT>
{

private:
  int m_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> A_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;

public:
  Randomized_GEVP_test(HDSA::Vector<RealT> &vec, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator)
      : HDSA::Randomized_GEVP<RealT>(), random_number_generator_(random_number_generator)
  {
    m_ = vec.Dimension();
    RealT h = 1.0 / static_cast<RealT>(m_ - 1);

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    A_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);

    S_->Set_Entry(0, 0, 1.0 / h);
    S_->Set_Entry(0, 1, -1.0 / h);
    for (int i = 1; i < m_ - 1; i++)
    {
      S_->Set_Entry(i, i, 2.0 / h);
      S_->Set_Entry(i, i - 1, -1.0 / h);
      S_->Set_Entry(i, i + 1, -1.0 / h);
    }
    S_->Set_Entry(m_ - 1, m_ - 2, -1.0 / h);
    S_->Set_Entry(m_ - 1, m_ - 1, 1.0 / h);

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

    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = (1.e-2) * (*S_)(i, j) + (*M_)(i, j);
        A_->Set_Entry(i, j, val);
      }
    }
  }

  void Apply_Operator(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const HDSA::Std_Vector<RealT> &vec_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(vec_in);
    HDSA::Std_Vector<RealT> &vec_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(vec_out);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, vec_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*A_, *x, *b);
    for (int k = 0; k < m_; k++)
    {
      vec_out_std.Set_Entry(k, (*x)(k, 0));
    }
  }

  void Apply_Weighting_Operator(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
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

  void Apply_Weighting_Operator_Inverse(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    const HDSA::Std_Vector<RealT> &vec_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(vec_in);
    HDSA::Std_Vector<RealT> &vec_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(vec_out);
    for (int k = 0; k < m_; k++)
    {
      b->Set_Entry(k, 0, vec_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*M_, *x, *b);
    for (int k = 0; k < m_; k++)
    {
      vec_out_std.Set_Entry(k, (*x)(k, 0));
    }
  }

  void Generate_Random_Samples(HDSA::MultiVector<RealT> &samples) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> R = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Linear_Algebra::Cholesky_Factorization(*M_, *R);

    int num_samples = samples.Number_of_Vectors();
    for (int i = 0; i < num_samples; i++)
    {

      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      HDSA::Ptr<HDSA::Vector<RealT>> vec_in = samples[0]->Clone();
      HDSA::Std_Vector<RealT> &vec_in_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*vec_in);
      vec_in_std.Randomize_Standard_Normal();
      for (int k = 0; k < m_; k++)
      {
        b->Set_Entry(k, 0, vec_in_std(k));
      }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
      HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*x, *b, *R);
      HDSA::Std_Vector<RealT> &vec_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(*samples[i]);
      for (int k = 0; k < m_; k++)
      {
        vec_out_std.Set_Entry(k, (*x)(k, 0));
      }
    }
  }
};

#endif
