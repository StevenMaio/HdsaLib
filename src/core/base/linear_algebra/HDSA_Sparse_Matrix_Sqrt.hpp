#ifndef HDSA_SPARSE_MATRIX_SQRT_HPP
#define HDSA_SPARSE_MATRIX_SQRT_HPP

#include "HDSA_Linear_Operator.hpp"
#include "HDSA_Operator_Sqrt.hpp"
#include "HDSA_Incomplete_Chol_Factor.hpp"

namespace HDSA
{

  template <class RealT>
  class Sparse_Matrix_Sqrt : public HDSA::Operator_Sqrt<RealT>
  {
  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A_;
    HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> L_;
    bool use_incomplete_factorization_;

  public:
    Sparse_Matrix_Sqrt(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &A, bool use_incomplete_factorization = false) : A_(A)
    {
      if (use_incomplete_factorization)
      {
        L_ = HDSA::makePtr<HDSA::Incomplete_Chol_Factor<RealT>>(A_); 
      }
      use_incomplete_factorization_ = use_incomplete_factorization;
    }

    virtual ~Sparse_Matrix_Sqrt()
    {
    }

    void Preconditioner_Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      if (use_incomplete_factorization_)
      {
        L_->Apply_Inverse(vec_out,vec_in);
      }
      else
      {
        vec_out.Set(vec_in);
      }
    }

    void Preconditioner_Transpose_Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      if (use_incomplete_factorization_)
      {
        L_->Apply_Inverse_Transpose(vec_out,vec_in);
      }
      else
      {
        vec_out.Set(vec_in);
      }
    }

    void Preconditioner_Inverse_Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      if (use_incomplete_factorization_)
      {
        L_->Apply(vec_out,vec_in);
      }
      else
      {
        vec_out.Set(vec_in);
      }
    }

    virtual void Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      A_->Apply(vec_out, vec_in);
    }
  };

}

#endif
