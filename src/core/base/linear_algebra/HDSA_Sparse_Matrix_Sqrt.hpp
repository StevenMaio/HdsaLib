/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_SPARSE_MATRIX_SQRT_HPP
#define HDSA_SPARSE_MATRIX_SQRT_HPP

#include "HDSA_Linear_Operator.hpp"
#include "HDSA_Matrix_Sqrt.hpp"
#include "HDSA_Incomplete_Chol_Factor.hpp"

namespace HDSA
{

  template <class RealT>
  class Sparse_Matrix_Sqrt : public HDSA::Matrix_Sqrt<RealT>
  {
  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A_;
    const HDSA::Ptr<HDSA::Linear_Operator<RealT>> A_op_;
    HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> L_;
    bool use_incomplete_factorization_;
    bool use_op_;

  public:
    Sparse_Matrix_Sqrt(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &A, const std::string solver_type_message = "") : HDSA::Matrix_Sqrt<RealT>(solver_type_message), A_(A), A_op_(nullptr)
    {
      use_incomplete_factorization_ = false;
      use_op_ = false;
    }

    Sparse_Matrix_Sqrt(const HDSA::Ptr<HDSA::Linear_Operator<RealT>> &A_op, const std::string solver_type_message = "") : HDSA::Matrix_Sqrt<RealT>(solver_type_message), A_(nullptr), A_op_(A_op)
    {
      use_incomplete_factorization_ = false;
      use_op_ = true;
    }

    virtual ~Sparse_Matrix_Sqrt()
    {
    }

    void Set_Incomplete_Factor(HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> &L)
    {
      L_ = L;
      use_incomplete_factorization_ = true;
    }

    void Disable_Incomplete_Factorization(void)
    {
      use_incomplete_factorization_ = false;
    }

    void Preconditioner_Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      if (use_incomplete_factorization_)
      {
        L_->Apply_Inverse(vec_out, vec_in);
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
        L_->Apply_Inverse_Transpose(vec_out, vec_in);
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
        L_->Apply(vec_out, vec_in);
      }
      else
      {
        vec_out.Set(vec_in);
      }
    }

    virtual void Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      if (use_op_)
      {
        A_op_->Apply(vec_out, vec_in);
      }
      else
      {
        A_->Apply(vec_out, vec_in);
      }
    }
  };

}

#endif
