#ifndef HDSA_SPARSE_MATRIX_SOLVER_HPP
#define HDSA_SPARSE_MATRIX_SOLVER_HPP

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Amesos2_Factory.hpp"
#include "HDSA_Sparse_Matrix.hpp"
#include "HDSA_Incomplete_Chol_Factor.hpp"
#include "HDSA_Linear_Algebra.hpp"

namespace HDSA
{
  template <class RealT,
            class LO = Tpetra::Map<>::local_ordinal_type,
            class GO = Tpetra::Map<>::global_ordinal_type,
            class Node = Tpetra::Map<>::node_type>
  class Sparse_Matrix_Solver
  {

  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A_;
    bool use_direct_;
    HDSA::Ptr<Amesos2::Solver<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>> solver_;
    HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> L_;
    bool use_incomplete_factorization_;
    bool verbose_;

  public:
    Sparse_Matrix_Solver(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &A, bool use_direct = true, bool verbose = false) : A_(A), use_direct_(use_direct)
    {
      if (use_direct_)
      {
        solver_ = Amesos2::create<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>("KLU2", A_->Get_Tpetra_Matrix());
        solver_->symbolicFactorization();
        solver_->numericFactorization();
      }
      use_incomplete_factorization_ = false;
      verbose_ = verbose;
    }

    virtual ~Sparse_Matrix_Solver()
    {
    }

    void Set_Incomplete_Factor(HDSA::Ptr<HDSA::Incomplete_Chol_Factor<RealT>> &L)
    {
      if (use_direct_)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA::Sparse_Matrix_Solver: An incomplete factorization was set but the solver is configured for sparse direct solves" << std::endl);
      }
      L_ = L;
      use_incomplete_factorization_ = true;
    }

    bool Use_Incomplete_Factor(void) const
    {
      return use_incomplete_factorization_;
    }

    void Apply_A_Inverse(HDSA::Vector<RealT> &x, const HDSA::Vector<RealT> &b)
    {
      if (use_direct_)
      {
        HDSA::Tpetra_Vector<RealT> &ex = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(x);
        const HDSA::Tpetra_Vector<RealT> &eb = dynamic_cast<const HDSA::Tpetra_Vector<RealT> &>(b);
        solver_->setX(ex.getVector());
        solver_->setB(eb.getVector());
        solver_->solve();
      }
      else
      {
        RealT tol = 1.0E-10;
        std::string solver = "GMRES";
        if (A_->Is_Symmetric())
        {
          solver = "CG";
        }
        HDSA::Ptr<A_Operator<RealT>> A_op = HDSA::makePtr<A_Operator<RealT>>(this);
        if (use_incomplete_factorization_)
        {
          HDSA::Ptr<HDSA::Vector<RealT>> b_prec = b.Clone();
          HDSA::Ptr<HDSA::Vector<RealT>> x_prec = x.Clone();
          L_->Apply_Inverse(*b_prec, b);
          HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(*x_prec, *b_prec, *A_op, tol, solver, verbose_);
          L_->Apply_Inverse_Transpose(x, *x_prec);
        }
        else
        {
          HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(x, b, *A_op, tol, solver, verbose_);
        }
      }
    }

    template <class ScalarType>
    class A_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      const Sparse_Matrix_Solver<ScalarType> *A_invert_;

    public:
      A_Operator(const Sparse_Matrix_Solver<ScalarType> *A_invert) : A_invert_(A_invert)
      {
      }

      ~A_Operator()
      {
      }

      void Apply(HDSA::Vector<ScalarType> &y, const HDSA::Vector<ScalarType> &x) const
      {
        if (A_invert_->Use_Incomplete_Factor())
        {
          HDSA::Ptr<HDSA::Vector<RealT>> vec_tmp1 = y.Clone();
          A_invert_->L_->Apply_Inverse_Transpose(*vec_tmp1, x);
          HDSA::Ptr<HDSA::Vector<RealT>> vec_tmp2 = y.Clone();
          A_invert_->A_->Apply(*vec_tmp2, *vec_tmp1);
          A_invert_->L_->Apply_Inverse(y, *vec_tmp2);
        }
        else
        {
          A_invert_->A_->Apply(y, x);
        }
      }
    };
  };
}
#endif
