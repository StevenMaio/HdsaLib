#ifndef HDSA_SPARSE_MATRIX_SOLVER_HPP
#define HDSA_SPARSE_MATRIX_SOLVER_HPP

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Amesos2_Factory.hpp"

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

  public:
    Sparse_Matrix_Solver(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &A, bool use_direct = true) : A_(A), use_direct_(use_direct)
    {
      if (use_direct_)
      {
        solver_ = Amesos2::create<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>("KLU2", A_->Get_Tpetra_Matrix());
        solver_->symbolicFactorization();
        solver_->numericFactorization();
      }
    }

    virtual ~Sparse_Matrix_Solver()
    {
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
        std::string solver = "CG";
        bool verbose = false;
        HDSA::Ptr<A_Operator<RealT>> A_op = HDSA::makePtr<A_Operator<RealT>>(this);
        HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(x, b, *A_op, tol, solver, verbose);
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

      void matvec(HDSA::Vector<ScalarType> &y, const HDSA::Vector<ScalarType> &x) const
      {
        A_invert_->A_->Apply(y, x);
      }
    };
  };
}
#endif
