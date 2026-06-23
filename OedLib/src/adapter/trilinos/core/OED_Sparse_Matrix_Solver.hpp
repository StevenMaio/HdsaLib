/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef OED_SPARSE_MATRIX_SOLVER_HPP
#define OED_SPARSE_MATRIX_SOLVER_HPP

#include "Tpetra_CrsMatrix_decl.hpp"
#include "Amesos2_Factory.hpp"
#include "OED_Sparse_Matrix.hpp"
#include "OED_Incomplete_Chol_Factor.hpp"
#include "OED_Linear_Algebra.hpp"

namespace OED::Trilinos_Adapter
{
  template <class RealT,
            class LO = Tpetra::Map<>::local_ordinal_type,
            class GO = Tpetra::Map<>::global_ordinal_type,
            class Node = Tpetra::Map<>::node_type>
  class Sparse_Matrix_Solver
  {

  private:
    const OED::Ptr<Sparse_Matrix<RealT>> A_;
    const OED::Ptr<OED::Linear_Operator<RealT>> A_op_;
    bool use_direct_;
    OED::Ptr<Amesos2::Solver<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>> solver_;
    OED::Ptr<Incomplete_Chol_Factor<RealT>> L_;
    bool use_incomplete_factorization_;
    int verbosity_;
    std::ostream &out_stream_;
    const std::string solver_type_message_;
    bool use_op_;

  public:
    Sparse_Matrix_Solver(const OED::Ptr<Sparse_Matrix<RealT>> &A, bool use_direct = true, int verbosity = 0, std::ostream &out_stream = std::cout, const std::string solver_type_message = "") : A_(A), A_op_(nullptr), use_direct_(use_direct), verbosity_(verbosity), out_stream_(out_stream), solver_type_message_(solver_type_message)
    {
      if (use_direct_)
      {
        solver_ = Amesos2::create<Tpetra::CrsMatrix<>, Tpetra::MultiVector<>>("KLU2", A_->Get_Tpetra_Matrix());
        solver_->symbolicFactorization();
        solver_->numericFactorization();
      }
      use_incomplete_factorization_ = false;
      use_op_ = false;
    }

    Sparse_Matrix_Solver(const OED::Ptr<OED::Linear_Operator<RealT>> &A_op, int verbosity = 0, std::ostream &out_stream = std::cout, const std::string solver_type_message = "") : A_(nullptr), A_op_(A_op), use_direct_(false), verbosity_(verbosity), out_stream_(out_stream), solver_type_message_(solver_type_message)
    {
      use_incomplete_factorization_ = false;
      use_op_ = true;
    }

    virtual ~Sparse_Matrix_Solver()
    {
    }

    void Set_Incomplete_Factor(OED::Ptr<Incomplete_Chol_Factor<RealT>> &L)
    {
      if (use_direct_)
      {
        OED_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in OED::Sparse_Matrix_Solver: An incomplete factorization was set but the solver is configured for sparse direct solves" << std::endl);
      }
      L_ = L;
      use_incomplete_factorization_ = true;
    }

    bool Use_Incomplete_Factor(void) const
    {
      return use_incomplete_factorization_;
    }

    bool Use_Operator(void) const
    {
      return use_op_;
    }

    std::string Apply_A_Inverse(OED::Vector<RealT> &x, const OED::Vector<RealT> &b)
    {
      std::string output_message;
      std::string output_message_solver;
      if (use_direct_)
      {
        OED::Tpetra_Vector<RealT> &ex = dynamic_cast<OED::Tpetra_Vector<RealT> &>(x);
        const OED::Tpetra_Vector<RealT> &eb = dynamic_cast<const OED::Tpetra_Vector<RealT> &>(b);
        solver_->setX(ex.getVector());
        solver_->setB(eb.getVector());
        solver_->solve();
      }
      else
      {
        RealT tol = 1.0E-10;
        std::string solver = "GMRES";
        if ((A_ && A_->Is_Symmetric()) || (A_op_ && A_op_->Is_Symmetric())) 
        {
          solver = "CG";
        }
        OED::Ptr<A_Operator<RealT>> A_op = OED::makePtr<A_Operator<RealT>>(this);
        if (use_incomplete_factorization_)
        {
          OED::Ptr<OED::Vector<RealT>> b_prec = b.Clone();
          OED::Ptr<OED::Vector<RealT>> x_prec = x.Clone();
          L_->Apply_Inverse(*b_prec, b);
          output_message_solver = OED::Linear_Algebra::Iterative_Linear_Solve<RealT>(*x_prec, *b_prec, *A_op, tol, solver, verbosity_, out_stream_);
          L_->Apply_Inverse_Transpose(x, *x_prec);
        }
        else
        {
          output_message_solver = OED::Linear_Algebra::Iterative_Linear_Solve<RealT>(x, b, *A_op, tol, solver, verbosity_, out_stream_);
        }

        output_message = solver_type_message_ + "::" + output_message_solver;
      }
      return output_message;
    }

    template <class ScalarType>
    class A_Operator : public OED::Linear_Operator<ScalarType>
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

      void Apply(OED::Vector<ScalarType> &y, const OED::Vector<ScalarType> &x) const
      {

        if (A_invert_->Use_Incomplete_Factor())
        {
          OED::Ptr<OED::Vector<RealT>> vec_tmp1 = y.Clone();
          A_invert_->L_->Apply_Inverse_Transpose(*vec_tmp1, x);
          OED::Ptr<OED::Vector<RealT>> vec_tmp2 = y.Clone();
          if (A_invert_->Use_Operator())
          {
            A_invert_->A_op_->Apply(*vec_tmp2, *vec_tmp1);
          }
          else
          {
            A_invert_->A_->Apply(*vec_tmp2, *vec_tmp1);
          }
          A_invert_->L_->Apply_Inverse(y, *vec_tmp2);
        }
        else
        {
          if (A_invert_->Use_Operator())
          {
            A_invert_->A_op_->Apply(y, x);
          }
          else
          {
            A_invert_->A_->Apply(y, x);
          }
        }
      }
    };
  };
}
#endif
