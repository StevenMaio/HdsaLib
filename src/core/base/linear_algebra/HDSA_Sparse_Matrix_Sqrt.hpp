#ifndef HDSA_SPARSE_MATRIX_SQRT_HPP
#define HDSA_SPARSE_MATRIX_SQRT_HPP

#include "HDSA_Linear_Operator.hpp"
#include "HDSA_Operator_Sqrt.hpp"

namespace HDSA
{

  template <class RealT>
  class Sparse_Matrix_Sqrt : public HDSA::Operator_Sqrt<RealT>
  {
  private:
    const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> A_;

  public:
    Sparse_Matrix_Sqrt(const HDSA::Ptr<HDSA::Sparse_Matrix<RealT>> &A) : A_(A)
    {
    }

    virtual ~Sparse_Matrix_Sqrt()
    {
    }

    virtual void Apply(HDSA::Vector<RealT> &vec_out, const HDSA::Vector<RealT> &vec_in) const
    {
      A_->Apply(vec_out, vec_in);
    }
  };

}

#endif
