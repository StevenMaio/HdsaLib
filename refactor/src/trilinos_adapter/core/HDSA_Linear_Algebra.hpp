#ifndef HDSA_LINEAR_ALGEBRA_HPP
#define HDSA_LINEAR_ALGEBRA_HPP

#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialDenseVector.hpp"
#include "Teuchos_LAPACK.hpp"
#include "Teuchos_SerialDenseSolver.hpp"
#include "Teuchos_SerialSpdDenseSolver.hpp"

namespace HDSA
{

namespace Linear_Algebra
{

  // Compute the SVD of A=U*S*V^T where U and V are orthogonal matrix and S is a diagaonal matrix, stored as a nx1 matrix
  template <class RealT>
  void SVD(const HDSA::Dense_Matrix<RealT> & A, HDSA::Dense_Matrix<RealT> & U, HDSA::Dense_Matrix<RealT> & VT, HDSA::Dense_Matrix<RealT> & S)
  {
    int m = A.numRows();
    int n = A.numCols();
    Teuchos::LAPACK<int, RealT> lapack;
    char JOBU = 'S';
    char JOBVT = 'S';
    HDSA::Ptr<Teuchos::SerialDenseVector<int, RealT> > S_vec = HDSA::makePtr<Teuchos::SerialDenseVector<int, RealT> >(n);
    int LWORK = std::max(1,std::max(3*std::min(m,n)+std::max(m,n),5*std::min(m,n))) + 1;
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > WORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(LWORK,1);
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > RWORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(LWORK,1);
    int info;
    lapack.GESVD(JOBU,JOBVT,m,n,A.Get_Teuchos_Matrix()->values(),m,(*S_vec).values(),U.Get_Teuchos_Matrix()->values(),m,VT.Get_Teuchos_Matrix()->values(),n,(*WORK).values(),LWORK,(*RWORK).values(),&info);
    // Yields the decomposition A = U*diag(S)*VT, note VT is the transpose of V
    for(int k = 0; k < n; k++)
      {
	S.Replace_Element(k,0,(*S_vec)(k));
      }
  }

  // Solve upper triangular system R*x=b for upper triangluar matrix R
  template <class RealT>
  void Upper_Tri_Solve(HDSA::Dense_Matrix<RealT> & x, const HDSA::Dense_Matrix<RealT> & b, const HDSA::Dense_Matrix<RealT> & R)
  {
    int n = R.numRows();
    for(int c = 0; c < x.numCols(); c++)
      {
	for(int k = n-1; k>=0; k--)
	  {
	    RealT val = b(k,c);
	    for(int j = k+1; j < n; j++)
	      {
		val -= x(j,c)*R(k,j);
	      }
	    val = val/R(k,k);
	    x.Replace_Element(k,c,val);
	  }
      }
  }

  // Compute Cholesky factorization of A, A=R^T*R
  template <class RealT>
  void Cholesky_Factorization(const HDSA::Dense_Matrix<RealT> & A, HDSA::Dense_Matrix<RealT> & R)
  {
    int n = A.numRows();
    Teuchos::SerialSpdDenseSolver<int, RealT> Chol_Solve;
    HDSA::Ptr<Teuchos::SerialSymDenseMatrix<int, RealT> > C = HDSA::makePtr<Teuchos::SerialSymDenseMatrix<int, RealT> >(n,n); 
    for(int i = 0; i < n; i++)
      {
	for(int j = 0; j < n; j++)
	  {
	    (*C)(i,j) = A(i,j);
	  }
      }
    Chol_Solve.setMatrix(C);
    Chol_Solve.factor();
    HDSA::Ptr<Teuchos::SerialSymDenseMatrix<int, RealT> > Rc = Chol_Solve.getFactoredMatrix(); // R should be upper triangular. The R returned is symmetric, its upper half is what we need.
    for(int i = 0; i < n; i++)
      {
	for(int j = i; j < n; j++)
	  {
	    R.Replace_Element(i,j,(*Rc)(i,j));
	  }
      }
  }

  // Solve the symmetric linear system via a direct method
  template <class RealT>
  void Symmetric_Direct_Linear_Solve(const HDSA::Dense_Matrix<RealT> & A, HDSA::Dense_Matrix<RealT> & x, const HDSA::Dense_Matrix<RealT> & b) 
  {
    int n = b.numRows();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(A,*R);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,x.numCols());
    for(int c = 0; c < x.numCols(); c++)
      {
	for(int i = 0; i < n; i++)
	  {
	    RealT val = b(i,c);
	    for(int j = 0; j < i; j++)
	      {
		val -= (*y)(j,c)*(*R)(j,i);
	      }
	    val = val/(*R)(i,i);
	    y->Replace_Element(i,c,val);
	  }
      }
    HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(x, *y, *R);
  }
   
  // Compute eigenvalue decomposition A=V*S*V^T for a symmetric matrix A, store eigenvalues in a size nx1 matrix S
  template <class RealT>
  void Symmetric_Eig_Decomposition(const HDSA::Dense_Matrix<RealT> & A, HDSA::Dense_Matrix<RealT> & V, HDSA::Dense_Matrix<RealT> & S) 
  {
    int n = A.numRows();
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > B = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(n,n);
    for(int i = 0; i < n; i++)
      {
	for(int j = 0; j < n; j++)
	  {
	    (*B)(i,j) = A(i,j); 
	  }
      }
    HDSA::Ptr<Teuchos::SerialDenseVector<int, RealT> > S_rev = HDSA::makePtr<Teuchos::SerialDenseVector<int, RealT> >(n); // SYEV outputs eigenvalues if reverse order
    Teuchos::LAPACK<int, RealT> lapack;
    char JOBZ = 'V';
    char UPLO = 'U';
    HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > WORK = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(3*n,3*n);
    int lwork = (*WORK).stride();
    int info;
    lapack.SYEV(JOBZ,UPLO,n,(*B).values(),n,(*S_rev).values(),(*WORK).values(),lwork,&info);
    for(int j = 0; j < n; j++)
      {
	S.Replace_Element(j,0,(*S_rev)(n-1-j));
	for(int i = 0; i < n; i++)
	  {
	    V.Replace_Element(i,j,(*B)(i,n-1-j));
	  }
      }
  }
 
}

}

#endif

