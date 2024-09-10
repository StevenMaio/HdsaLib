#ifndef HDSA_DENSE_MATRIX_HPP
#define HDSA_DENSE_MATRIX_HPP

#include "Teuchos_SerialDenseMatrix.hpp"

namespace HDSA
{

template <class RealT>
class Dense_Matrix {

  HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > A_;

public:

  // Null constructor
  Dense_Matrix(void)
  { }

  // Constructor given matrix dimensions
  Dense_Matrix(int m, int n) 
  { 
    A_ = HDSA::makePtr<Teuchos::SerialDenseMatrix<int, RealT> >(m,n);
  }

  ~Dense_Matrix()
  { }

  // Number of rows
  int numRows(void) const
  {
    return A_->numRows();
  }

  // Number of columns  
  int numCols(void) const
  {
    return A_->numCols();
  }

  // Access the (i,j) element
  RealT operator () (int i, int j) const
  {
    return (*A_)(i,j);
  }

  // Overwrite the (i,j) element
  void Replace_Element(int i, int j, RealT val)
  {
    (*A_)(i,j) = val;
  }

 // Multiply this*B (a matrix multiply) with options to transpose this and/or B
  void Multiply(HDSA::Dense_Matrix<RealT> & C, const HDSA::Dense_Matrix<RealT> & B, bool A_Trans = false, bool B_Trans = false) const 
  { 
    if(!A_Trans && !B_Trans)
      {
	// No transposes
	C.Get_Teuchos_Matrix()->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*A_,*B.Get_Teuchos_Matrix(),0.0);
      }
    else if(A_Trans && !B_Trans)
      {
	// Transpose A and not B
	C.Get_Teuchos_Matrix()->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*A_,*B.Get_Teuchos_Matrix(),0.0);
      }
    else if(!A_Trans && B_Trans)
      {
	// Transpose B and not A
	C.Get_Teuchos_Matrix()->multiply(Teuchos::NO_TRANS,Teuchos::TRANS,1.0,*A_,*B.Get_Teuchos_Matrix(),0.0);
      }
    else
      {
	// Transpose both A and B
	C.Get_Teuchos_Matrix()->multiply(Teuchos::TRANS,Teuchos::TRANS,1.0,*A_,*B.Get_Teuchos_Matrix(),0.0);
      }
  }

  void zeros(void)
  {
    A_->putScalar(0.0);
  }

  HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > Get_Teuchos_Matrix(void) const
  {
    return A_;
  }

  void Write_to_File(const std::string & name) const
  {
    int m = this->numRows();
    int n = this->numCols();
    std::ofstream fout;
    fout.open(name);
    for(int i = 0; i < m; i++)
      {
	for(int j = 0; j < n; j++)
	{
	  fout << std::setprecision(16) << (*this)(i,j) << "  ";
	}
	fout << "  " << std::endl;
      }
    fout.close();
  }

};

}

#endif
