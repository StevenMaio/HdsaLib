#ifndef HDSA_DENSE_MATRIX_HPP
#define HDSA_DENSE_MATRIX_HPP

#include <fstream>

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

  // Clone the given matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Clone(void) const
  {
    int m = numRows();
    int n = numCols();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<Dense_Matrix<RealT> >(m,n);
    return C;
  }

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

  // Access the (i,j) element reference
  RealT* Get_Element_Ptr(int i, int j) const
  {
    return &(*A_)(i,j);
  }

  // Overwrite the (i,j) element
  void Replace_Element(int i, int j, RealT val)
  {
    (*A_)(i,j) = val;
  }
 
  // Multiply this*B (a matrix multiply) with options to transpose this and/or B
  void Multiply(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & C, HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & B, bool A_Trans = false, bool B_Trans = false) const 
  {
    int Am = numRows();
    int An = numCols();
    int Bm = B->numRows();
    int Bn = B->numCols();
    
    if(!A_Trans && !B_Trans)
      {
	// No transposes
	C->Get_Teuchos_Matrix()->multiply(Teuchos::NO_TRANS,Teuchos::NO_TRANS,1.0,*A_,*B->Get_Teuchos_Matrix(),0.0);
      }
    else if(A_Trans && !B_Trans)
      {
	// Transpose A and not B
	C->Get_Teuchos_Matrix()->multiply(Teuchos::TRANS,Teuchos::NO_TRANS,1.0,*A_,*B->Get_Teuchos_Matrix(),0.0);
      }
    else if(!A_Trans && B_Trans)
      {
	// Transpose B and not A
	C->Get_Teuchos_Matrix()->multiply(Teuchos::NO_TRANS,Teuchos::TRANS,1.0,*A_,*B->Get_Teuchos_Matrix(),0.0);
      }
    else
      {
	// Transpose both A and B
	C->Get_Teuchos_Matrix()->multiply(Teuchos::TRANS,Teuchos::TRANS,1.0,*A_,*B->Get_Teuchos_Matrix(),0.0);
      }
  }

  // Multiply this*x (a matrix-vector multiply) with options to transpose this
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Multiply(const HDSA::Ptr<HDSA::Vector<RealT> > & x, bool A_Trans = false) const 
  {
    int Am = numRows();
    int An = numCols();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Ax;
    if(A_Trans)
      {
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(Am,1);
	x_mat->Write_Vector_to_Column(0,x);
	Ax = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(An,1);
	Multiply(Ax, x_mat, A_Trans);
      }
    else
      {
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(An,1);
	x_mat->Write_Vector_to_Column(0,x);
	Ax = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(Am,1);
	Multiply(Ax, x_mat, A_Trans);
      }
    return Ax;
  }

  // Multiply this*x (a matrix-vector multiply) with options to transpose this
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Multiply(const HDSA::Vector<RealT> & x, bool A_Trans = false) const 
  {
    int Am = numRows();
    int An = numCols();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Ax;
    if(A_Trans)
      {
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(Am,1);
	x_mat->Write_Vector_to_Column(0,x);
	Ax = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(An,1);
	Multiply(Ax, x_mat, A_Trans);
      }
    else
      {
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(An,1);
	x_mat->Write_Vector_to_Column(0,x);
	Ax = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(Am,1);
	Multiply(Ax, x_mat, A_Trans);
      }
    return Ax;
  }

  // Write the jth column to an HDSA::Vector vec
  void Write_Column_to_Vector(int j, HDSA::Ptr<HDSA::Vector<RealT> > & vec)
  {
    Write_Column_to_Vector(j,*vec);
  }

  // Write the jth column to an HDSA::Vector vec
  void Write_Column_to_Vector(int j, HDSA::Vector<RealT> & vec)
  {
    if(vec.Get_enforce_zeros())
      {
	int dim = vec.Get_map_full_to_reduced().size();
	int nonzero_dim = vec.Get_map_reduced_to_full().size();
	for(int i = 0; i < nonzero_dim; i++)
	  {
	    vec.Replace_Element(vec.Get_map_reduced_to_full()[i],(*A_)(i,j));
	  }
	// Handles the case of a Joint_Vector
	int dim_diff = vec.dimension()-dim;
	if(dim_diff>0)
	  {
	    for(int i = 0; i < dim_diff; i++)
	      {
		vec.Replace_Element(dim+i,(*A_)(i+nonzero_dim,j));
	      }
	  }
      }
    else
      {
	for(int i = 0; i < vec.dimension(); i++)
	  {
	    vec.Replace_Element(i,(*A_)(i,j));
	  }
      }
  }
  
  // Write an HDSA::Vector vec to the jth column, use Nom option to specify enforce_z_zeros
  void Write_Vector_to_Column(int j, const HDSA::Ptr<HDSA::Vector<RealT> > & vec)
  {
    Write_Vector_to_Column(j,*vec);
  }

  // Write an HDSA::Vector vec to the jth column
  void Write_Vector_to_Column(int j, const HDSA::Vector<RealT> & vec)
  {
    if(vec.Get_enforce_zeros())
      {
    	int dim = vec.Get_map_full_to_reduced().size();
    	int nonzero_dim = vec.Get_map_reduced_to_full().size();
    	for(int i = 0; i < nonzero_dim; i++)
    	  {
    	    Replace_Element(i,j,vec(vec.Get_map_reduced_to_full()[i]));
    	  }
    	// Handles the case of a Joint_Vector
    	int dim_diff = vec.dimension()-dim;
    	if(dim_diff>0)
    	  {
    	    for(int i = 0; i < dim_diff; i++)
    	      {
    		Replace_Element(nonzero_dim+i,j,vec(dim+i));
    	      }
    	  }
      }
    else
      {
    	for(int i = 0; i < vec.dimension(); i++)
    	  {
    	    Replace_Element(i,j,vec(i));
    	  }
      }
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

  HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > Get_Teuchos_Matrix(void)
  {
    return A_;
  }

  void Copy_from_Teuchos_Matrix(const HDSA::Ptr<Teuchos::SerialDenseMatrix<int, RealT> > & A)
  {
    int m = A->numRows();
    int n = A->numCols();
    for(int i = 0; i < m; i++)
      {
  	for(int j = 0; j < n; j++)
  	  {
  	    Replace_Element(i,j,(*A)(i,j));
  	  }
      }
  }

};

}

#endif
