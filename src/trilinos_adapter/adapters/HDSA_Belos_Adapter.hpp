#ifndef HDSA_BELOS_ADAPTER_HPP
#define HDSA_BELOS_ADAPTER_HPP

#include "BelosMultiVec.hpp"
#include "BelosOperator.hpp"

template <class RealT>
class HDSA_Belos_Vector : public Belos::MultiVec<RealT>
{
private:
  const int NumVecs_;

public:
  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > vec;

  HDSA_Belos_Vector(const HDSA::Vector<RealT> & vec_in, const int & NumVecs): 
    NumVecs_(NumVecs)
  {
    vec.resize(NumVecs_);
    for(int k = 0; k < NumVecs_; k++)
      {
	// Create vector and set to zero
	vec[k] = vec_in.clone();
      }
  }

  HDSA_Belos_Vector(const int & NumVecs): 
    NumVecs_(NumVecs)
  {
    vec.resize(NumVecs_);
  }
  
  //! Destructor
  ~HDSA_Belos_Vector()
  {}
  
  //! Returns a clone of the current vector.
  HDSA_Belos_Vector* Clone(const int NumberVecs) const
  {
    HDSA_Belos_Vector* tmp = new HDSA_Belos_Vector(*vec[0], NumberVecs);
    return tmp;
  }
  
  // Returns a clone of the current multi-vector.
  HDSA_Belos_Vector* CloneCopy() const
  {
    HDSA_Belos_Vector* tmp = new HDSA_Belos_Vector(*vec[0], NumVecs_);
    for(int k = 0; k < NumVecs_; k++)
      {
	tmp->vec[k]->set(*vec[k]);
      }
    return tmp;
  }
  
  //! Returns a clone copy of specified vectors.
  HDSA_Belos_Vector* CloneCopy(const std::vector< int > &index) const
  {
    int size = index.size();
    HDSA_Belos_Vector* tmp = new HDSA_Belos_Vector(*vec[0], size);
    for(int k = 0; k < size; k++)
      {
	tmp->vec[k]->set(*vec[index[k]]);
      }
    return tmp;
  }
 
  //! Returns a view of current vector (shallow copy)
  HDSA_Belos_Vector* CloneViewNonConst(const std::vector< int > &index) 
  {
    int size = index.size();
    HDSA_Belos_Vector* tmp = new HDSA_Belos_Vector(size);
    for(int i = 0; i < size; i++)
      {
  	tmp->vec[i]=vec[index[i]];
      }
    return tmp;
  }
  
  //! Returns a view of current vector (shallow copy), const version.
  const HDSA_Belos_Vector* CloneView(const std::vector< int > &index) const
  {
    int size = index.size();
    HDSA_Belos_Vector* tmp = new HDSA_Belos_Vector(size);
    for(int i = 0; i < size; i++)
      {
  	tmp->vec[i]=vec[index[i]];
      }
    return tmp;
  }
  
  ptrdiff_t GetGlobalLength () const
  {
    return vec[0]->dimension();
  }
  
  int GetNumberVecs () const
  {
    return NumVecs_;
  }

  // Update *this with alpha * A * B + beta * (*this). 
  void MvTimesMatAddMv (RealT alpha, const Belos::MultiVec<RealT> &A, 
                        const Teuchos::SerialDenseMatrix<int, RealT> &B, 
                        RealT beta)
  {
    assert (GetGlobalLength() == A.GetGlobalLength());
    assert (B.numRows() == A.GetNumberVecs());
    assert (B.numCols() == NumVecs_);

    HDSA_Belos_Vector* MyA;
    MyA = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(A)); 

    for(int k = 0; k < NumVecs_; k++)
      {
	vec[k]->scale(beta);
	for(int i = 0; i < B.numRows(); i++)
	  {
	    vec[k]->axpy(alpha*B(i,k),*MyA->vec[i]);
	  }
      }

  }

 // Replace *this with alpha * A + beta * B. 
  void MvAddMv (RealT alpha, const Belos::MultiVec<RealT>& A, 
                RealT beta,  const Belos::MultiVec<RealT>& B)
  {
    assert (NumVecs_ == A.GetNumberVecs());
    assert (NumVecs_ == B.GetNumberVecs()); 
    assert (GetGlobalLength() == A.GetGlobalLength());
    assert (GetGlobalLength() == B.GetGlobalLength());

    HDSA_Belos_Vector* MyA;
    MyA = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(A)); 
    HDSA_Belos_Vector* MyB;
    MyB = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(B)); 

    HDSA_Belos_Vector<RealT> Mytmp(*vec[0],NumVecs_);

    for(int k = 0; k < NumVecs_; k++)
      {
	Mytmp.vec[k]->set(*(MyA->vec[k]));
	Mytmp.vec[k]->scale(alpha);
	Mytmp.vec[k]->axpy(beta,*(MyB->vec[k]));
	vec[k]->set(*Mytmp.vec[k]);
      }
  }
  
  // Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this). 
  void MvTransMv (RealT alpha, const Belos::MultiVec<RealT>& A, 
                  Teuchos::SerialDenseMatrix< int, RealT >& B) const
  {
    assert (A.GetGlobalLength() == GetGlobalLength());
    assert (NumVecs_ == B.numCols());
    assert (A.GetNumberVecs() == B.numRows());

    HDSA_Belos_Vector* MyA;
    MyA = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(A)); 
    for(int i = 0; i < MyA->NumVecs_; i++)
      {
	for(int j = 0; j < NumVecs_; j++)
	  {
	    B(i,j) = alpha*(MyA->vec[i]->dot(*vec[j]));
	  }
      }
  }
  
  // Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
  void MvDot (const Belos::MultiVec<RealT>& A, std::vector<RealT> &b) const
  {
    assert (NumVecs_ == A.GetNumberVecs());

    HDSA_Belos_Vector* MyA;
    MyA = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(A)); 
    for(int k = 0; k < NumVecs_; k++)
      {
	b[k] = MyA->vec[k]->dot(*vec[k]);
      }
  }  

  // Scale the vectors by alpha
  void MvScale( RealT alpha )
  {
    for(int k = 0; k < NumVecs_; k++)
      {
	vec[k]->scale(alpha);
      }
  }

  // Scale the i-th vector by alpha[i]
  void MvScale( const std::vector<RealT>& alpha )
  {
    for(int k = 0; k < NumVecs_; k++)
      {
	vec[k]->scale(alpha[k]);
      }
  }
  
  void MvNorm ( std::vector<typename Teuchos::ScalarTraits<RealT>::magnitudeType> &normvec, Belos::NormType type = Belos::TwoNorm) const 
  {
    for(int k = 0; k < NumVecs_; k++)
      {
	normvec[k] = std::sqrt(vec[k]->dot(*vec[k]));
      }
  }
  
  // Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
  // A are copied to a subset of vectors in *this indicated by the indices given 
  // in index.
  void SetBlock (const Belos::MultiVec<RealT>& A, 
                 const std::vector<int> &index)
  {
    HDSA_Belos_Vector* MyA;
    MyA = dynamic_cast<HDSA_Belos_Vector*>(&const_cast<Belos::MultiVec<RealT> &>(A)); 
       for(unsigned int k = 0; k < index.size(); k++)
	 {
	   vec[index[k]]->set(*(MyA->vec[k]));
	 }
  }
  
  // Replace each element of the vectors in *this with alpha.
  void  MvInit (RealT alpha)
  {
    for(int k = 0; k < NumVecs_; k++)
      {
        vec[k]->setScalar(alpha);
      }
  }

  // Fill the vectors in *this with random numbers.
  void  MvRandom ()
  {
    for(int k = 0; k < NumVecs_; k++)
      {
	vec[k]->randomize_standard_normal();
      }
  }
  
  void MvPrint (std::ostream &os) const
  {
     os << "Object HDSA_Belos_Vector" << std::endl;
     os << "Number of rows = " << GetGlobalLength() << std::endl;
     os << "Number of vecs = " << NumVecs_ << std::endl;
  }
  
};

// Overload Belos Operator to take matrix vector products for the linear system solve
template <class RealT>
class HDSA_Belos_Operator : public Belos::Operator<RealT>
{
  const HDSA::Linear_Operator<RealT>* A_;
  
public:
  HDSA_Belos_Operator(const HDSA::Linear_Operator<RealT>* A): A_(A)
  { }
  
  //! Dtor
  ~HDSA_Belos_Operator()
  {}
  
  void Apply(const Belos::MultiVec<RealT>& X, 
	     Belos::MultiVec<RealT>& Y,
	     Belos::ETrans trans=Belos::NOTRANS) const
  {
    const HDSA_Belos_Vector<RealT>* MyX;
    MyX = dynamic_cast<const HDSA_Belos_Vector<RealT>*>(&X);    
    HDSA_Belos_Vector<RealT>* MyY;
    MyY = dynamic_cast<HDSA_Belos_Vector<RealT>*>(&Y); 
    
    for(int k = 0; k < MyX->GetNumberVecs(); k++)
      {
	A_->matvec(*MyY->vec[k],*MyX->vec[k]);
      }
    
  }
  
};

#endif

