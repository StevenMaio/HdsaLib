#ifndef HDSA_ANASAZI_INTERFACE_HPP
#define HDSA_ANASAZI_INTERFACE_HPP

#include "AnasaziMultiVec.hpp"
#include "AnasaziOperator.hpp"
#include <random>

template <class RealT>
class HDSA_Anasazi_Vector : public Anasazi::MultiVec<RealT>
{
private:
  const int NumVecs_;
  unsigned seed_;
  std::default_random_engine generator_;
  std::uniform_real_distribution<> distribution_;
  bool is_joint_;
  int theta_dim_;
  int nonzero_z_dim_;

public:
  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > vec;
  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > z_vec;
  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_vec;

  HDSA_Anasazi_Vector(const HDSA::Ptr<HDSA::Vector<RealT> > & vec_in, const int & NumVecs, const bool is_joint = false): 
    NumVecs_(NumVecs), is_joint_(is_joint)
  {
    if(!is_joint_)
      {
	vec.resize(NumVecs_);
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k] = vec_in->Clone();
	  }
      }
    else
      {
	vec.resize(NumVecs_);
	z_vec.resize(NumVecs_);
	theta_vec.resize(NumVecs_); 
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k] = vec_in->Clone();
	    HDSA::Joint_Vector<RealT>* eveck = dynamic_cast<HDSA::Joint_Vector<RealT>* >(&(*vec[k]));
	    z_vec[k]  = eveck->Get_Component_Vector_1();
	    theta_vec[k] =  eveck->Get_Component_Vector_2();
	  }
	theta_dim_ = theta_vec[0]->dimension();
	nonzero_z_dim_ = z_vec[0]->Get_nonzero_dim();
      }
    MvInit(0.0);

    seed_ = 395;
    generator_.seed(seed_);
    distribution_ = std::uniform_real_distribution<>(0.0,1.0);
  }

  HDSA_Anasazi_Vector(const HDSA::Ptr<HDSA::Vector<RealT> > & z_vec_in, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_vec_in, const int & NumVecs): 
    NumVecs_(NumVecs)
  {
    is_joint_ = true;
    theta_dim_ = theta_vec_in->dimension();
    nonzero_z_dim_ = z_vec_in->Get_nonzero_dim();

    vec.resize(NumVecs_);
    z_vec.resize(NumVecs_);
    theta_vec.resize(NumVecs_); 
    for(int k = 0; k < NumVecs_; k++)
      {
	vec[k] = HDSA::makePtr<HDSA::Joint_Vector<RealT> >(z_vec_in,theta_vec_in);
	HDSA::Joint_Vector<RealT>* eveck = dynamic_cast<HDSA::Joint_Vector<RealT>* >(&(*vec[k]));
	z_vec[k]  = eveck->Get_Component_Vector_1();
	theta_vec[k] =  eveck->Get_Component_Vector_2();
      }
    MvInit(0.0);
    
    seed_ = 395;
    generator_.seed(seed_);
    distribution_ = std::uniform_real_distribution<>(0.0,1.0);
  }

  //! Destructor
  ~HDSA_Anasazi_Vector()
  { }
  
  //! Returns a clone of the current vector.
  HDSA_Anasazi_Vector* Clone(const int NumberVecs) const
  {
    HDSA_Anasazi_Vector* tmp;
    if(!is_joint_)
      {
	tmp = new HDSA_Anasazi_Vector(vec[0], NumberVecs);
      }
    else
      {
	tmp = new HDSA_Anasazi_Vector(z_vec[0],theta_vec[0], NumberVecs);
      }
    return tmp;
  }
  
  // Returns a clone of the current multi-vector.
  HDSA_Anasazi_Vector* CloneCopy() const
  {
    HDSA_Anasazi_Vector* tmp;
    if(!is_joint_)
      {
	tmp = new HDSA_Anasazi_Vector(vec[0], NumVecs_);
	for(int k = 0; k < NumVecs_; k++)
	  {
	    tmp->vec[k]->set(*vec[k]);
	  }
      }
    else
      {
	tmp = new HDSA_Anasazi_Vector(z_vec[0],theta_vec[0], NumVecs_);
	for(int k = 0; k < NumVecs_; k++)
	  {
	    tmp->z_vec[k]->set(*z_vec[k]);
	    tmp->theta_vec[k]->set(*theta_vec[k]);
	  }
      }
    return tmp;
  }
  
  //! Returns a clone copy of specified vectors.
  HDSA_Anasazi_Vector* CloneCopy(const std::vector< int > &index) const
  {
    int size = index.size();
    HDSA_Anasazi_Vector* tmp;
    if(!is_joint_)
      {
	tmp = new HDSA_Anasazi_Vector(vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->vec[k]->set(*vec[index[k]]);
	  }
      }
    else
      {
	tmp = new HDSA_Anasazi_Vector(z_vec[0],theta_vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->z_vec[k]->set(*z_vec[index[k]]);
	    tmp->theta_vec[k]->set(*theta_vec[index[k]]);
	  }
      }
    return tmp;
  }
  
  //! Returns a view of current vector (shallow copy)
  HDSA_Anasazi_Vector* CloneViewNonConst(const std::vector< int > &index) 
  {
    int size = index.size();
    HDSA_Anasazi_Vector* tmp;
    if(!is_joint_)
      {
	tmp = new HDSA_Anasazi_Vector(vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->vec[k] = vec[index[k]];
	  }
      }
    else
      {
	tmp = new HDSA_Anasazi_Vector(z_vec[0],theta_vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->vec[k] = vec[index[k]];
	    tmp->z_vec[k] = z_vec[index[k]];
	    tmp->theta_vec[k] = theta_vec[index[k]];
	  }
      }
    return tmp;
  }
  
  //! Returns a view of current vector (shallow copy), const version.
  const HDSA_Anasazi_Vector* CloneView(const std::vector< int > &index) const
  {
    int size = index.size();
    HDSA_Anasazi_Vector* tmp;
    if(!is_joint_)
      {
	tmp = new HDSA_Anasazi_Vector(vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->vec[k] = vec[index[k]];
	  }
      }
    else
      {
	tmp = new HDSA_Anasazi_Vector(z_vec[0],theta_vec[0], size);
	for(int k = 0; k < size; k++)
	  {
	    tmp->vec[k] = vec[index[k]];
	    tmp->z_vec[k] = z_vec[index[k]];
	    tmp->theta_vec[k] = theta_vec[index[k]];
	  }
      }
    return tmp;
  }
  
  ptrdiff_t GetGlobalLength () const
  {
    ptrdiff_t k;
    if(!is_joint_)
      {
	k = vec[0]->dimension();
      }
    else
      {
	k = nonzero_z_dim_ + theta_dim_;
      }
    return k;
  }
  
  int GetNumberVecs () const
  {
    return NumVecs_;
  }

 void MvTimesMatAddMv (RealT alpha, const Anasazi::MultiVec<RealT> &A, 
                        const Teuchos::SerialDenseMatrix<int, RealT> &B, 
                        RealT beta)
  {
    assert (GetGlobalLength() == A.GetGlobalLength());
    assert (B.numRows() == A.GetNumberVecs());
    assert (B.numCols() == NumVecs_);

    if(!is_joint_)
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	HDSA_Anasazi_Vector<RealT> Mytmp(vec[0],NumVecs_);
	
	int dim = vec[0]->dimension();
	RealT val = 0;

	for(int i = 0; i < dim; i++)
	  {
	    for(int j = 0; j < NumVecs_; j++)
	      {
		// Populate ith entry of Mytmp.vec[j] by dot product of ith entries of MyA opt vectors with jth column of B
		for(int k = 0; k < MyA->NumVecs_; k++)
		  {
		    val += (*MyA->vec[k])(i)*B(k,j);
		  }
		Mytmp.vec[j]->Replace_Element(i,val);
		val = 0;
	      }
	  }
	
	// Take linear combination of Mytmp and *this
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k]->scale(beta);
	    vec[k]->axpy(alpha,*(Mytmp.vec[k]));
	  }
      }
    else
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	HDSA_Anasazi_Vector<RealT> Mytmp(z_vec[0], theta_vec[0], NumVecs_);

	RealT val = 0;
    
	for(int i = 0; i < nonzero_z_dim_; i++)
	  {
	    for(int j = 0; j < NumVecs_; j++)
	      {
		// Populate ith entry of Mytmp.z_vec[j] by dot product of ith entries of MyA opt vectors with jth column of B
		for(int k = 0; k < MyA->NumVecs_; k++)
		  {
		    val += (*MyA->z_vec[k])(MyA->z_vec[k]->Get_map_reduced_to_full(i))*B(k,j);
		  }
		Mytmp.z_vec[j]->Replace_Element(Mytmp.z_vec[j]->Get_map_reduced_to_full(i),val);
		val = 0;
	      }
	  }
	
	for(int i = nonzero_z_dim_; i < nonzero_z_dim_ + theta_dim_; i++)
	  {
	    for(int j = 0; j < NumVecs_; j++)
	      {
		// Populate (i-nonzero_z_dim_)th entry of Mytmp.theta_vec[j] by dot product of (i-nonzero_z_dim_)th entries of MyA parameter vectors with jth column of B
		for(int k = 0; k < MyA->NumVecs_; k++)
		  {
		    val += (*MyA->theta_vec[k])(i-nonzero_z_dim_)*B(k,j);
		  }
		Mytmp.theta_vec[j]->Replace_Element(i-nonzero_z_dim_,val);
		val = 0;
	      }
	  }
	
	// Take linear combination of Mytmp and *this
	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->scale(beta);
	    z_vec[k]->axpy(alpha,*Mytmp.z_vec[k]);
	    theta_vec[k]->scale(beta);
	    theta_vec[k]->axpy(alpha,*Mytmp.theta_vec[k]);
	  }
      }
    
  }


  // Replace *this with alpha * A + beta * B. 
  void MvAddMv (RealT alpha, const Anasazi::MultiVec<RealT>& A, 
                RealT beta,  const Anasazi::MultiVec<RealT>& B)
  {
    assert (NumVecs_ == A.GetNumberVecs());
    assert (NumVecs_ == B.GetNumberVecs()); 
    assert (GetGlobalLength() == A.GetGlobalLength());
    assert (GetGlobalLength() == B.GetGlobalLength());

    if(!is_joint_)
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	HDSA_Anasazi_Vector* MyB;
	MyB = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(B)); 
	
	HDSA_Anasazi_Vector<RealT> Mytmp(vec[0],NumVecs_);

	for(int k = 0; k < NumVecs_; k++)
	  {
	    Mytmp.vec[k]->set(*(MyA->vec[k]));
	    Mytmp.vec[k]->scale(alpha);
	    Mytmp.vec[k]->axpy(beta,*(MyB->vec[k]));
	    vec[k]->set(*Mytmp.vec[k]);
	  }
      }
    else
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	HDSA_Anasazi_Vector* MyB;
	MyB = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(B)); 

	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->set(*MyA->z_vec[k]);
	    z_vec[k]->scale(alpha);
	    z_vec[k]->axpy(beta,*MyB->z_vec[k]);
	    theta_vec[k]->set(*MyA->theta_vec[k]);
	    theta_vec[k]->scale(alpha);
	    theta_vec[k]->axpy(beta,*MyB->theta_vec[k]);
	  }
      }
  }
  
  // Compute a dense matrix B through the matrix-matrix multiply alpha * A^H * (*this). 
  void MvTransMv (RealT alpha, const Anasazi::MultiVec<RealT>& A, 
                  Teuchos::SerialDenseMatrix< int, RealT >& B) const
  {
    assert (A.GetGlobalLength() == GetGlobalLength());
    assert (NumVecs_ == B.numCols());
    assert (A.GetNumberVecs() == B.numRows());

    if(!is_joint_)
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(int i = 0; i < MyA->NumVecs_; i++)
	  {
	    for(int j = 0; j < NumVecs_; j++)
	      {
		B(i,j) = alpha*(MyA->vec[i]->dot(*vec[j]));
	      }
	  }
      }
    else
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(int i = 0; i < MyA->NumVecs_; i++)
	  {
	    for(int j = 0; j < NumVecs_; j++)
	      {
		B(i,j) = alpha*( MyA->z_vec[i]->dot(*z_vec[j]) + MyA->theta_vec[i]->dot(*theta_vec[j]) );
	      }
	  }
      }
  }
  
  
  // Compute a vector b where the components are the individual dot-products, i.e.b[i] = A[i]^H*this[i] where A[i] is the i-th column of A. 
  void MvDot (const Anasazi::MultiVec<RealT>& A, std::vector<RealT> &b) const
  {
    assert (NumVecs_ <= (int)b.size());
    assert (NumVecs_ == A.GetNumberVecs());
    assert (GetGlobalLength() == A.GetGlobalLength());

    if(!is_joint_)
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(int k = 0; k < NumVecs_; k++)
	  {
	    b[k] = MyA->vec[k]->dot(*vec[k]);
	  }
      }
    else
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(int k = 0; k < NumVecs_; k++)
	  {
	    b[k] = MyA->z_vec[k]->dot(*z_vec[k]) + MyA->theta_vec[k]->dot(*theta_vec[k]);
	  }
      }
  }  

  // Scale the vectors by alpha
  void MvScale( RealT alpha )
  {
    if(!is_joint_)
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k]->scale(alpha);
	  }
      }
    else
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->scale(alpha);
	    theta_vec[k]->scale(alpha);
	  }
      }
  }

  // Scale the i-th vector by alpha[i]
  void MvScale( const std::vector<RealT>& alpha )
  {
    if(!is_joint_)
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k]->scale(alpha[k]);
	  }
      }
    else
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->scale(alpha[k]);
	    theta_vec[k]->scale(alpha[k]);
	  }
      }
  }
  
  void MvNorm ( std::vector<typename Teuchos::ScalarTraits<RealT>::magnitudeType> &normvec ) const 
  {
    if(!is_joint_)
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    normvec[k] = vec[k]->norm();
	  }
      }
    else
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    RealT val = z_vec[k]->dot(*z_vec[k]);
	    val += theta_vec[k]->dot(*theta_vec[k]);
	    normvec[k] = std::sqrt(val);
	  }
      }
  }
  
  // Copy the vectors in A to a set of vectors in *this. The numvecs vectors in 
  // A are copied to a subset of vectors in *this indicated by the indices given 
  // in index.
  void SetBlock (const Anasazi::MultiVec<RealT>& A, 
                 const std::vector<int> &index)
  {
    if(!is_joint_)
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(unsigned int k = 0; k < index.size(); k++)
	  {
	    vec[index[k]]->set(*MyA->vec[k]);
	  }
      }
    else
      {
	HDSA_Anasazi_Vector* MyA;
	MyA = dynamic_cast<HDSA_Anasazi_Vector*>(&const_cast<Anasazi::MultiVec<RealT> &>(A)); 
	for(unsigned int k = 0; k < index.size(); k++)
	  {
	    z_vec[index[k]]->set(*MyA->z_vec[k]);
	    theta_vec[index[k]]->set(*MyA->theta_vec[k]);
	  }
      }
  }
  
  // Replace each element of the vectors in *this with alpha.
  void  MvInit (RealT alpha) {
    if(!is_joint_)
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k]->zero();
	    if(vec[k]->Get_enforce_zeros())
	      {
		for(unsigned int i = 0; i < vec[k]->Get_map_reduced_to_full().size(); i++)
		  {
		    vec[k]->Replace_Element(vec[k]->Get_map_reduced_to_full()[i],alpha);
		  }
	      }
	    else
	      {
		vec[k]->setScalar(alpha);
	      }
	  }
      }
    else
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->zero();
	    for(int i = 0; i < nonzero_z_dim_; i++)
	      {
		z_vec[k]->Replace_Element(z_vec[k]->Get_map_reduced_to_full(i),alpha);
	      }
	    theta_vec[k]->setScalar(alpha);
	  }
      }
  }

  // Fill the vectors in *this with random numbers.
  void  MvRandom ()
  {
    if(!is_joint_)
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    vec[k]->zero();
	    if(vec[k]->Get_enforce_zeros())
	      {
		for(unsigned int i = 0; i < vec[k]->Get_map_reduced_to_full().size(); i++)
		  {
		    vec[k]->Replace_Element(vec[k]->Get_map_reduced_to_full()[i],distribution_(generator_));
		  }
	      }
	    else
	      {
		for(int i = 0; i < vec[k]->dimension(); i++)
		  {
		    vec[k]->Replace_Element(i,distribution_(generator_));
		  }
	      }
	  }
      }
    else
      {
	for(int k = 0; k < NumVecs_; k++)
	  {
	    z_vec[k]->zero();
	    for(int i = 0; i < nonzero_z_dim_; i++)
	      {
		z_vec[k]->Replace_Element(z_vec[k]->Get_map_reduced_to_full(i),distribution_(generator_));
	      }
	    
	    for(int i = 0; i < theta_dim_; i++)
	      {
		theta_vec[k]->Replace_Element(i, distribution_(generator_));
	      }
	  }
      }
  }
  

  void MvPrint (std::ostream &os) const
  {
     os << "Object HDSA_Anasazi_Vector" << std::endl;
     os << "Number of rows = " << GetGlobalLength() << std::endl;
     os << "Number of vecs = " << NumVecs_ << std::endl;
  }
  
};

template <class RealT>
class HDSA_Anasazi_Operator : public Anasazi::Operator<RealT>
{
private:
  HDSA::Ptr<HDSA::Linear_Operator<RealT> > A_;
  
public:
  HDSA_Anasazi_Operator(const HDSA::Ptr<HDSA::Linear_Operator<RealT> > & A): A_(A)
  { }
  
  ~HDSA_Anasazi_Operator()
  {}
  
  void Apply(const Anasazi::MultiVec<RealT> & x, Anasazi::MultiVec<RealT> & y) const
  {
    const HDSA_Anasazi_Vector<RealT>* MyX;
    MyX = dynamic_cast<const HDSA_Anasazi_Vector<RealT>* >(&x); 
    HDSA_Anasazi_Vector<RealT>* MyY;
    MyY = dynamic_cast<HDSA_Anasazi_Vector<RealT>* >(&y);	
    for(int k = 0; k < MyX->GetNumberVecs(); k++)
      {
	A_->matvec(MyY->vec[k],MyX->vec[k]);
      }
  }
};

#endif 

