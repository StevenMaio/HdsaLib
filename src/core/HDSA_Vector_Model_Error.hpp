#ifndef HDSA_VECTOR_MODEL_ERROR_HPP
#define HDSA_VECTOR_MODEL_ERROR_HPP

namespace HDSA
{

template <class RealT>
class Vector_Model_Error : public Vector<RealT>{

  std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_;
  int m_;
  int n_;

public:

  Vector_Model_Error(int m, const HDSA::Ptr<HDSA::Vector<RealT> > & z): m_(m)
  {
    n_ = z->dimension();
    theta_.resize(m_);
    for(int i = 0; i < m_; i++)
      {
	theta_[i] = z->Clone();
      }
  }

  virtual ~Vector_Model_Error()
  { }

  void map_two_to_one(int & k, int i, int j) const
  {
    k = i*n_+j;
  }

  void  map_one_to_two(int & i, int & j, int k) const
  {
    i = k/n_;
    j = k - i*n_;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Get_Row(int i) const
  {
    return theta_[i];
  }

  void Set_Row(const HDSA::Vector<RealT> & vec, int i)
  {
    theta_[i]->set(vec);
  }

  // Access the (i,j) element
  RealT operator () (int i, int k) const 
  {
    return (*theta_[i])(k);
  }

  // Replace the kth element of the vector by val
  void Replace_Element(int i, int k, RealT val) 
  {
    theta_[i]->Replace_Element(k,val);
  }

  // Access the (i,j) element
  RealT operator () (int k) const 
  {
    int i = 0;
    int j = 0;
    map_one_to_two(i,j,k);
    return (*theta_[i])(j);
  }

  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val) 
  {
    int i = 0;
    int j = 0;
    map_one_to_two(i,j,k);
    theta_[i]->Replace_Element(j,val);
  }

  // Get the data on this processor
  std::vector<RealT> Get_Data_on_Processor(void) const
  {
    std::vector<RealT> data;
    return data;
  }

  // Get the indices on this processor
  std::vector<int> Get_Indices_on_Processor(void) const 
  {
    std::vector<int> indices;
    return indices;
  }

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Vector_Model_Error<RealT> >(m_,theta_[0]);
    return vec;
  }

  // add x to this
  void plus( const HDSA::Vector<RealT> & x )
  {
    const Vector_Model_Error<RealT> &ex = dynamic_cast<const Vector_Model_Error<RealT>&>(x);
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->plus(*ex.Get_Row(i));
      }
  }

  // scale this by val
  void scale( const RealT val ) 
  {
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->scale(val);
      }
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const 
  {
    const Vector_Model_Error<RealT> &ex = dynamic_cast<const Vector_Model_Error<RealT>&>(x);
    RealT val = 0.0;
    for(int i = 0; i < m_; i++)
      {
	val += theta_[i]->dot(*ex.Get_Row(i));
      }
    return val;
  }

  // compute the norm of this
  RealT norm(void) const 
  {
    RealT val = 0.0;
    for(int i = 0; i < m_; i++)
      {
	val += theta_[i]->dot(*theta_[i]);
      }
    val = std::sqrt(val);
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) 
  {
    const Vector_Model_Error<RealT> &ex = dynamic_cast<const Vector_Model_Error<RealT>&>(x);
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->axpy(alpha,*ex.Get_Row(i));
      }
  }
 
  // set this=0
  void zero(void) 
  {
    this->setScalar(0.0);
  }

  // set this= ith canonical basis vector
  void basis( const int i )
  {
    this->zero();
    this->Replace_Element(i,1.0);
  }
 
  // return vector dimension
  int dimension() const 
  {
    return m_*n_;
  }

  // set this=x
  void set( const HDSA::Vector<RealT> &x ) 
  {
    const Vector_Model_Error<RealT> &ex = dynamic_cast<const Vector_Model_Error<RealT>&>(x);
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->set(*ex.Get_Row(i));
      }
  }

  // set this=val elementwise
  void setScalar( const RealT val ) 
  {
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->setScalar(val);
      }
  }

  // set entries of this to random numbers in [l,u]
  void randomize( const RealT l = 0.0, const RealT u = 1.0 ) 
  {
    for(int i = 0; i < m_; i++)
      {
	theta_[i]->randomize(l,u);
      }
  }
 
};

}

#endif
