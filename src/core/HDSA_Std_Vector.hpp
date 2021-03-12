#ifndef HDSA_STDVECTOR_HPP
#define HDSA_STDVECTOR_HPP

#include <algorithm>
#include <cstdlib>
#include <random>

template <class RealT>
class Std_Vector : public HDSA::Vector<RealT> {
  
private:
  int dim_;
  HDSA::Ptr<std::vector<RealT> > vec_;
  unsigned seed_;
  std::default_random_engine generator_;
  std::uniform_real_distribution<RealT> distribution_;
  
public:
  
  Std_Vector(int dim): dim_(dim)
  {
    vec_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
    seed_ = 832;
    generator_.seed(seed_);
    distribution_ = std::uniform_real_distribution<RealT>(0.0,1.0);
  }
  
  ~Std_Vector()
  { }
  
  // Access the (i,j) element
  RealT operator () (int k) const 
  {
    return (*vec_)[k];
  }

  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val)
  {
    (*vec_)[k] = val;
  } 

  // Get the data on this processor
  std::vector<RealT> Get_Data_on_Processor(void) const
  {
    return *vec_;
  }

  // Get the indices on this processor
  std::vector<int> Get_Indices_on_Processor(void) const 
  {
    std::vector<int> indices = std::vector<int>(vec_->size());
    for(unsigned int k = 0; k < vec_->size(); k++)
      {
	indices[k] = k;
      }
    return indices;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Std_Vector<RealT> >(dim_);
    return vec;
  }

  // add x to this
  void plus( const HDSA::Vector<RealT> & x )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] += x(k);
      }
  }

  // scale this by val
  void scale( const RealT val )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] *= val;
      }
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    RealT val = 0.0;
    for(int k = 0; k < dim_; k++)
      {
	val += (*vec_)[k]*(x(k));
      }
    return val;
  }

  // compute the norm of this
  RealT norm(void) const 
  {
    return std::sqrt(this->dot(*this));
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] += alpha*(x(k));
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
    return dim_;
  }

  // set this=x
  void set( const HDSA::Vector<RealT> &x )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = x(k);
      }
  }

  // set this=val elementwise
  void setScalar( const RealT val )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = val;
      }
  }

  // set entries of this to random numbers in [l,u]
  void randomize( const RealT l = 0.0, const RealT u = 1.0 )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = distribution_(generator_)*(u-l) + l;
      }
  }

  // Access underlying std::vector
  const HDSA::Ptr<std::vector<RealT> > get_std_vec(void) const
  {
    return vec_;
  }
  
};


#endif
