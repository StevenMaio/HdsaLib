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
  std::normal_distribution<RealT> distribution_;

public:  
  Std_Vector(int dim): dim_(dim)
  {
    vec_ = HDSA::makePtr<std::vector<RealT> >(dim,0.0);
    seed_ = std::rand();
    generator_.seed(seed_);
    distribution_ = std::normal_distribution<RealT>(0.0,1.0); 
  }
  
  ~Std_Vector()
  { }

  //////////////////////////////////////////////////////////////////////////////////
  // Overloading pure virtual functions in HDSA::Vector base class
  //////////////////////////////////////////////////////////////////////////////////
  
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Std_Vector<RealT> >(dim_);
    return vec;
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    RealT val = 0.0;
    const Std_Vector<RealT> x_std = dynamic_cast<const Std_Vector<RealT>&>(x);
    for(int k = 0; k < dim_; k++)
      {
	val += (*vec_)[k]*(x_std(k));
      }
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const Std_Vector<RealT> x_std = dynamic_cast<const Std_Vector<RealT>&>(x);
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] += alpha*(x_std(k));
      }
  }

  // return vector dimension
  int dimension() const
  {
    return dim_;
  }

  // set this=val elementwise
  void setScalar( const RealT val )
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = val;
      }
  }

  void randomize_standard_normal( ) 
  {
    for(int k = 0; k < dim_; k++)
      {
	(*vec_)[k] = distribution_(generator_);
      }
  }

  //////////////////////////////////////////////////////////////////////////////////
  // Function specific to this class for convenience
  //////////////////////////////////////////////////////////////////////////////////

  // Access underlying std::vector
  const HDSA::Ptr<std::vector<RealT> > get_std_vec(void) const
  {
    return vec_;
  }
  
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

};


#endif
