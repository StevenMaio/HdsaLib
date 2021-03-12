#ifndef HDSA_ROL_VECTOR_ADAPTER_HPP
#define HDSA_ROL_VECTOR_ADAPTER_HPP

#include <algorithm>
#include <cstdlib>

#include "ROL_Vector.hpp"

template <class RealT>
class ROL_Vector : public HDSA::Vector<RealT> {
  
protected:
  HDSA::Ptr<ROL::Vector<RealT> > rol_vec_;
  
public:
  
  ROL_Vector(void)
  { }
  
  ~ROL_Vector()
  { }
  
  // Access the kth element of the vector
  virtual RealT operator () (int k) const = 0;

  // Replace the kth element of the vector by val
  virtual void Replace_Element(int k, RealT val) = 0;

  // Get the data on this processor
  virtual std::vector<RealT> Get_Data_on_Processor(void) const = 0;

  // Clone the vector
  virtual HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const = 0;

  // add x to this
  void plus( const HDSA::Vector<RealT> & x )
  {
    const ROL_Vector<RealT> &ex = dynamic_cast<const ROL_Vector<RealT>&>(x);
    rol_vec_->plus(*ex.get_rol_vec());
  }

  // scale this by val
  void scale( const RealT val )
  {
    rol_vec_->scale(val);
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const ROL_Vector<RealT> &ex = dynamic_cast<const ROL_Vector<RealT>&>(x);
    return rol_vec_->dot(*ex.get_rol_vec());
  }

  // compute the norm of this
  RealT norm(void) const 
  {
    return rol_vec_->norm();
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const ROL_Vector<RealT> &ex = dynamic_cast<const ROL_Vector<RealT>&>(x);
    rol_vec_->axpy(alpha,*ex.get_rol_vec());
  }
 
  // set this=0
  void zero(void)
  {
    rol_vec_->zero();
  }
 
  // set this= ith canonical basis vector
  void basis( const int i )
  {
    rol_vec_->zero();
    Replace_Element(i,1.0);
  }

  // return vector dimension
  int dimension() const
  {
    return rol_vec_->dimension();
  }

  // set this=x
  void set( const HDSA::Vector<RealT> &x )
  {
    const ROL_Vector<RealT> &ex = dynamic_cast<const ROL_Vector<RealT>&>(x);
    rol_vec_->set(*ex.get_rol_vec());
  }

  // set this=val elementwise
  void setScalar( const RealT val )
  {
    rol_vec_->setScalar(val);
  }

  // set entries of this to random numbers in [l,u]
  void randomize( const RealT l = 0.0, const RealT u = 1.0 )
  {
    rol_vec_->randomize(l,u);
  }
  
  // Access underlying ROL::Vector
  const HDSA::Ptr<ROL::Vector<RealT> > get_rol_vec(void) const
  {
    return rol_vec_;
  }
  
};


#endif
