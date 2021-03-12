#ifndef HDSA_ROL_INTERFACE_HPP
#define HDSA_ROL_INTERFACE_HPP

#include "ROL_Vector.hpp"

template <class RealT>
class HDSA_Rol_Vector : public ROL::Vector<RealT>
{

public:
  HDSA::Ptr<HDSA::Vector<RealT> > vec;

  HDSA_Rol_Vector(const HDSA::Ptr<HDSA::Vector<RealT> > & vec_in)
  {
    vec = vec_in->Clone();
    vec->set(*vec_in);
  }

  virtual ~HDSA_Rol_Vector() {}

  void plus( const ROL::Vector<RealT> & x )
  {
    HDSA_Rol_Vector* My_x;
    My_x = dynamic_cast<HDSA_Rol_Vector*>(&const_cast<ROL::Vector<RealT> &>(x));
    vec->plus(*My_x->vec);   
  }

  void scale( const RealT alpha )
  {
    vec->scale(alpha);
  }

  RealT dot( const ROL::Vector<RealT> & x ) const
  {
    HDSA_Rol_Vector* My_x;
    My_x = dynamic_cast<HDSA_Rol_Vector*>(&const_cast<ROL::Vector<RealT> &>(x));
    return vec->dot(*My_x->vec);   
  }

  RealT norm(void) const 
  {
    return vec->norm();
  }

  HDSA::Ptr<ROL::Vector<RealT> > clone() const
  {
    HDSA::Ptr<ROL::Vector<RealT> > clone_vec = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(vec);
    return clone_vec;
  }

  HDSA::Ptr<ROL::Vector<RealT> > basis( const int i ) const 
  {
    HDSA::Ptr<ROL::Vector<RealT> > clone_vec = HDSA::makePtr<HDSA_Rol_Vector<RealT> >(vec);
    clone_vec->basis(i);
    return clone_vec;  
  }

  int dimension() const
  {
    return vec->dimension();
  }

  void setScalar( const RealT C ) 
  {
    vec->setScalar(C);
  }

  void randomize( const RealT l = 0.0, const RealT u = 1.0 ) 
  {
    vec->randomize(l,u);
  }


};

#endif
