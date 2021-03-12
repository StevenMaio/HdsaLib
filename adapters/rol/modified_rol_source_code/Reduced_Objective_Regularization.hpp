#ifndef REDUCED_OBJECTIVE_REGULARIZATION_HPP
#define REDUCED_OBJECTIVE_REGULARIZATION_HPP

#include "ROL_Objective_SimOpt.hpp"

template <class Real>
class Reduced_Objective_Regularization : public ROL::Objective<Real> {
private:
  ROL::Ptr<ROL::Objective_SimOpt<Real> > obj_;
  ROL::Ptr<ROL::Vector<Real> > state_;


public:

  Reduced_Objective_Regularization(const ROL::Ptr<ROL::Objective_SimOpt<Real> > &obj, const ROL::Ptr<ROL::Vector<Real> > & state) 
    : obj_(obj)
  { 
    state_ = state->clone();
  }

  Real value( const ROL::Vector<Real> &z, Real &tol ) {
    return obj_->value(*state_,z,tol);
  }

  void gradient( ROL::Vector<Real> &g, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->gradient_2(g,*state_,z,tol);
  }

  void hessVec( ROL::Vector<Real> &hv, const ROL::Vector<Real> &v, const ROL::Vector<Real> &z, Real &tol ) {
    obj_->hessVec_22(hv,v,*state_,z,tol);
  }

}; // class Reduced_Objective_Regularization


#endif
