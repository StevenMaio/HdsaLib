#ifndef PDE_LTIOBJECTIVE_TS_HPP
#define PDE_LTIOBJECTIVE_TS_HPP

#include "ROL_DynamicObjective.hpp"

template <class Real>
class LTI_Objective_TS : public ROL::DynamicObjective<Real> {
private:
  const ROL::Ptr<Objective_SimOpt_TS<Real> > obj_;
  Real theta_, T_;
  mutable ROL::Ptr<ROL::Vector<Real> > zdual_, udual_;
  mutable bool isInit_;

  void initialize(const ROL::Vector<Real> &z, const ROL::Vector<Real> &u) const {
    if (!isInit_) {
      zdual_ = z.dual().clone();
      udual_ = u.dual().clone();
      isInit_ = true;
    }
  }

public:
  LTI_Objective_TS(ROL::ParameterList                    &parlist,
		   const ROL::Ptr<Objective_SimOpt_TS<Real> > & obj)
    : obj_ (obj),
      isInit_(false) {
    theta_ = parlist.sublist("Time Discretization").get("Theta",    1.0);
    T_     = parlist.sublist("Time Discretization").get("End Time", 1.0);
  }

  Real value( const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    // Integrated Objective
    Real valo(0), valn(0);
    obj_->update(uo,z,timeOld);
    valo = obj_->value(uo,z,tol);
    obj_->update(un,z,timeNew);
    valn = obj_->value(un,z,tol);
    
    return dt*((one-theta_)*valo + theta_*valn);
  }

  void gradient_uo( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    obj_->update(uo,z,timeOld);
    obj_->gradient_1(g,uo,z,tol);
    g.scale(dt*(one-theta_));
  }

  void gradient_un( ROL::Vector<Real> &g,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    g.zero();
    obj_->update(un,z,timeNew);
    obj_->gradient_1(*udual_,un,z,tol);
    g.axpy(dt*theta_, *udual_);
  }

  void gradient_z( ROL::Vector<Real> &g,
             const ROL::Vector<Real> &uo,
             const ROL::Vector<Real> &un,
             const ROL::Vector<Real> &z,
             const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    g.zero();
    obj_->update(uo,z,timeOld);
    obj_->gradient_2(*zdual_,uo,z,tol);
    g.axpy(dt*(one-theta_), *zdual_);
    obj_->update(un,z,timeNew);
    obj_->gradient_2(*zdual_,un,z,tol);
    g.axpy(dt*theta_,*zdual_);
  }

  void hessVec_uo_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z,timeOld);
    obj_->hessVec_11(*udual_,v,uo,z,tol);
    hv.axpy(dt*(one-theta_), *udual_);
  }

  void hessVec_uo_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    hv.zero();
  }

  void hessVec_uo_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z,timeOld);
    obj_->hessVec_12(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
  }

  void hessVec_un_uo( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    hv.zero();
  }

  void hessVec_un_un( ROL::Vector<Real> &hv,
                const ROL::Vector<Real> &v,
                const ROL::Vector<Real> &uo,
                const ROL::Vector<Real> &un,
                const ROL::Vector<Real> &z,
                const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    obj_->update(un,z,timeNew);
    obj_->hessVec_11(*udual_,v,un,z,tol);
    hv.axpy(dt*theta_, *udual_);
  }

  void hessVec_un_z( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    obj_->update(un,z,timeNew);
    obj_->hessVec_12(*udual_,v,un,z,tol);
    hv.axpy(dt*theta_, *udual_);
  }

  void hessVec_z_uo( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    obj_->update(uo,z,timeOld);
    obj_->hessVec_21(hv,v,uo,z,tol);
    hv.scale(dt*(one-theta_));
  }

  void hessVec_z_un( ROL::Vector<Real> &hv,
               const ROL::Vector<Real> &v,
               const ROL::Vector<Real> &uo,
               const ROL::Vector<Real> &un,
               const ROL::Vector<Real> &z,
               const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    obj_->update(un,z,timeNew);
    obj_->hessVec_21(*zdual_,v,un,z,tol);
    hv.axpy(dt*theta_,*zdual_);
  }

  void hessVec_z_z( ROL::Vector<Real> &hv,
              const ROL::Vector<Real> &v,
              const ROL::Vector<Real> &uo,
              const ROL::Vector<Real> &un,
              const ROL::Vector<Real> &z,
              const ROL::TimeStamp<Real> &ts ) const {
    initialize(z,un);
    const Real one(1);
    Real tol(std::sqrt(ROL::ROL_EPSILON<Real>()));
    Real timeOld = ts.t[0], timeNew = ts.t[1];
    Real dt = timeNew - timeOld;
    hv.zero();
    obj_->update(uo,z,timeOld);
    obj_->hessVec_22(*zdual_,v,uo,z,tol);
    hv.axpy(dt*(one-theta_),*zdual_);
    obj_->update(un,z,timeNew);
    obj_->hessVec_22(*zdual_,v,un,z,tol);
    hv.axpy(dt*theta_,*zdual_);
  }

public:
  void setParameter(const std::vector<Real> &param) {
    obj_->setParameter(param);
  }

}; // class LTI_Objective

#endif
