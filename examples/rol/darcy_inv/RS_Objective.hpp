#ifndef RS_OBJECTIVE_HPP
#define RS_OBJECTIVE_HPP

template<typename RealT>
class RS_Objective : public ROL::LinearCombinationObjective<RealT> {
private:
  std::vector<ROL::Ptr<ROL::Objective<RealT>>> obj_;

public:
  RS_Objective(const std::vector<RealT> &weights,const std::vector<ROL::Ptr<ROL::Objective<RealT>>> &obj): 
  ROL::LinearCombinationObjective<RealT>(weights,obj), obj_(obj)
  {}

  void precond( ROL::Vector<RealT> &Pv, const ROL::Vector<RealT> &v, const ROL::Vector<RealT> &x, RealT &tol ) override {
    obj_[1]->precond(Pv,v,x,tol);
    Pv.scale(1.e-3);
  }

};

#endif