#ifndef HDSA_OPT_PROBLEM_OBJECTS_HPP
#define HDSA_OPT_PROBLEM_OBJECTS_HPP

// The user should derive off this class and define the members and pure virtual functions below

namespace HDSA
{

template <class RealT>
class Opt_Problem_Objects {

public:

  // z and theta must be implemented
  HDSA::Ptr<HDSA::Vector<RealT> > z;
  HDSA::Ptr<HDSA::Vector<RealT> > theta;

  // rs_obj must be implemented for a reduced space formulation
  HDSA::Ptr<HDSA::RS_Objective<RealT> > rs_obj;

  // if rs_obj is not implemented, then fs_obj, con, u, and lambda must be implemented for a full space formulation
  HDSA::Ptr<HDSA::FS_Objective<RealT> > fs_obj;
  HDSA::Ptr<HDSA::Constraint<RealT> > con;
  HDSA::Ptr<HDSA::Vector<RealT> > u;
  HDSA::Ptr<HDSA::Vector<RealT> > lambda;

  // Need to implement all objects with the exception of lambda for the model error HDSA

  Opt_Problem_Objects()
  { }

  virtual ~Opt_Problem_Objects()
  { }

  virtual HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) = 0;

  virtual void Solve_Optimization_Problem() = 0;

  virtual void Load_Optimal_Solution() = 0;

  virtual void Write_Optimal_Solution() = 0;
  
}; 

}

#endif
