#ifndef HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_HPP
#define HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_HPP

#include "HDSA_PC_Auxillary_Parameter_Trajectory.hpp"

namespace HDSA
{

  template <class RealT>
  class PC_Sensitivity_Operator_Interface
  {

  public:
    PC_Sensitivity_Operator_Interface()
    {
    }

    virtual ~PC_Sensitivity_Operator_Interface()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual void Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::PC_Auxillary_Parameter_Trajectory<RealT> &theta_traj, RealT &time_index) const = 0;

    virtual void Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::PC_Auxillary_Parameter_Trajectory<RealT> &theta_traj, RealT &time_index) const = 0;

    virtual void Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z, const HDSA::PC_Auxillary_Parameter_Trajectory<RealT> &theta_traj, RealT &time_index) const = 0;
 
  };

}

#endif
