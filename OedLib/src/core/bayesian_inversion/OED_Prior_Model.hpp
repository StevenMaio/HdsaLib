//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_PRIOR_MODEL_HPP
#define OEDLIB_PRIOR_MODEL_HPP

#include "OED_Likelihood_Interface.hpp"
#include "../base/vectors/OED_Vector.hpp"

namespace OED
{
  template <class RealT>
  class Prior_Interface
  {
  public:
    virtual ~Prior_Interface() {}

    virtual void Prior_Precision_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) = 0;
    virtual void Prior_Covariance_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) = 0;
    virtual void Get_Prior_Mean(Vector<RealT> &z_out) = 0;
    virtual void Prior_Factor_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) = 0;
    virtual int Param_Dimension() = 0;

    // TODO: move this at some point to a potentially new class
    virtual void Mass_Matrix_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) = 0;
    virtual void Mass_Matrix_Inverse_Apply(Vector<RealT> &z_out, Vector<RealT> &z_in) = 0;
  };
}

#endif //OEDLIB_PRIOR_MODEL_HPP
