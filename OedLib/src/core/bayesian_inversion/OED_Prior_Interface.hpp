//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_PRIOR_INTERFACE_HPP
#define OEDLIB_PRIOR_INTERFACE_HPP

#include <memory>

#include "OED_Error_Model_Interface.hpp"
#include "OED_Vector.hpp"

#include "OED_Ptr.hpp"

namespace OED
{
  template <class RealT>
  class Prior_Interface
  {
  public:
    virtual ~Prior_Interface() {}

    virtual void Prior_Precision_Apply(Vector<RealT> &m_out, Vector<RealT> &m_in) = 0;
    virtual void Prior_Covariance_Apply(Vector<RealT> &m_out, Vector<RealT> &m_in) = 0;
    virtual void Get_Prior_Mean(Vector<RealT> &m_out) = 0;
    virtual void Prior_Covariance_Factor_Apply(Vector<RealT> &m_out, Vector<RealT> &m_in) = 0;
    virtual int Param_Dimension() = 0;

    // TODO: move this at some point to a potentially new class
    virtual void Mass_Matrix_Apply(Vector<RealT> &m_out, Vector<RealT> &m_in) = 0;
    virtual void Mass_Matrix_Inverse_Apply(Vector<RealT> &m_out, Vector<RealT> &m_in) = 0;

    virtual Ptr<Vector<RealT>> Sample_Vector() = 0;
  };
}

#endif //OEDLIB_PRIOR_INTERFACE_HPP
