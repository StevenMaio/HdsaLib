//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LIKELIHOOD_MODEL_HPP
#define OEDLIB_LIKELIHOOD_MODEL_HPP

#include "OED_Vector.hpp"

using OED::Vector;

namespace OED
{
  template<class RealT>
  class Error_Model_Interface
  {
  public:
    virtual ~Error_Model_Interface()
    {
    }

    virtual void Noise_Precision_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) = 0;

    virtual void Noise_Covariance_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) = 0;

    virtual int Data_Dimension() = 0;

    // TODO: implement other methods
  };
}

#endif //OEDLIB_LIKELIHOOD_MODEL_HPP
