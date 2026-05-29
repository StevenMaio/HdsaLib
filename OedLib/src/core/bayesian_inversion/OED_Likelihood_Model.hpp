//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LIKELIHOOD_MODEL_HPP
#define OEDLIB_LIKELIHOOD_MODEL_HPP

#include "../base/vectors/OED_Vector.hpp"

using OED::Vector;

namespace OED
{
  template<class RealT>
  class Likelihood_Model
  {
  public:
    virtual ~Likelihood_Model()
    {
    }

    virtual void Noise_Precision_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) = 0;

    virtual void Noise_Covariance_Apply(Vector<RealT> &d_out, Vector<RealT> &d_in) = 0;

    virtual void Observation_Operator_Apply(Vector<RealT> &d_out, Vector<RealT> &u_in) = 0;

    virtual void Observation_Operator_Transpose_Apply(Vector<RealT> &u_out, Vector<RealT> &d_in) = 0;

    virtual void Get_Observed_Data(Vector<RealT> &d) = 0;

    virtual int State_Dimension() = 0;

    virtual int Data_Dimension() = 0;

    // TODO: implement other methods
  };
}

#endif //OEDLIB_LIKELIHOOD_MODEL_HPP
