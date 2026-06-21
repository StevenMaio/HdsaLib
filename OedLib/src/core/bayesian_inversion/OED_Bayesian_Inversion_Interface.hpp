//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_BAYESIAN_INVERSION_HPP
#define OEDLIB_BAYESIAN_INVERSION_HPP

#include <memory>

#include "OED_Error_Model_Interface.hpp"
#include "OED_Prior_Interface.hpp"
#include "OED_Model_Interface.hpp"
#include "OED_Observation_Operator_Interface.hpp"

#include "OED_Ptr.hpp"

namespace OED
{
  template <class RealT>
  class Bayesian_Inversion_Interface
  {
  private:
    Ptr<Model_Interface<RealT>> model_;
    Ptr<Observation_Operator_Interface<RealT>> obs_operator_;
    Ptr<Prior_Interface<RealT>> prior_;
    Ptr<Error_Model_Interface<RealT>> error_model_;

  public:

    Bayesian_Inversion_Interface(
        Ptr<Model_Interface<RealT>> model,
        Ptr<Observation_Operator_Interface<RealT>> obs_operator,
        Ptr<Prior_Interface<RealT>> prior,
        Ptr<Error_Model_Interface<RealT>> error_model
    )
      : model_(model), obs_operator_(obs_operator),
        prior_(prior), error_model_(error_model) {}

    Ptr<Error_Model_Interface<RealT>> &Error_Model()
    {
      return this->error_model_;
    }

    Ptr<Prior_Interface<RealT>> &Prior()
    {
      return this->prior_;
    }

    Ptr<Model_Interface<RealT>> &Model()
    {
      return this->model_;
    }

    Ptr<Observation_Operator_Interface<RealT>> &Observation_Operator()
    {
      return this->obs_operator_;
    }

    virtual ~Bayesian_Inversion_Interface() {}

    // Need these to build the relevant things
    // TODO: change these to use smart pointers later
    virtual Ptr<Vector<RealT>> Get_Empty_Parameter_Vector() = 0;
    virtual Ptr<Vector<RealT>> Get_Empty_State_Vector() = 0;
    virtual Ptr<Vector<RealT>> Get_Empty_Data_Vector() = 0;    // TODO: need to eventually think about how to deal with reducing size of data vector
  };
}

#endif //OEDLIB_BAYESIAN_INVERSION_HPP
