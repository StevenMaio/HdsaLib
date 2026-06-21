//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_TEST_LINEAR_BAYESIAN_INVERSION_HPP
#define OEDLIB_TEST_LINEAR_BAYESIAN_INVERSION_HPP

#include "OED_Std_Vector.hpp"

#include "OED_Model_Interface.hpp"
#include "OED_Error_Model_Interface.hpp"
#include "OED_Observation_Operator_Interface.hpp"
#include "OED_Prior_Interface.hpp"
#include "OED_Bayesian_Inversion_Interface.hpp"
#include "OED_Ptr.hpp"

using namespace OED;

namespace OED_TEST
{
  template <class RealT>
  class Test_Linear_Bayesian_Inversion : public Bayesian_Inversion_Interface<RealT>
  {
  private:
    Ptr<Vector<RealT>> data_;

  public:
    Test_Linear_Bayesian_Inversion(
        Ptr<Model_Interface<RealT>> model,
        Ptr<Observation_Operator_Interface<RealT>> obs_operator,
        Ptr<Prior_Interface<RealT>> prior,
        Ptr<Error_Model_Interface<RealT>> error_model
    )
      : Bayesian_Inversion_Interface<RealT>(model, obs_operator, prior, error_model)
    {
      this->data_ = this->Get_Empty_Data_Vector();
    }

    Ptr<Vector<RealT>> Get_Empty_Parameter_Vector() override
    {
      int param_dim = this->Prior()->Param_Dimension();
      Ptr<Std_Vector<RealT>> m = std::make_shared<Std_Vector<RealT>>(param_dim);
      return m;
    }

    Ptr<Vector<RealT>> Get_Empty_State_Vector() override
    {
      int state_dim = this->Model()->State_Dimension();
      Ptr<Std_Vector<RealT>> u = std::make_shared<Std_Vector<RealT>>(state_dim);
      return u;
    };

    Ptr<Vector<RealT>> Get_Empty_Data_Vector() override
    {
      int data_dim = this->Error_Model()->Data_Dimension();
      Ptr<Std_Vector<RealT>> d = std::make_shared<Std_Vector<RealT>>(data_dim);
      return d;
    };

    void Compute_MAP_Point(Ptr<Vector<double>> &m_out)
    {
      // TODO: need to build the linear equation I want to solve...
    }
  };
}

#endif // OEDLIB_TEST_LINEAR_BAYESIAN_INVERSION_HPP
