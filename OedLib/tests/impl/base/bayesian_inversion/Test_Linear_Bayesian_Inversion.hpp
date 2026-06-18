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

using namespace OED;

namespace OED_TEST
{
  template <class RealT>
  class Test_Linear_Bayesian_Inversion : public Bayesian_Inversion_Interface<RealT>
  {
  private:
    std::shared_ptr<Vector<RealT>> data_;

  public:
    Test_Linear_Bayesian_Inversion(
        std::shared_ptr<Model_Interface<RealT>> model,
        std::shared_ptr<Observation_Operator_Interface<RealT>> obs_operator,
        std::shared_ptr<Prior_Interface<RealT>> prior,
        std::shared_ptr<Error_Model_Interface<RealT>> error_model
    )
      : Bayesian_Inversion_Interface<RealT>(model, obs_operator, prior, error_model)
    {
      this->data_ = this->Get_Empty_Data_Vector();
    }

    std::shared_ptr<Vector<RealT>> Get_Empty_Parameter_Vector() override
    {
      int param_dim = this->Prior()->Param_Dimension();
      std::shared_ptr<Std_Vector<RealT>> m = std::make_shared<Std_Vector<RealT>>(param_dim);
      return m;
    }

    std::shared_ptr<Vector<RealT>> Get_Empty_State_Vector() override
    {
      int state_dim = this->Model()->State_Dimension();
      std::shared_ptr<Std_Vector<RealT>> u = std::make_shared<Std_Vector<RealT>>(state_dim);
      return u;
    };

    std::shared_ptr<Vector<RealT>> Get_Empty_Data_Vector() override
    {
      int data_dim = this->Error_Model()->Data_Dimension();
      std::shared_ptr<Std_Vector<RealT>> d = std::make_shared<Std_Vector<RealT>>(data_dim);
      return d;
    };

    void Compute_MAP_Point(std::shared_ptr<Vector<double>> &m_out)
    {
      // TODO: need to build the linear equation I want to solve...
    }
  };
}

#endif // OEDLIB_TEST_LINEAR_BAYESIAN_INVERSION_HPP
