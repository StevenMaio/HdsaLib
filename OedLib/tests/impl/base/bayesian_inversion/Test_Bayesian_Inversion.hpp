//
// Created by Steven Maio on 5/27/26.
//

#ifndef OEDLIB_TEST_BAYESIAN_INVERSION_HPP
#define OEDLIB_TEST_BAYESIAN_INVERSION_HPP

#include "../vectors/OED_Test_Vector.hpp"
#include "../../../../src/core/bayesian_inversion/OED_Bayesian_Inversion_Interface.hpp"

namespace OED_TEST
{
  class Test_Bayesian_Inversion : public OED::Bayesian_Inversion_Interface<double>
  {
  public:
    Test_Bayesian_Inversion(
        std::shared_ptr<OED::Likelihood_Interface<double>> likelihood,
        std::shared_ptr<OED::Prior_Interface<double>> prior,
        std::shared_ptr<OED::Constraint_Interface<double>> constraint
    )
      : Bayesian_Inversion_Interface(likelihood, prior, constraint) {}

    std::shared_ptr<Vector<double>> Get_Empty_Parameter_Vector() override
    {
      int param_dim = this->Prior()->Param_Dimension();
      std::shared_ptr<Test_Vector<double>> m = std::make_shared<Test_Vector<double>>(param_dim);
      return m;
    };

    std::shared_ptr<Vector<double>> Get_Empty_State_Vector() override
    {
      int state_dim = this->Constraint()->State_Dimension();
      std::shared_ptr<Test_Vector<double>> u = std::make_shared<Test_Vector<double>>(state_dim);
      return u;
    };

    std::shared_ptr<Vector<double>> Get_Empty_Data_Vector() override
    {
      int data_dim = this->Likelihood()->Data_Dimension();
      std::shared_ptr<Test_Vector<double>> d = std::make_shared<Test_Vector<double>>(data_dim);
      return d;
    };
  };
}

#endif //OEDLIB_TEST_BAYESIAN_INVERSION_HPP
