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
    Test_Bayesian_Inversion(OED::Likelihood_Model<double> &likelihood, OED::Prior_Model<double> &prior, OED::Constraint<double> &constraint)
      : Bayesian_Inversion_Interface(likelihood, prior, constraint) {}

    // TODO: will deal with this typing issue later
    Vector<double> *Get_Empty_Parameter_Vector() override
    {
      int param_dim = this->Prior().Param_Dimension();
      // I KNOW I'M COMMITTING A CRIME HERE
      auto *m = new Test_Vector<double>(param_dim);
      return m;
    };

    Vector<double> *Get_Empty_State_Vector() override
    {
      int state_dim = this->Constraint().State_Dimension();
      // I KNOW I'M COMMITTING A CRIME HERE
      auto *u = new Test_Vector<double>(state_dim);
      return u;
    };

    Vector<double> *Get_Empty_Data_Vector() override
    {
      int data_dim = this->Likelihood().Data_Dimension();
      // I KNOW I'M COMMITTING A CRIME HERE
      auto *d = new Test_Vector<double>(data_dim);
      return d;
    };
  };
}

#endif //OEDLIB_TEST_BAYESIAN_INVERSION_HPP
