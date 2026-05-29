//
// Created by Steven Maio on 5/25/26.
//

#ifndef OEDLIB_TEST_LINEAR_OPERATOR_HPP
#define OEDLIB_TEST_LINEAR_OPERATOR_HPP
#include "../../../../src/core/base/linear_algebra/OED_Linear_Operator.hpp"

#include <Eigen/src/SparseLU/SparseLU.h>

namespace OED_TEST
{
  template <class RealT>
  class Test_Linear_Operator : public OED::Linear_Operator<RealT>
  {
  };
}

#endif //OEDLIB_TEST_LINEAR_OPERATOR_HPP
