//
// Created by Steven Maio on 5/25/26.
//
#include <iostream>
#include <vector>
#include <Eigen/Dense>

#include "Test_Bayesian_Inversion.hpp"
#include "OED_Lazy_Greedy.hpp"
#include "OED_Linear_OED_D_Opt.hpp"
#include "OED_Test_Vector.hpp"
#include "Poisson_Constraint.hpp"
#include "Poisson_Likelihood.hpp"
#include "Poisson_Prior.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using OED_TEST::Test_Vector;

int main()
{
  int dim = 100;
  double norm_scale = 5.0 / 6;
  double grad_scale = 1.0 / 30;

  VectorXd m(dim);
  VectorXd u(dim);
  for (int i = 0; i < dim; i++)
  {
    m(i) = 2.0;
  }

  Test_Vector<double> param(dim);
  Test_Vector<double> state(dim);
  param.Vec() = m;

  auto constraint = std::make_shared<OED_TEST::Poisson_Constraint>(dim);
  constraint->State_Solve(state, param);

  u = state.Vec();
  std::cout << u << std::endl;
  std::vector<int> obs_vec;
  double noise_std = 1e-2;
  for (int i = 0; i < 12; i++)
  {
    obs_vec.push_back(i * 9);
  }
  int data_dim = obs_vec.size();

  auto likelihood = std::make_shared<OED_TEST::Poisson_Likelihood>(dim, noise_std, obs_vec);
  Test_Vector<double> data(data_dim);
  likelihood->Observation_Operator_Apply(data, state);
  std::cout << data.Vec() << std::endl;

  // test the noise precision apply -- looks good to me
  Test_Vector<double> test(data_dim);
  likelihood->Noise_Precision_Apply(test, data);
  std::cout << test.Vec() << std::endl;

  // Construct the rows of F
  auto prior = std::make_shared<OED_TEST::Poisson_Prior>(constraint, norm_scale, grad_scale);
  std::shared_ptr<OED::Bayesian_Inversion_Interface<double>> inversion_problem
      = std::make_shared<OED_TEST::Test_Bayesian_Inversion>(likelihood, prior, constraint);

  auto oed_problem = std::make_shared<OED::Linear_OED_D_Opt<double>>(inversion_problem);
  int budget = 5;
  OED::Active_Sensors design = OED::Lazy_Greedy_Solve(*oed_problem, data_dim, budget);
  std::cout << oed_problem->Evaluate(design) << std::endl;
  design.Print_Sensors();
}
