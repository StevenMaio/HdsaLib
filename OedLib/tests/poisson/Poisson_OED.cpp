//
// Created by Steven Maio on 5/25/26.
//
#include <iostream>
#include <vector>
#include "Eigen/Dense"

#include "Test_Linear_Bayesian_Inversion.hpp"
#include "OED_Lazy_Greedy.hpp"
#include "OED_Linear_OED_D_Opt.hpp"
#include "OED_Std_Vector.hpp"
#include "Poisson_Model.hpp"
#include "OED_Gaussian_Error.hpp"
#include "Poisson_Prior.hpp"
#include "Poisson_Obs.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using OED::Std_Vector;

int main()
{
  int dim = 100;
  double norm_scale = 5.0 / 6;
  double grad_scale = 1.0 / 30;

  VectorXd m(dim);
  VectorXd u(dim);
  m.setConstant(2.0);

  auto param = std::make_shared<Std_Vector<double>>(dim);
  auto state = std::make_shared<Std_Vector<double>>(dim);
  param->Vec() = m;

  std::cout << "Param:" << std::endl << m << std::endl << std::endl;

  auto model = std::make_shared<OED_TEST::Poisson_Model<double>>(dim);
  model->State_Solve(*state, *param);

  u = state->Vec();
  std::cout << "State:" << std::endl << u << std::endl;
  std::vector<int> obs_vec;
  double noise_std = 1e-2;
  for (int i = 0; i < 12; i++)
  {
    obs_vec.push_back(i * 9);
  }
  int data_dim = obs_vec.size();

  auto obs_op = std::make_shared<OED_TEST::Poisson_Observation_Operator<double>>(dim, obs_vec);
  auto error_model = std::make_shared<OED::Gaussian_Error<double>>(data_dim, noise_std);
  auto data = std::make_shared<Std_Vector<double>>(data_dim);
  obs_op->Observation_Operator_Apply(*data, *state);
  std::cout << data->Vec() << std::endl;

  // test the noise precision apply -- looks good to me
  Std_Vector<double> test(data_dim);
  error_model->Noise_Precision_Apply(test, *data);
  std::cout << test.Vec() << std::endl;

  // Construct the rows of F
  auto prior = std::make_shared<OED_TEST::Poisson_Prior<double>>(model, norm_scale, grad_scale);
  auto inversion_problem = std::make_shared<OED_TEST::Test_Linear_Bayesian_Inversion<double>>(model, obs_op, prior, error_model);

  auto oed_problem = std::make_shared<OED::Linear_OED_D_Opt<double>>(inversion_problem);
  int budget = 5;
  OED::Active_Sensors design = OED::Lazy_Greedy_Solve(*oed_problem, data_dim, budget);
  std::cout << oed_problem->Evaluate(design) << std::endl;
  design.Print_Sensors();
}
