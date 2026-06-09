#include <iostream>
#include <vector>
#include "Eigen/Dense"

#include "Test_Linear_Bayesian_Inversion.hpp"
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
  // Set up the problem
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
  auto prior = std::make_shared<OED_TEST::Poisson_Prior>(constraint, norm_scale, grad_scale);
  std::shared_ptr<OED::Bayesian_Inversion_Interface<double>> inversion_problem
      = std::make_shared<OED_TEST::Test_Linear_Bayesian_Inversion>(likelihood, prior, constraint);

  // TODO: need to do actual inversion
}
