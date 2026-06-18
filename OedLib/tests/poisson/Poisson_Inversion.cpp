#include <iostream>
#include <vector>
#include "Eigen/Dense"

#include "Test_Linear_Bayesian_Inversion.hpp"
#include "OED_Lazy_Greedy.hpp"
#include "OED_Linear_OED_D_Opt.hpp"
#include "OED_Std_Vector.hpp"
#include "Poisson_Model.hpp"
#include "Poisson_Likelihood.hpp"
#include "Poisson_Prior.hpp"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using OED::Std_Vector;

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

  Std_Vector<double> param(dim);
  Std_Vector<double> state(dim);
  param.Vec() = m;

  auto model = std::make_shared<OED_TEST::Poisson_Model<double>>(dim);
  model->State_Solve(state, param);

  u = state.Vec();
  std::cout << u << std::endl;
  std::vector<int> obs_vec;
  double noise_std = 1e-2;
  for (int i = 0; i < 12; i++)
  {
    obs_vec.push_back(i * 9);
  }
  int data_dim = obs_vec.size();

  auto likelihood = std::make_shared<OED_TEST::Poisson_Likelihood<double>>(dim, noise_std, obs_vec);
  auto prior = std::make_shared<OED_TEST::Poisson_Prior<double>>(model, norm_scale, grad_scale);
  auto inversion_problem = std::make_shared<OED_TEST::Test_Linear_Bayesian_Inversion<double>>(likelihood, prior, model);

  // create data and intialize map point
  auto data = inversion_problem->Get_Empty_Data_Vector();
  likelihood->Observation_Operator_Apply(*data, (OED::Vector<double> &) state);
  // inversion_problem->Set_Data(data);

  auto map_estimate = inversion_problem->Get_Empty_Parameter_Vector();
  inversion_problem->Compute_MAP_Point(map_estimate);
  // TODO: need to do actual inversion

  // TODO: create the forward map

  // TODO: maybe do a weighted inner product decomposition?
}
