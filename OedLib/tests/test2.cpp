//
// Created by Steven Maio on 5/25/26.
//
#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>

#include "impl/base/vectors/OED_Test_Vector.hpp"
#include "poisson/Poisson_Constraint.hpp"
#include "poisson/Poisson_Likelihood.hpp"
#include "poisson/Poisson_Prior.hpp"

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

  OED_TEST::Poisson_Constraint constraint(dim);
  constraint.State_Solve(state, param);

  u = state.Vec();
  std::cout << u << std::endl;
  std::vector<int> obs_vec;
  double noise_std = 1e-2;
  for (int i = 0; i < 12; i++)
  {
    obs_vec.push_back(i * 9);
  }
  int data_dim = obs_vec.size();

  OED_TEST::Poisson_Likelihood likelihood(dim, noise_std, obs_vec);
  Test_Vector<double> data(data_dim);
  likelihood.Observation_Operator_Apply(data, state);
  std::cout << data.Vec() << std::endl;

  // test the noise precision apply -- looks good to me
  Test_Vector<double> test(data_dim);
  likelihood.Noise_Precision_Apply(test, data);
  std::cout << test.Vec() << std::endl;

  // Construct the rows of F
  MatrixXd F(data_dim, dim);
  OED_TEST::Poisson_Prior prior(constraint, norm_scale, grad_scale);

  for (int i = 0; i < data_dim; i++)
  {
    Test_Vector<double> x(data_dim);
    Test_Vector<double> row(dim);
    x.Vec()(i) = 1.0;
    likelihood.Observation_Operator_Transpose_Apply(row, x);
    constraint.c_u_Transpose_Inverse_Apply(row, row, state, param);
    constraint.c_z_Transpose_Apply(row, row, state, param);
    prior.Mass_Matrix_Inverse_Apply(row, row);
    prior.Prior_Factor_Apply(row, row);
    row.Vec() *= - 1 / noise_std;
    for (int j = 0; j < dim; j++)
    {
      F(i, j) = row.Vec()(j);
    }
  }

  // Build matrix to move log det into data space
  MatrixXd Fm_cov(data_dim, data_dim);
  for (int i = 0; i < data_dim; i++)
  {
    // TODO: this suggests we need to be able to create zero vectors of some arbitrary size...
    Test_Vector<double> x(data_dim);
    Test_Vector<double> row(dim);
    x.Vec()(i) = 1.0;
    likelihood.Observation_Operator_Transpose_Apply(row, x);
    constraint.c_u_Transpose_Inverse_Apply(row, row, state, param);
    constraint.c_z_Transpose_Apply(row, row, state, param);
    row.Vec() *= - 1;
    prior.Mass_Matrix_Inverse_Apply(row, row);
    prior.Prior_Covariance_Apply(row, row);
    constraint.State_Solve(row, row);
    likelihood.Observation_Operator_Apply(x, row);
    for (int j = 0; j < data_dim; j++)
    {
      Fm_cov(i, j) = x.Vec()(j);
    }
  }

  std::cout << Fm_cov << std::endl;
}
