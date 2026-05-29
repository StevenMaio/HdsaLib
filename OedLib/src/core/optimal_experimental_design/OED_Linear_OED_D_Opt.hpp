//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LINEAR_OED_D_OPT_HPP
#define OEDLIB_LINEAR_OED_D_OPT_HPP

#include "../bayesian_inversion/OED_Bayesian_Inversion_Interface.hpp"
#include "../base/oed/OED_Discrete_Design_Criterion.hpp"
#include "../base/vectors/OED_Vector.hpp"
#include "../base/linear_algebra/OED_Dense_Matrix.hpp"

namespace OED
{
  template <class RealT>
  class Linear_OED_D_Opt : public Discrete_Design_Criterion
  {
  private:
    Bayesian_Inversion_Interface<RealT> &inversion_problem_;

    Dense_Matrix<RealT> Fm_cov_;
    Dense_Matrix<RealT> noise_cov_;

    public:
    Linear_OED_D_Opt(Bayesian_Inversion_Interface<RealT> &inversion_problem)
      : inversion_problem_(inversion_problem),
        Fm_cov_(inversion_problem.Likelihood().Data_Dimension(), inversion_problem.Likelihood().Data_Dimension()),
        noise_cov_(inversion_problem.Likelihood().Data_Dimension(), inversion_problem.Likelihood().Data_Dimension())
    {
      Constraint<RealT> &constraint = inversion_problem.Constraint();
      Likelihood_Model<RealT> &likelihood = inversion_problem.Likelihood();
      Prior_Model<RealT> &prior = inversion_problem.Prior();

      int data_dim = inversion_problem.Likelihood().Data_Dimension();

      // Construct Fm covariance matrix
      Vector<RealT> *temp_ptr = inversion_problem.Get_Empty_Data_Vector(); // This is only usable because we have a linear inverse problem
      Vector<RealT> &temp = *temp_ptr;
      for (int i = 0; i < data_dim; i++)
      {
        // TODO: fix this choice -- use smart pointers eventually
        Vector<RealT> *col_ptr = inversion_problem.Get_Empty_Data_Vector();
        Vector<RealT> *row_ptr = inversion_problem.Get_Empty_State_Vector();
        Vector<RealT> &col = *col_ptr;
        Vector<RealT> &row = *row_ptr;
        col.Set_Entry(i, 1.0);
        likelihood.Observation_Operator_Transpose_Apply(row, col);
        // TODO: maybe create some kind of Adjoint solve function
        constraint.c_u_Transpose_Inverse_Apply(row, row, temp, temp);
        constraint.c_z_Transpose_Apply(row, row, temp, temp);
        prior.Mass_Matrix_Inverse_Apply(row, row);
        prior.Prior_Covariance_Apply(row, row);
        constraint.State_Solve(row, row);
        likelihood.Observation_Operator_Apply(col, row);
        for (int j = 0; j < data_dim; j++)
        {
          this->Fm_cov_.Set_Entry(i, j, -col.Get_Entry(j));
        }
        delete col_ptr;
        delete row_ptr;
      }
      delete temp_ptr;
      // TODO: delete this later
      std::cout << this->Fm_cov_.Data() << std::endl;

      // TODO: construct noise covariance matrix
      for (int i = 0; i < data_dim; i++)
      {
        Vector<RealT> *row_ptr = inversion_problem.Get_Empty_Data_Vector();
        auto &row = *row_ptr;
        row.Set_Entry(i, 1.0);
        likelihood.Noise_Covariance_Apply(row, row);
        for (int j = 0; j < data_dim; j++)
        {
          this->noise_cov_.Set_Entry(j, i, row.Get_Entry(j));
        }
        delete row_ptr;
      }
      // TODO: delete this later
      std::cout << this->noise_cov_.Data() << std::endl;
    }

    double Evaluate(Active_Sensors &sensors) override
    {
      // TODO: fix this as well at some point
      Dense_Matrix<RealT> A = this->Fm_cov_.Select_Subsquare_Matrix(sensors.Selection());
      Dense_Matrix<RealT> sub_noise_cov = this->noise_cov_.Select_Subsquare_Matrix(sensors.Selection());
      sub_noise_cov.Right_Inverse_Multiply(A, A);
      for (int i = 0; i < sensors.Selection().size(); i++)
      {
        // A(i, i) += 1;
        A.Set_Entry(i, i, A(i, i) + 1);
      }
      return 0.5 * std::log(A.Compute_Determinant());
    };

    double Compute_Marginal_Gain(Active_Sensors &sensors, int v) override
    {
      // TODO: implement this
      return 0;
    };
  };
}

#endif //OEDLIB_LINEAR_OED_D_OPT_H
