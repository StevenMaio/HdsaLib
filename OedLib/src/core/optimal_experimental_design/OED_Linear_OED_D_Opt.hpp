//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LINEAR_OED_D_OPT_HPP
#define OEDLIB_LINEAR_OED_D_OPT_HPP

#include <memory>

#include "OED_Bayesian_Inversion_Interface.hpp"
#include "OED_Discrete_Design_Criterion.hpp"
#include "OED_Vector.hpp"
#include "OED_Dense_Matrix.hpp"

namespace OED
{
  template <class RealT>
  class Linear_OED_D_Opt : public Discrete_Design_Criterion
  {
    private:
    std::shared_ptr<Bayesian_Inversion_Interface<RealT>> inversion_problem_;

    std::shared_ptr<Dense_Matrix<RealT>> forward_cov_;
    std::shared_ptr<Dense_Matrix<RealT>> noise_cov_;

    public:
    Linear_OED_D_Opt(std::shared_ptr<Bayesian_Inversion_Interface<RealT>> &inversion_problem)
      : inversion_problem_(inversion_problem)
    {
      int data_dim = inversion_problem->Likelihood()->Data_Dimension();
      this->forward_cov_ = std::make_shared<Dense_Matrix<RealT>>(data_dim, data_dim);
      this->noise_cov_ = std::make_shared<Dense_Matrix<RealT>>(data_dim, data_dim);

      this->Construct_Forward_Covariance();
      this->Construct_Noise_Covariance();

      // TODO: delete these later
      std::cout << this->forward_cov_->Data() << std::endl;
      std::cout << this->noise_cov_->Data() << std::endl;
    }

    double Evaluate(Active_Sensors &sensors) override
    {
      std::shared_ptr<Dense_Matrix<RealT>> A = this->forward_cov_->Select_Subsquare_Matrix(sensors.Selection());
      std::shared_ptr<Dense_Matrix<RealT>> sub_noise_cov = this->noise_cov_->Select_Subsquare_Matrix(sensors.Selection());
      sub_noise_cov->Right_Inverse_Multiply(*A, *A);
      for (int i = 0; i < sensors.Selection().size(); i++)
      {
        A->Set_Entry(i, i, (*A)(i, i) + 1);
      }
      return 0.5 * std::log(A->Compute_Determinant());
    };

    double Compute_Marginal_Gain(Active_Sensors &sensors, int v) override
    {
      // TODO: implement this
      return 0;
    };
    private:

    void Construct_Forward_Covariance()
    {
      auto &inversion_problem = this->inversion_problem_;
      std::shared_ptr<Likelihood_Interface<RealT>> &likelihood = inversion_problem->Likelihood();
      std::shared_ptr<Constraint_Interface<RealT>> &constraint = inversion_problem->Constraint();
      std::shared_ptr<Prior_Interface<RealT>> &prior = inversion_problem->Prior();
      int data_dim = likelihood->Data_Dimension();

      std::shared_ptr<Vector<RealT>> temp = inversion_problem->Get_Empty_Data_Vector(); // This is only usable because we have a linear inverse problem
      std::shared_ptr<Vector<RealT>> col = inversion_problem->Get_Empty_Data_Vector();
      std::shared_ptr<Vector<RealT>> row = inversion_problem->Get_Empty_State_Vector();
      for (int i = 0; i < data_dim; i++)
      {
        col->Set_Entry(i, 1.0);
        likelihood->Observation_Operator_Transpose_Apply(*row, *col);
        // TODO: maybe create some kind of Adjoint solve function
        constraint->c_u_Transpose_Inverse_Apply(*row, *row, *temp, *temp);
        constraint->c_z_Transpose_Apply(*row, *row, *temp, *temp);
        prior->Mass_Matrix_Inverse_Apply(*row, *row);
        prior->Prior_Covariance_Apply(*row, *row);
        constraint->State_Solve(*row, *row);
        likelihood->Observation_Operator_Apply(*col, *row);
        for (int j = 0; j < data_dim; j++)
        {
          this->forward_cov_->Set_Entry(j, i, -col->Get_Entry(j));
        }
        row->Zeros();
        col->Zeros();
      }

    }

    inline void Construct_Noise_Covariance()
    {
      auto &inversion_problem = this->inversion_problem_;
      std::shared_ptr<Likelihood_Interface<RealT>> &likelihood = inversion_problem->Likelihood();
      std::shared_ptr<Vector<RealT>> row = inversion_problem->Get_Empty_State_Vector();
      int data_dim = likelihood->Data_Dimension();

      for (int i = 0; i < data_dim; i++)
      {
        row->Set_Entry(i, 1.0);
        likelihood->Noise_Covariance_Apply(*row, *row);
        for (int j = 0; j < data_dim; j++)
        {
          this->noise_cov_->Set_Entry(j, i, row->Get_Entry(j));
        }
        row->Zeros();
      }
    }
  };
}

#endif //OEDLIB_LINEAR_OED_D_OPT_H
