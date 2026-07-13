//
// Created by Steven Maio on 5/22/26.
//

#ifndef OEDLIB_LINEAR_OED_D_OPT_HPP
#define OEDLIB_LINEAR_OED_D_OPT_HPP

#include <memory>

#include "OED_Bayesian_Inversion_Interface.hpp"
#include "OED_Transient_Observation_Operator.hpp"
#include "OED_Discrete_Design_Criterion.hpp"
#include "OED_Vector.hpp"
#include "OED_Dense_Matrix.hpp"

#include "OED_Ptr.hpp"

namespace OED
{
  template <class RealT>
  class Linear_OED_D_Opt : public Discrete_Design_Criterion
  {
    private:
    Ptr<Bayesian_Inversion_Interface<RealT>> inversion_problem_;

    Ptr<Dense_Matrix<RealT>> forward_cov_;
    Ptr<Dense_Matrix<RealT>> noise_cov_;

    public:
    Linear_OED_D_Opt(Ptr<Bayesian_Inversion_Interface<RealT>> inversion_problem)
      : inversion_problem_(inversion_problem)
    {
      int data_dim = inversion_problem->Error_Model()->Data_Dimension();
      this->forward_cov_ = OED::makePtr<Dense_Matrix<RealT>>(data_dim, data_dim);
      this->noise_cov_ = OED::makePtr<Dense_Matrix<RealT>>(data_dim, data_dim);

      this->Construct_Forward_Covariance();
      this->Construct_Noise_Covariance();

      // TODO: delete these later
      std::cout << this->forward_cov_->Data() << std::endl;
      std::cout << this->noise_cov_->Data() << std::endl;
    }

    double Evaluate(Active_Sensors &sensors) override
    {
      std::vector<int> selection = this->Compute_Selection_Indices(sensors);
      Ptr<Dense_Matrix<RealT>> A = this->forward_cov_->Select_Subsquare_Matrix(selection);
      Ptr<Dense_Matrix<RealT>> sub_noise_cov = this->noise_cov_->Select_Subsquare_Matrix(selection);
      sub_noise_cov->Right_Inverse_Multiply(*A, *A);
      for (int i = 0; i < selection.size(); i++)
      {
        A->Set_Entry(i, i, (*A)(i, i) + 1);
      }
      return 0.5 * std::log(A->Compute_Determinant());
    }

    // double Evaluate(Active_Sensors &sensors) override
    // {
    //   std::vector<int> selection = this->Compute_Selection_Indices(sensors);
    //   // TODO: determine indices based on sensors and problem type
    //   Ptr<Dense_Matrix<RealT>> A = this->forward_cov_->Select_Subsquare_Matrix(selection);
    //   Ptr<Dense_Matrix<RealT>> sub_noise_cov = this->noise_cov_->Select_Subsquare_Matrix(selection);
    //   A->Scaled_Plus(sub_noise_cov);
    // }

    double Compute_Marginal_Gain(Active_Sensors &sensors, int v) override
    {
      // TODO: implement this
      return 0;
    }

    private:

    inline void Construct_Forward_Covariance()
    {
      Ptr<Bayesian_Inversion_Interface<RealT>> &inversion_problem = this->inversion_problem_;
      Ptr<Model_Interface<RealT>> &model = inversion_problem->Model();
      Ptr<Observation_Operator_Interface<RealT>> &obs_op = inversion_problem->Observation_Operator();
      Ptr<Prior_Interface<RealT>> &prior = inversion_problem->Prior();
      Ptr<Error_Model_Interface<RealT>> &error_model = inversion_problem->Error_Model();
      int data_dim = error_model->Data_Dimension();

      // Required to use methods, but not actually used here
      Ptr<Vector<RealT>> d = obs_op->Get_Empty_Data_Vector();
      Ptr<Vector<RealT>> u_temp = model->Get_Empty_State_Vector();
      Ptr<Vector<RealT>> u1 = model->Get_Empty_State_Vector();
      Ptr<Vector<RealT>> m_temp = model->Get_Empty_Parameter_Vector();
      Ptr<Vector<RealT>> m1 = model->Get_Empty_Parameter_Vector();
      Ptr<Vector<RealT>> m2 = model->Get_Empty_Parameter_Vector();

      u_temp->Zeros();
      m_temp->Zeros();
      for (int i = 0; i < data_dim; i++)
      {
        // clear vectors
        d->Zeros();
        u1->Zeros();
        m1->Zeros();
        m2->Zeros();

        d->Set_Entry(i, 1.0);
        obs_op->Observation_Operator_Transpose_Apply(*u1, *d);
        model->State_Transpose_Apply(*m1, *u1, *u_temp, *m_temp);
        prior->Mass_Matrix_Inverse_Apply(*m2, *m1);
        prior->Prior_Covariance_Apply(*m1, *m2);
        model->State_Solve(*u1, *m1);
        obs_op->Observation_Operator_Apply(*d, *u1);
        for (int j = 0; j < data_dim; j++)
        {
          this->forward_cov_->Set_Entry(j, i, d->Get_Entry(j));
        }
      }
      std::cout << "Linear_OED_D_Opt::Construct_Forward_Covariance finished constructing forward covariance matrix" << std::endl;;
    }

    inline void Construct_Noise_Covariance()
    {
      auto &inversion_problem = this->inversion_problem_;
      Ptr<Error_Model_Interface<RealT>> &error_model = inversion_problem->Error_Model();
      Ptr<Vector<RealT>> d_in = inversion_problem->Get_Empty_State_Vector();
      Ptr<Vector<RealT>> d_out = inversion_problem->Get_Empty_State_Vector();
      int data_dim = error_model->Data_Dimension();

      for (int i = 0; i < data_dim; i++)
      {
        d_in->Zeros();
        d_out->Zeros();
        d_in->Set_Entry(i, 1.0);
        error_model->Noise_Covariance_Apply(*d_out, *d_in);
        for (int j = 0; j < data_dim; j++)
        {
          this->noise_cov_->Set_Entry(j, i, d_out->Get_Entry(j));
        }
      }
      std::cout << "Linear_OED_D_Opt::Construct_Forward_Covariance finished constructing noise covariance matrix" << std::endl;;
    }

    private:

    inline std::vector<int> Compute_Selection_Indices(Active_Sensors &sensors)
    {
      std::vector<int> selection;
      if (!this->inversion_problem_->Is_Transient())
      {
        selection = sensors.Selection();
      }
      else
      {
        auto &trans_obs = dynamic_cast<Transient_Observation_Operator<RealT> &>(*this->inversion_problem_->Observation_Operator());
        int num_t = trans_obs.Get_Num_T();
        int num_sensors = trans_obs.Num_Sensors();
        for (int v: sensors.Selection())
        {
          for (int i = 0; i < num_t; i++)
          {
            selection.push_back(v + i * num_sensors);
          }
        }
      }
      return selection;
    }
  };
}

#endif //OEDLIB_LINEAR_OED_D_OPT_H
