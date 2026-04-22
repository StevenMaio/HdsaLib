/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_PC_ADV_DIFF_CONSTRAINT_HPP
#define HDSA_PC_ADV_DIFF_CONSTRAINT_HPP

#include "HDSA_Std_Vector.hpp"
#include "HDSA_Dense_Matrix.hpp"
#include "HDSA_Linear_Algebra.hpp"

template <class RealT>
class Adv_Diff_Constraint
{

public:
  int m_;                                         // Mesh resolution
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x_;        // Mesh nodes on [0,1]
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S_;        // Stiffness matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> M_;        // Mass matrix
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> forcing_;  // forcing_vector
  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> robin_bc_; // Robin boundary condition matrix
  RealT alpha_;                                   // Robin coefficient

  void Dense_Linear_Solve(const HDSA::Dense_Matrix<RealT> &A, HDSA::Dense_Matrix<RealT> &x, const HDSA::Dense_Matrix<RealT> &b) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> U = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> S = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Linear_Algebra::SVD<RealT>(A, *U, *VT, *S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    U->Multiply(*tmp, b, true, false);
    for (int k = 0; k < m_; k++)
    {
      RealT val = (*tmp)(k, 0) / (*S)(k, 0);
      tmp->Set_Entry(k, 0, val);
    }
    VT->Multiply(x, *tmp, true, false);
  }

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u);
    const HDSA::Std_Vector<RealT> &z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly(*z_std.get_std_vec());
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = Velocity_Assembly(*theta_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> A = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = (*D)(i, j) + (*V)(i, j) + alpha_ * (*robin_bc_)(i, j);
        A->Set_Entry(i, j, val);
      }
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    Dense_Linear_Solve(*A, *u_tmp, *forcing_);

    for (int k = 0; k < m_; k++)
    {
      u_std.Set_Entry(k, (*u_tmp)(k, 0));
    }
  }

  void c_u_Transpose_Inverse_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &u_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u_in);
    // const HDSA::Std_Vector<RealT>& u_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(u);
    const HDSA::Std_Vector<RealT> &z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly(*z_std.get_std_vec());
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = Velocity_Assembly(*theta_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> c_u_trans = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = (*D)(i, j) + (*V)(i, j) + alpha_ * (*robin_bc_)(i, j);
        c_u_trans->Set_Entry(j, i, val);
      }
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      u_tmp1->Set_Entry(k, 0, u_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    Dense_Linear_Solve(*c_u_trans, *u_tmp2, *u_tmp1);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp2)(k, 0));
    }
  }

  void c_u_Inverse_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &u_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u_in);
    // const HDSA::Std_Vector<RealT>& u_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(u);
    const HDSA::Std_Vector<RealT> &z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly(*z_std.get_std_vec());
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = Velocity_Assembly(*theta_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> c_u = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    for (int i = 0; i < m_; i++)
    {
      for (int j = 0; j < m_; j++)
      {
        RealT val = (*D)(i, j) + (*V)(i, j) + alpha_ * (*robin_bc_)(i, j);
        c_u->Set_Entry(i, j, val);
      }
    }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      u_tmp1->Set_Entry(k, 0, u_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    Dense_Linear_Solve(*c_u, *u_tmp2, *u_tmp1);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp2)(k, 0));
    }
  }

  void c_z_Transpose_Apply(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &z_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(z_out);
    const HDSA::Std_Vector<RealT> &u_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u_in);
    const HDSA::Std_Vector<RealT> &u_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly_z_Jacobian(*u_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      u_tmp->Set_Entry(k, 0, u_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> z_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    D->Multiply(*z_tmp, *u_tmp, true);

    for (int k = 0; k < m_; k++)
    {
      z_out_std.Set_Entry(k, (*z_tmp)(k, 0));
    }
  }

  void c_z_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &z_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z_in);
    const HDSA::Std_Vector<RealT> &u_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly_z_Jacobian(*u_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> z_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      z_tmp->Set_Entry(k, 0, z_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    D->Multiply(*u_tmp, *z_tmp);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp)(k, 0));
    }
  }

  void c_uu_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    u_out.Zeros();
  }

  void c_uz_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &z_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z_in);
    // const HDSA::Std_Vector<RealT>& u_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    const HDSA::Std_Vector<RealT> &lambda_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(lambda);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly_z_Jacobian(*lambda_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> z_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      z_tmp->Set_Entry(k, 0, z_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    D->Multiply(*u_tmp, *z_tmp);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp)(k, 0));
    }
  }

  void c_zu_Apply(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &z_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(z_out);
    const HDSA::Std_Vector<RealT> &u_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u_in);
    // const HDSA::Std_Vector<RealT>& u_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    const HDSA::Std_Vector<RealT> &lambda_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(lambda);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = Diff_Assembly_z_Jacobian(*lambda_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      u_tmp->Set_Entry(k, 0, u_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> z_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    D->Multiply(*z_tmp, *u_tmp, true);

    for (int k = 0; k < m_; k++)
    {
      z_out_std.Set_Entry(k, (*z_tmp)(k, 0));
    }
  }

  void c_zz_Apply(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    z_out.Zeros();
  }

  void c_theta_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &theta_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta_in);
    const HDSA::Std_Vector<RealT> &u_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = Velocity_Assembly_theta_Jacobian(*u_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> theta_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      theta_tmp->Set_Entry(k, 0, theta_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    V->Multiply(*u_tmp, *theta_tmp);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp)(k, 0));
    }
  }

  void c_ztheta_Apply(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    z_out.Zeros();
  }

  void c_utheta_Apply(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &lambda, const HDSA::Vector<RealT> &theta) const
  {
    HDSA::Std_Vector<RealT> &u_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(u_out);
    const HDSA::Std_Vector<RealT> &theta_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta_in);
    // const HDSA::Std_Vector<RealT>& u_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(u);
    // const HDSA::Std_Vector<RealT>& z_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(z);
    const HDSA::Std_Vector<RealT> &lambda_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(lambda);
    // const HDSA::Std_Vector<RealT>& theta_std = dynamic_cast<const HDSA::Std_Vector<RealT>&>(theta);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = Velocity_Assembly_utheta_Hessian(*lambda_std.get_std_vec());

    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> theta_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      theta_tmp->Set_Entry(k, 0, theta_in_std(k));
    }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> u_tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    V->Multiply(*u_tmp, *theta_tmp);

    for (int k = 0; k < m_; k++)
    {
      u_out_std.Set_Entry(k, (*u_tmp)(k, 0));
    }
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Diff_Assembly(const std::vector<RealT> &z) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    RealT h = (*x_)(1, 0) - (*x_)(0, 0);
    std::vector<RealT> phi_down_Dot = {-1.0 / h, -1.0 / h};
    std::vector<RealT> phi_up_Dot = {1.0 / h, 1.0 / h};
    RealT val = 0.0;

    for (int i = 1; i <= m_; ++i)
    {
      if (i > 1)
      {
        RealT x1 = (i - 2) * h + (h / 2) * (-1.0 / std::sqrt(3.0) + 1.0);
        RealT x2 = (i - 2) * h + (h / 2) * (1.0 / std::sqrt(3.0) + 1.0);
        std::vector<RealT> perm = {Diffusion_Coeff(x1, z), Diffusion_Coeff(x2, z)};
        val = (h / 2) * (phi_up_Dot[0] * phi_down_Dot[0] * perm[0] + phi_up_Dot[1] * phi_down_Dot[1] * perm[1]);
        D->Set_Entry(i - 2, i - 1, val);
        val = (h / 2) * (phi_up_Dot[0] * phi_up_Dot[0] * perm[0] + phi_up_Dot[1] * phi_up_Dot[1] * perm[1]);
        D->Set_Entry(i - 1, i - 1, val);
      }
      if (i < m_)
      {
        RealT x1 = (i - 1) * h + (h / 2) * (-1.0 / std::sqrt(3.0) + 1.0);
        RealT x2 = (i - 1) * h + (h / 2) * (1.0 / std::sqrt(3.0) + 1.0);
        std::vector<RealT> perm = {Diffusion_Coeff(x1, z), Diffusion_Coeff(x2, z)};
        val = (*D)(i - 1, i - 1) + (h / 2) * (phi_down_Dot[0] * phi_down_Dot[0] * perm[0] + phi_down_Dot[1] * phi_down_Dot[1] * perm[1]);
        D->Set_Entry(i - 1, i - 1, val);
        val = (h / 2) * (phi_up_Dot[0] * phi_down_Dot[0] * perm[0] + phi_up_Dot[1] * phi_down_Dot[1] * perm[1]);
        D->Set_Entry(i, i - 1, val);
      }
    }
    return D;
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Diff_Assembly_z_Jacobian(const std::vector<RealT> &u) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> D_diff = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    RealT h = (*x_)(1, 0) - (*x_)(0, 0);
    RealT x1 = (h / 2) * (-1 / std::sqrt(3) + 1);
    RealT x2 = (h / 2) * (1 / std::sqrt(3) + 1);
    std::vector<RealT> phi_down = {x2 / h, x1 / h};
    std::vector<RealT> phi_up = {x1 / h, x2 / h};
    std::vector<RealT> phi_down_Dot = {-1.0 / h, -1.0 / h};
    std::vector<RealT> phi_up_Dot = {1.0 / h, 1.0 / h};
    RealT val = 0.0;

    for (int i = 0; i < m_; ++i)
    {
      if (i > 0)
      {
        std::vector<RealT> u_prime = {u[i - 1] * phi_down_Dot[0] + u[i] * phi_up_Dot[0],
                                      u[i - 1] * phi_down_Dot[1] + u[i] * phi_up_Dot[1]};
        val = (h / 2) * (phi_up[0] * phi_down_Dot[0] * u_prime[0] + phi_up[1] * phi_down_Dot[1] * u_prime[1]);
        D_diff->Set_Entry(i - 1, i, val);
        val = (h / 2) * (phi_up[0] * phi_up_Dot[0] * u_prime[0] + phi_up[1] * phi_up_Dot[1] * u_prime[1]);
        D_diff->Set_Entry(i, i, val);
      }
      if (i < m_ - 1)
      {
        std::vector<RealT> u_prime = {u[i] * phi_down_Dot[0] + u[i + 1] * phi_up_Dot[0],
                                      u[i] * phi_down_Dot[1] + u[i + 1] * phi_up_Dot[1]};
        val = (*D_diff)(i, i) + (h / 2) * (phi_down[0] * phi_down_Dot[0] * u_prime[0] + phi_down[1] * phi_down_Dot[1] * u_prime[1]);
        D_diff->Set_Entry(i, i, val);
        val = (h / 2) * (phi_down[0] * phi_up_Dot[0] * u_prime[0] + phi_down[1] * phi_up_Dot[1] * u_prime[1]);
        D_diff->Set_Entry(i + 1, i, val);
      }
    }

    return D_diff;
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Velocity_Assembly(const std::vector<RealT> &theta) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    RealT h = (*x_)(1, 0) - (*x_)(0, 0);
    std::vector<RealT> phi_down_Dot = {-1.0 / h, -1.0 / h};
    std::vector<RealT> phi_up_Dot = {1.0 / h, 1.0 / h};
    RealT x1 = (h / 2) * (-1.0 / std::sqrt(3) + 1);
    RealT x2 = (h / 2) * (1.0 / std::sqrt(3) + 1);
    std::vector<RealT> phi_down = {x2 / h, x1 / h};
    std::vector<RealT> phi_up = {x1 / h, x2 / h};
    RealT val = 0.0;

    for (int i = 1; i <= m_; ++i)
    {
      if (i > 1)
      {
        x1 = (i - 2) * h + (h / 2) * (-1.0 / std::sqrt(3) + 1);
        x2 = (i - 2) * h + (h / 2) * (1.0 / std::sqrt(3) + 1);
        std::vector<RealT> vel = {Velocity_Coeff(x1, theta), Velocity_Coeff(x2, theta)};
        val = (*V)(i - 1, i - 1) + (h / 2.0) * (phi_up_Dot[0] * phi_down[0] * vel[0] + phi_up_Dot[1] * phi_down[1] * vel[1]);
        V->Set_Entry(i - 2, i - 1, val);
        val = (*V)(i - 1, i - 1) + (h / 2.0) * (phi_up_Dot[0] * phi_up[0] * vel[0] + phi_up_Dot[1] * phi_up[1] * vel[1]);
        V->Set_Entry(i - 1, i - 1, val);
      }
      if (i < m_)
      {
        x1 = (i - 1) * h + (h / 2) * (-1.0 / std::sqrt(3) + 1);
        x2 = (i - 1) * h + (h / 2) * (1.0 / std::sqrt(3) + 1);
        std::vector<RealT> vel = {Velocity_Coeff(x1, theta), Velocity_Coeff(x2, theta)};
        val = (*V)(i - 1, i - 1) + (h / 2.0) * (phi_down_Dot[0] * phi_down[0] * vel[0] + phi_down_Dot[1] * phi_down[1] * vel[1]);
        V->Set_Entry(i - 1, i - 1, val);
        val = (h / 2.0) * (phi_down_Dot[0] * phi_up[0] * vel[0] + phi_down_Dot[1] * phi_up[1] * vel[1]);
        V->Set_Entry(i, i - 1, val);
      }
    }

    return V;
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Velocity_Assembly_theta_Jacobian(const std::vector<double> &u) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V_diff = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    RealT h = (*x_)(1, 0) - (*x_)(0, 0);
    RealT x1 = (h / 2) * (-1 / std::sqrt(3) + 1);
    RealT x2 = (h / 2) * (1 / std::sqrt(3) + 1);
    std::vector<RealT> phi_down = {x2 / h, x1 / h};
    std::vector<RealT> phi_up = {x1 / h, x2 / h};
    std::vector<RealT> phi_down_Dot = {-1.0 / h, -1.0 / h};
    std::vector<RealT> phi_up_Dot = {1.0 / h, 1.0 / h};
    RealT val = 0.0;

    for (int i = 0; i < m_; ++i)
    {
      if (i > 0)
      {
        std::vector<RealT> u_prime = {u[i - 1] * phi_down_Dot[0] + u[i] * phi_up_Dot[0],
                                      u[i - 1] * phi_down_Dot[1] + u[i] * phi_up_Dot[1]};
        val = (h / 2) * (phi_up[0] * phi_down[0] * u_prime[0] + phi_up[1] * phi_down[1] * u_prime[1]);
        V_diff->Set_Entry(i - 1, i, val);
        val = (h / 2) * (phi_up[0] * phi_up[0] * u_prime[0] + phi_up[1] * phi_up[1] * u_prime[1]);
        V_diff->Set_Entry(i, i, val);
      }
      if (i < m_ - 1)
      {
        std::vector<RealT> u_prime = {u[i] * phi_down_Dot[0] + u[i + 1] * phi_up_Dot[0],
                                      u[i] * phi_down_Dot[1] + u[i + 1] * phi_up_Dot[1]};
        val = (*V_diff)(i, i) + (h / 2) * (phi_down[0] * phi_down[0] * u_prime[0] + phi_down[1] * phi_down[1] * u_prime[1]);
        V_diff->Set_Entry(i, i, val);
        val = (h / 2) * (phi_down[0] * phi_up[0] * u_prime[0] + phi_down[1] * phi_up[1] * u_prime[1]);
        V_diff->Set_Entry(i + 1, i, val);
      }
    }

    return V_diff;
  }

  HDSA::Ptr<HDSA::Dense_Matrix<RealT>> Velocity_Assembly_utheta_Hessian(const std::vector<double> &lambda) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> V_hess = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    RealT h = (*x_)(1, 0) - (*x_)(0, 0);
    RealT x1 = (h / 2) * (-1 / std::sqrt(3) + 1);
    RealT x2 = (h / 2) * (1 / std::sqrt(3) + 1);
    std::vector<RealT> phi_down = {x2 / h, x1 / h};
    std::vector<RealT> phi_up = {x1 / h, x2 / h};
    std::vector<RealT> phi_down_Dot = {-1.0 / h, -1.0 / h};
    std::vector<RealT> phi_up_Dot = {1.0 / h, 1.0 / h};
    RealT val = 0.0;

    for (int i = 0; i < m_; ++i)
    {
      if (i > 0)
      {
        std::vector<RealT> lam = {lambda[i - 1] * phi_down[0] + lambda[i] * phi_up[0],
                                  lambda[i - 1] * phi_down[1] + lambda[i] * phi_up[1]};
        val = (h / 2) * (phi_up[0] * phi_down_Dot[0] * lam[0] + phi_up[1] * phi_down_Dot[1] * lam[1]);
        V_hess->Set_Entry(i - 1, i, val);
        val = (h / 2) * (phi_up[0] * phi_up_Dot[0] * lam[0] + phi_up[1] * phi_up_Dot[1] * lam[1]);
        V_hess->Set_Entry(i, i, val);
      }
      if (i < m_ - 1)
      {
        std::vector<RealT> lam = {lambda[i] * phi_down[0] + lambda[i + 1] * phi_up[0],
                                  lambda[i] * phi_down[1] + lambda[i + 1] * phi_up[1]};
        val = (*V_hess)(i, i) + (h / 2) * (phi_down[0] * phi_down_Dot[0] * lam[0] + phi_down[1] * phi_down_Dot[1] * lam[1]);
        V_hess->Set_Entry(i, i, val);
        val = (h / 2) * (phi_down[0] * phi_up_Dot[0] * lam[0] + phi_down[1] * phi_up_Dot[1] * lam[1]);
        V_hess->Set_Entry(i + 1, i, val);
      }
    }

    return V_hess;
  }

  RealT Diffusion_Coeff(const RealT &x_pt, const std::vector<RealT> &z) const
  {
    RealT val = Interpolation_1D(x_pt, z);
    return val;
  }

  RealT Velocity_Coeff(const RealT &x_pt, const std::vector<RealT> &theta) const
  {
    RealT val = Interpolation_1D(x_pt, theta);
    return val;
  }

  RealT Interpolation_1D(const RealT &x_pt, const std::vector<RealT> &y) const
  {
    int i = 0;
    if (x_pt >= (*x_)(m_ - 2, 0))
    {
      i = m_ - 2;
    }
    else
    {
      while (x_pt > (*x_)(i + 1, 0))
      {
        i++;
      }
    }
    RealT xL = (*x_)(i, 0), yL = y[i], xR = (*x_)(i + 1, 0), yR = y[i + 1];
    RealT dydx = (yR - yL) / (xR - xL);
    return yL + dydx * (x_pt - xL);
  }

  Adv_Diff_Constraint(int m)
  {
    m_ = m;
    RealT h = 1.0 / static_cast<RealT>(m_ - 1);
    x_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      x_->Set_Entry(k, 0, static_cast<RealT>(k) / static_cast<RealT>(m_ - 1));
    }

    S_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    M_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);
    robin_bc_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, m_);

    S_->Set_Entry(0, 0, 1.0 / h);
    S_->Set_Entry(0, 1, -1.0 / h);
    for (int i = 1; i < m_ - 1; i++)
    {
      S_->Set_Entry(i, i, 2.0 / h);
      S_->Set_Entry(i, i - 1, -1.0 / h);
      S_->Set_Entry(i, i + 1, -1.0 / h);
    }
    S_->Set_Entry(m_ - 1, m_ - 2, -1.0 / h);
    S_->Set_Entry(m_ - 1, m_ - 1, 1.0 / h);

    M_->Set_Entry(0, 0, (1.0 / 3.0) * h);
    M_->Set_Entry(0, 1, (1.0 / 6.0) * h);
    for (int i = 1; i < m_ - 1; i++)
    {
      M_->Set_Entry(i, i, (2.0 / 3.0) * h);
      M_->Set_Entry(i, i - 1, (1.0 / 6.0) * h);
      M_->Set_Entry(i, i + 1, (1.0 / 6.0) * h);
    }
    M_->Set_Entry(m_ - 1, m_ - 2, (1.0 / 6.0) * h);
    M_->Set_Entry(m_ - 1, m_ - 1, (1.0 / 3.0) * h);

    robin_bc_->Set_Entry(0, 0, 1.0);
    robin_bc_->Set_Entry(m_ - 1, m_ - 1, 1.0);

    alpha_ = 1.0;

    forcing_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT>> tmp = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(m_, 1);
    for (int k = 0; k < m_; k++)
    {
      RealT val = std::exp(-50.0 * std::pow((*x_)(k, 0) - 0.5, 2.0));
      tmp->Set_Entry(k, 0, val);
    }
    M_->Multiply(*forcing_, *tmp);
  }

  virtual ~Adv_Diff_Constraint()
  {
  }
};

#endif
