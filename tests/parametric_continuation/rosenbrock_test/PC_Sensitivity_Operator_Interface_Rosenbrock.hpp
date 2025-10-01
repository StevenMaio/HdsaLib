#ifndef HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_ROSENBROCK_HPP
#define HDSA_PC_SENSITIVITY_OPERATOR_INTERFACE_ROSENBROCK_HPP

#include "HDSA_Std_Vector.hpp"

template <class RealT>
class PC_Sensitivity_Operator_Interface_Rosenbrock : public HDSA::PC_Sensitivity_Operator_Interface<RealT>
{

private:
  HDSA::Ptr<Rosenbrock<RealT>> rosenbrock_;

public:
  PC_Sensitivity_Operator_Interface_Rosenbrock(HDSA::Ptr<Rosenbrock<RealT>> &rosenbrock) : rosenbrock_(rosenbrock)
  {
  }

  virtual ~PC_Sensitivity_Operator_Interface_Rosenbrock()
  {
  }

  void Gradient(HDSA::Vector<RealT> &grad, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    const HDSA::Std_Vector<RealT> &z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);
    HDSA::Std_Vector<RealT> &grad_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(grad);

    int d = grad.dimension();
    RealT val = 0.0;
    std::vector<RealT> grad_tmp = std::vector<RealT>(d);
    std::vector<std::vector<RealT>> hess_tmp;
    hess_tmp.resize(d);
    for (int i = 0; i < d; i++)
    {
      hess_tmp[i].resize(d);
    }
    rosenbrock_->J(val, grad_tmp, hess_tmp, *z_std.get_std_vec(), *theta_std.get_std_vec());

    for (int i = 0; i < d; i++)
    {
      grad_std.Set_Entry(i, grad_tmp[i]);
    }
  }

  void Apply_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    const HDSA::Std_Vector<RealT> &z_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z);
    const HDSA::Std_Vector<RealT> &theta_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta);
    const HDSA::Std_Vector<RealT> &z_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(z_in);
    HDSA::Std_Vector<RealT> &z_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(z_out);

    int d = z_in.dimension();
    RealT val = 0.0;
    std::vector<RealT> grad_tmp = std::vector<RealT>(d);
    std::vector<std::vector<RealT>> hess_tmp;
    hess_tmp.resize(d);
    for (int i = 0; i < d; i++)
    {
      hess_tmp[i].resize(d);
    }
    rosenbrock_->J(val, grad_tmp, hess_tmp, *z_std.get_std_vec(), *theta_std.get_std_vec());

    for (int i = 0; i < d; i++)
    {
      RealT tmp = 0.0;
      for (int j = 0; j < d; j++)
      {
        tmp += hess_tmp[i][j] * z_in_std(j);
      }
      z_out_std.Set_Entry(i, tmp);
    }
  }

  void Apply_B(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &theta_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta) const
  {
    const HDSA::Std_Vector<RealT> &theta_in_std = dynamic_cast<const HDSA::Std_Vector<RealT> &>(theta_in);
    HDSA::Std_Vector<RealT> &z_out_std = dynamic_cast<HDSA::Std_Vector<RealT> &>(z_out);

    std::vector<std::vector<RealT>> B = rosenbrock_->Compute_B();

    int d = z_out.dimension();
    for (int i = 0; i < d; i++)
    {
      RealT tmp = 0.0;
      for (int j = 0; j < d - 1; j++)
      {
        tmp += B[i][j] * theta_in_std(j);
      }
      z_out_std.Set_Entry(i, tmp);
    }
  }
};

#endif
