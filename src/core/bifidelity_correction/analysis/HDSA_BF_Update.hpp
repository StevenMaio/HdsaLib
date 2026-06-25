/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis

 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_BF_UPDATE_HPP
#define HDSA_BF_UPDATE_HPP

#include "HDSA_BF_Sol_Op_Interface.hpp"
#include "HDSA_MD_Opt_Prob_Interface.hpp"
#include "HDSA_Hessian_Inversion.hpp"

namespace HDSA
{

  template <class RealT>
  class BF_Update
  {

  private:
    HDSA::Ptr<HDSA::BF_Sol_Op_Interface<RealT>> sol_op_interface_;
    HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> opt_prob_interface_;
    HDSA::Ptr<HDSA::Hessian_Inversion<RealT>> hess_invert_;

  public:
    BF_Update(const HDSA::Ptr<HDSA::BF_Sol_Op_Interface<RealT>> &sol_op_interface, const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> &opt_prob_interface, int verbosity = 0, std::string hessian_solver = "CG") : sol_op_interface_(sol_op_interface), opt_prob_interface_(opt_prob_interface)
    {
      hess_invert_ = HDSA::makePtr<BF_Hessian_Inversion<RealT>>(opt_prob_interface, verbosity, hessian_solver);
    }

    ~BF_Update(void)
    {
    }

    HDSA::Ptr<HDSA::Vector<RealT>> Update(const HDSA::Vector<RealT> &u_lofi, const HDSA::Vector<RealT> &z_lofi) const
    {
      HDSA::Ptr<HDSA::Vector<RealT>> B = z_lofi.Clone();
      HDSA::Ptr<HDSA::Vector<RealT>> z_update = z_lofi.Clone();

      HDSA::Ptr<HDSA::Vector<RealT>> u_hifi = u_lofi.Clone();
      sol_op_interface_->State_Solve(*u_hifi, z_lofi);

      HDSA::Ptr<HDSA::Vector<RealT>> J_grad_u = u_lofi.Clone();
      opt_prob_interface_->Misfit_Gradient(*J_grad_u, u_lofi, z_lofi);

      sol_op_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*B, *J_grad_u, z_lofi);

      HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z_lofi.Clone();
      opt_prob_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*z_tmp, *J_grad_u, z_lofi);
      B->Scaled_Plus(-1.0, *z_tmp);

      HDSA::Ptr<HDSA::Vector<RealT>> discrep = u_lofi.Clone();
      discrep->Set(*u_hifi);
      discrep->Scaled_Plus(-1.0, u_lofi);
      HDSA::Ptr<HDSA::Vector<RealT>> u_tmp = u_lofi.Clone();
      opt_prob_interface_->Apply_Misfit_Hessian(*u_tmp, *discrep, u_lofi, z_lofi);
      z_tmp->Zeros();
      opt_prob_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*z_tmp, *u_tmp, z_lofi);
      B->Plus(*z_tmp);

      B->Scale(-1.0);
      hess_invert_->Apply_RS_Hessian_Inverse(*z_update, *B, z_lofi);
      z_update->Plus(z_lofi);
      return z_update;
    }

    template <class ScalarType>
    class BF_Hessian_Inversion : public HDSA::Hessian_Inversion<ScalarType>
    {
    private:
      HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType>> opt_prob_interface_;

    public:
      BF_Hessian_Inversion(const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType>> &opt_prob_interface, int verbosity, std::string solver) : HDSA::Hessian_Inversion<ScalarType>(verbosity, 1.e-6, solver), opt_prob_interface_(opt_prob_interface)
      {
      }

      ~BF_Hessian_Inversion()
      {
      }

      void Apply_RS_Hessian(HDSA::Vector<ScalarType> &z_out, const HDSA::Vector<ScalarType> &z_in, const HDSA::Vector<ScalarType> &z) const
      {
        opt_prob_interface_->Apply_RS_Hessian(z_out, z_in, z);
      }
    };
  };

}

#endif
