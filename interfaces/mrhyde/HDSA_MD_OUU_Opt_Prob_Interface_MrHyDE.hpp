#ifndef HDSA_MD_OUU_OPT_PROB_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OUU_OPT_PROB_INTERFACE_MRHYDE_HPP

template <class RealT>
class MD_OUU_Opt_Prob_Interface_MrHyDE : public HDSA::MD_OUU_Opt_Prob_Interface<RealT>
{

private:
  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> opt_prob_interface_mrhyde_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  HDSA::Ptr<ROL::SampleGenerator<RealT>> sampler_;
  int ens_size_;

public:
  MD_OUU_Opt_Prob_Interface_MrHyDE(HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> &opt_prob_interface_mrhyde, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, HDSA::Ptr<ROL::SampleGenerator<RealT>> &sampler, int ens_size) : HDSA::MD_OUU_Opt_Prob_Interface<RealT>(ens_size), opt_prob_interface_mrhyde_(opt_prob_interface_mrhyde), params_(params), sampler_(sampler), ens_size_(ens_size)
  {
  }

  virtual ~MD_OUU_Opt_Prob_Interface_MrHyDE()
  {
  }

  void Apply_Solution_Operator_z_Jacobian_Transpose_Per_Sample(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &z, int s) const
  {
    std::vector<RealT> pt_s = sampler_->getMyPoint(s);
    params_->updateParams(pt_s, "stochastic");

    HDSA::Ptr<ROL::Vector<RealT>> z_tmp = params_->getCurrentVector().clone();
    z_tmp->zero();
    MrHyDE_OptVector z_mrhyde = Teuchos::dyn_cast<MrHyDE_OptVector>(const_cast<ROL::Vector<RealT> &>(*z_tmp));
    params_->updateParams(z_mrhyde);

    opt_prob_interface_mrhyde_->Apply_Solution_Operator_z_Jacobian_Transpose(z_out, u_in, z);
  }

  void Apply_RS_Hessian_Per_Sample(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, int s) const
  {
    std::vector<RealT> pt_s = sampler_->getMyPoint(s);
    params_->updateParams(pt_s, "stochastic");

    HDSA::Ptr<ROL::Vector<RealT>> z_tmp = params_->getCurrentVector().clone();
    z_tmp->zero();
    MrHyDE_OptVector z_mrhyde = Teuchos::dyn_cast<MrHyDE_OptVector>(const_cast<ROL::Vector<RealT> &>(*z_tmp));
    params_->updateParams(z_mrhyde);

    opt_prob_interface_mrhyde_->Apply_RS_Hessian(z_out, z_in, z);
  }

  void Misfit_Gradient_Per_Sample(HDSA::Vector<RealT> &u_grad, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, int s) const
  {
    opt_prob_interface_mrhyde_->Misfit_Gradient(u_grad, u, z);
  }

  void Apply_Misfit_Hessian_Per_Sample(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z, int s) const
  {
    opt_prob_interface_mrhyde_->Apply_Misfit_Hessian(u_out, u_in, u, z);
  }
};

#endif
