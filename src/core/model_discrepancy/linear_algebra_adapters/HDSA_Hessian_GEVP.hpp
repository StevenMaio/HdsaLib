#ifndef HDSA_HESSIAN_GEVP_HPP
#define HDSA_HESSIAN_GEVP_HPP

namespace HDSA
{

  template <class RealT>
  class Hessian_GEVP : public HDSA::Randomized_GEVP<RealT> {

  private:
    HDSA::Ptr<HDSA::Vector<RealT> > z_;
    HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > opt_prob_interface_;
    HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;
    RealT normalization_coeff_;

  public:

    Hessian_GEVP(const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > & opt_prob_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface, const HDSA::Vector<RealT> & z): HDSA::Randomized_GEVP<RealT>(z)
    {
      z_ = z.clone();
      z_->set(z);
      opt_prob_interface_ = opt_prob_interface;
      z_prior_interface_ = z_prior_interface;
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z.clone();
      z_prior_interface_->Apply_W_z(*z_tmp,*z_);
      normalization_coeff_ = z_->dot(*z_tmp);
    }

    virtual ~Hessian_GEVP()
    { }

    void Compute_Hessian_GEVP(HDSA::MultiVector<RealT> & evecs, HDSA::Dense_Matrix<RealT> & evals, int & num_evals, int & oversampling)
    {
      Compute_GEVP(evecs, evals, num_evals, oversampling);
    }

    void Apply_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
    {
      opt_prob_interface_->Apply_RS_Hessian(vec_out,vec_in,*z_);
      vec_out.set(vec_in);
    }

    void Apply_Weighting_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
    {
      z_prior_interface_->Apply_W_z(vec_out,vec_in);
      vec_out.scale(1.0/normalization_coeff_);
    }

    void Apply_Weighting_Operator_Inverse(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const
    {
      z_prior_interface_->Apply_W_z_Inverse(vec_out,vec_in);
      vec_out.scale(normalization_coeff_);
    }
     
    void Apply_Weighting_Operator_Preconditioner_Factor(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in)
    {
      z_prior_interface_->Apply_W_z_Inverse_Factor(vec_out,vec_in);
      vec_out.scale(std::sqrt(normalization_coeff_));
    }

  };
    
}
  


#endif
