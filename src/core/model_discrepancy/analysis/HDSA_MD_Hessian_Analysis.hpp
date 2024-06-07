#ifndef HDSA_MD_HESSIAN_ANALYSIS_HPP
#define HDSA_MD_HESSIAN_ANALYSIS_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Hessian_Analysis {

  private:
    HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > opt_prob_interface_;
    HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;
    HDSA::Ptr<HDSA::Vector<RealT> > z_current_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > evecs_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > evals_;
    bool use_projector_;
    HDSA::Ptr<HDSA::Hessian_Inversion<RealT> > hess_invert_;

  public:
    MD_Hessian_Analysis(const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > & opt_prob_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface): opt_prob_interface_(opt_prob_interface), z_prior_interface_(z_prior_interface)
    { 
      use_projector_ = false;
      hess_invert_ = HDSA::makePtr<MD_Hessian_Inversion<RealT> >(opt_prob_interface);
    }

    virtual ~MD_Hessian_Analysis()
    { }

    void Compute_Hessian_GEVP(const HDSA::Vector<RealT> & z, const int & num_evals, const int & oversampling)
    {
      HDSA::Ptr<HDSA::Randomized_GEVP<RealT> > hessian_gevp = HDSA::makePtr<MD_Hessian_GEVP<RealT> >(opt_prob_interface_, z_prior_interface_, z);
      evecs_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_evals,z);
      evals_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_evals,1);
      z_current_ = z.clone();
      z_current_->set(z);
      hessian_gevp->Compute_GEVP(*evecs_,*evals_,num_evals,oversampling);
      use_projector_ = true;
      std::ofstream fout;
      std::string name = "hessian_evals.txt";
      fout.open(name);
      for(int i = 0; i < num_evals; i++)
	{
	  fout << std::setprecision(16) << (*evals_)(i,0) << "  ";
	}
      fout.close();
    }

    void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const
    {
      if(use_projector_)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z.clone();
	  z_tmp->set(z);
	  z_tmp->axpy(-1.0,*z_current_);
	  if(z_tmp->norm() > 0.0)
	    {
	      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
				       "The z input has changed. Need to recompute the GEVP." << std::endl);
	    }
	  else
	    {
	      Apply_Projected_RS_Hessian_Inverse(z_out,z_in);
	    }
	}
      else
	{
	  Apply_RS_Hessian_Inverse_Krylov(z_out,z_in,z);
	}

    }

    void Apply_Projected_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > coeffs = evecs_->MatVec(z_in);
      z_out.zeros();
      int r = evecs_->Number_of_Vectors();
      for(int k = 0; k < r; k++)
	{
	  RealT val = (*coeffs)(k,0)/(*evals_)(k,0);
	  z_out.axpy(val,*(*evecs_)[k]);
	}
    }

    void Apply_RS_Hessian_Inverse_Krylov(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const
    {
      hess_invert_->Apply_RS_Hessian_Inverse(z_out,z_in,z);
    }

    template <class ScalarType>
    class MD_Hessian_Inversion : public HDSA::Hessian_Inversion<ScalarType>
    {
    private:
      HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType> > opt_prob_interface_;

    public:
      MD_Hessian_Inversion(const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType> > opt_prob_interface): opt_prob_interface_(opt_prob_interface)
      { }

      ~MD_Hessian_Inversion()
      {}

      void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const
      {
	opt_prob_interface_->Apply_RS_Hessian(z_out,z_in,z);
      }

    };

    template <class ScalarType>
    class MD_Hessian_GEVP : public HDSA::Randomized_GEVP<ScalarType> {

    private:
      HDSA::Ptr<HDSA::Vector<ScalarType> > z_;
      HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType> > opt_prob_interface_;
      HDSA::Ptr<HDSA::MD_z_Prior_Interface<ScalarType> > z_prior_interface_;
      ScalarType normalization_coeff_;

    public:
      MD_Hessian_GEVP(const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<ScalarType> > & opt_prob_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<ScalarType> > & z_prior_interface,
		   const HDSA::Vector<ScalarType> & z):
	HDSA::Randomized_GEVP<ScalarType>(z)
      {
	z_ = z.clone();
	z_->set(z);
	opt_prob_interface_ = opt_prob_interface;
	z_prior_interface_ = z_prior_interface;
	HDSA::Ptr<HDSA::Vector<ScalarType> > z_tmp = z.clone();
	z_prior_interface_->Apply_W_z(*z_tmp,*z_);
	normalization_coeff_ = z_->dot(*z_tmp);
      }

      virtual ~MD_Hessian_GEVP()
      { }

      void Compute_Hessian_GEVP(HDSA::MultiVector<ScalarType> & evecs, HDSA::Dense_Matrix<ScalarType> & evals, int & num_evals, int & oversampling)
      {
	Compute_GEVP(evecs, evals, num_evals, oversampling);
      }

      void Apply_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	opt_prob_interface_->Apply_RS_Hessian(vec_out,vec_in,*z_);
      }

      void Apply_Weighting_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	z_prior_interface_->Apply_W_z(vec_out,vec_in);
	vec_out.scale(1.0/normalization_coeff_);
      }

      void Apply_Weighting_Operator_Inverse(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	z_prior_interface_->Apply_W_z_Inverse(vec_out,vec_in);
	vec_out.scale(normalization_coeff_);
      }

      void Generate_Random_Samples(HDSA::MultiVector<RealT> & samples) const 
      {
	z_prior_interface_->Sample_with_Covariance_W_z_Inverse(samples);
	if( samples[0]->norm() > 0.0 )
	  {
	    samples.scale(std::sqrt(normalization_coeff_));
	  }
	else
	  {
	    samples.randomize_standard_normal();
	  }
      }

    };

  };

}

#endif
