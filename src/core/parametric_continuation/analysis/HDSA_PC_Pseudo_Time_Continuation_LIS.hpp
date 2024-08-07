#ifndef HDSA_PC_PSEUDO_TIME_CONTINUATION_LIS_HPP
#define HDSA_PC_PSEUDO_TIME_CONTINUATION_LIS_HPP

namespace HDSA
{

  template <class RealT>
  class PC_Pseudo_Time_Continuation_LIS : public HDSA::PC_Pseudo_Time_Continuation<RealT> {

  private:
    HDSA::Ptr<HDSA::PC_LIS_Interface<RealT> > lis_interface_;
    HDSA::Ptr<HDSA::Vector<RealT> > z_bar_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_bar_;
    int rank_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > evecs_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > evals_;
    
  public:
    
    PC_Pseudo_Time_Continuation_LIS(const HDSA::Ptr<HDSA::Vector<RealT> > & z_bar, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_bar, const HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT> > & sen_op, const HDSA::Ptr<HDSA::PC_LIS_Interface<RealT> > & lis_interface): 
      HDSA::PC_Pseudo_Time_Continuation<RealT>(z_bar,theta_bar,sen_op), lis_interface_(lis_interface), z_bar_(z_bar), theta_bar_(theta_bar)
    {
      rank_ = 0;
    }

    virtual ~PC_Pseudo_Time_Continuation_LIS()
    { }

    void Compute_Hessian_GEVP(const HDSA::Vector<RealT> & z, const HDSA::Vector<RealT> & theta, const int & num_evals, const int & oversampling)
    {
      rank_ = num_evals;
      HDSA::Ptr<HDSA::Randomized_GEVP<RealT> > hessian_gevp = HDSA::makePtr<PC_LIS_Hessian_GEVP<RealT> >(lis_interface_,*z_bar_,*theta_bar_);
      evecs_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_evals,z_bar_);
      evals_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_evals,1);
      hessian_gevp->Compute_GEVP(*evecs_,*evals_,num_evals,oversampling);
      std::ofstream fout;
      std::string name = "hessian_evals.txt";
      fout.open(name);
      for(int i = 0; i < num_evals; i++)
	{
	  fout << std::setprecision(16) << (*evals_)(i,0) << "  ";
	}
      fout.close();
    }
    
  protected:

    // Overload this function if a better initialization is available
    void Apply_Initial_Inverse_BFGS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      lis_interface_->Apply_Prior_Covariance(z_out,z_in);
      if(rank_ > 0)
	{
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp = evecs_->MatVec(z_in);
	  for(int k = 0; k < rank_; k++)
	    {
	      RealT val = (*tmp)(k,0)*(*evals_)(k,0)/(1.0 + (*evals_)(k,0));
	      z_out.axpy(-val,*(*evecs_)[k]);
	    }
	}
    }

    template <class ScalarType>
    class PC_LIS_Hessian_GEVP : public HDSA::Randomized_GEVP<ScalarType> {

    private:
      HDSA::Ptr<HDSA::Vector<ScalarType> > z_bar_;
      HDSA::Ptr<HDSA::Vector<ScalarType> > theta_bar_;
      HDSA::Ptr<HDSA::PC_LIS_Interface<ScalarType> > lis_interface_;
      
    public:
      PC_LIS_Hessian_GEVP(const HDSA::Ptr<HDSA::PC_LIS_Interface<ScalarType> > & lis_interface, const HDSA::Vector<ScalarType> & z_bar, const HDSA::Vector<ScalarType> & theta_bar):
	HDSA::Randomized_GEVP<ScalarType>(z_bar)
      {
	z_bar_ = z_bar.clone();
	z_bar_->set(z_bar);
	theta_bar_ = theta_bar.clone();
	theta_bar_->set(theta_bar);
	lis_interface_ = lis_interface;
      }
      
      virtual ~PC_LIS_Hessian_GEVP()
      { }

      void Apply_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	lis_interface_->Apply_Misfit_Hessian(vec_out, vec_in, *z_bar_, *theta_bar_);
      }

      void Apply_Weighting_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	lis_interface_->Apply_Prior_Precision(vec_out, vec_in);
      }

      void Apply_Weighting_Operator_Inverse(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	lis_interface_->Apply_Prior_Covariance(vec_out, vec_in);
      }

      void Generate_Random_Samples(HDSA::MultiVector<RealT> & samples) const 
      {
	lis_interface_->Generate_Prior_Samples(samples);
      }

    };
    
  };

}

#endif
