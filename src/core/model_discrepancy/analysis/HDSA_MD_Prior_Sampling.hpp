#ifndef HDSA_MD_PRIOR_SAMPLING_HPP
#define HDSA_MD_PRIOR_SAMPLING_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Prior_Sampling {

  private:
    HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
    HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface_;
    HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;

  public:
    MD_Prior_Sampling(HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface, HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface): data_interface_(data_interface), u_prior_interface_(u_prior_interface), z_prior_interface_(z_prior_interface)
    { }

    virtual ~MD_Prior_Sampling()
    { }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Prior_Discrepancy_Samples_at_z_opt(const int & num_samples)
    {
      HDSA::Ptr<HDSA::MultiVector<RealT> > delta_samples = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_samples,*data_interface_->u_opt);
      u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*delta_samples);
      return delta_samples;
    }

    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Prior_Discrepancy_Samples(const HDSA::MultiVector<RealT> & z, const int & num_samples)
    {
      int N = z.Number_of_Vectors();
      
      std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > delta_samples;
      delta_samples.resize(num_samples);

      // Compute Sigma = Z^T*W_z_inv*Z
      HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->z_opt);
      HDSA::Ptr<HDSA::MultiVector<RealT> > W_z_inv_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->z_opt);
      for(int k = 0; k < N; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = (*Z)[k];
	  z_tmp1->set(*z[k]);
	  z_tmp1->axpy(-1.0,*data_interface_->z_opt);

	  HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = (*W_z_inv_Z)[k];
	  z_prior_interface_->Apply_W_z_Inverse(*z_tmp2,*z_tmp1);
	}
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Sigma = Z->MatMat(*W_z_inv_Z);

      // Factorize Sigma = R*R^T
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,N);
      HDSA::Linear_Algebra::Cholesky_Factorization(*Sigma,*R);

      // Loop over samples
      for(int i = 0; i < num_samples; i++)
	{
	 delta_samples[i] = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->u_opt);

	 // Generate random sample i 
	 HDSA::Ptr<HDSA::MultiVector<RealT> > u_rand = HDSA::makePtr<HDSA::MultiVector<RealT> >(N+1,*data_interface_->u_opt);
	 u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*u_rand);

	 for(int k = 0; k < N; k++)
	   {
	     HDSA::Ptr<HDSA::Vector<RealT> > u_out = (*delta_samples[i])[k];
	     u_out->set(*(*u_rand)[N]);
	     for(int j = 0; j < k; j++)
	       {
		 u_out->axpy((*R)(j,k),*(*u_rand)[j]);
	       }
	   }

	}

      return delta_samples;
    }

  };

}

#endif
