#ifndef HDSA_MD_PRIOR_SAMPLING_HPP
#define HDSA_MD_PRIOR_SAMPLING_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Prior_Sampling {

  private:
    const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
    const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface_;
    const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;
    int z_pert_subsample_factor_;

  public:
    MD_Prior_Sampling(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface):
     data_interface_(data_interface), u_prior_interface_(u_prior_interface), z_prior_interface_(z_prior_interface)
    { 
      z_pert_subsample_factor_ = 1;
    }

    virtual ~MD_Prior_Sampling()
    { }

    void Generate_Prior_Discrepancy_Sample_Data(int & num_samps)
    {
      Generate_Prior_Discrepancy_z_opt_Sample_Data(num_samps);
      Generate_Prior_Discrepancy_z_pert_Sample_Data();
    }

    void Generate_Prior_Discrepancy_z_opt_Sample_Data(int & num_samps)
    {
      std::cout << "Need to implement Generate_Prior_Discrepancy_z_opt_Sample_Data" << std::endl;
    }

    void Generate_Prior_Discrepancy_z_pert_Sample_Data(void)
    {
      std::cout << "Need to implement Generate_Prior_Discrepancy_z_pert_Sample_Data" << std::endl;
    }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Prior_Discrepancy_Samples_at_z_opt(const int & num_samples)
    {
      HDSA::Ptr<HDSA::MultiVector<RealT> > delta_samples = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_samples,*data_interface_->get_u_opt());
      u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*delta_samples);
      return delta_samples;
    }

    std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > Prior_Discrepancy_Samples(const HDSA::MultiVector<RealT> & z, const int & num_samples)
    {
      int N = z.Number_of_Vectors();
      
      std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > delta_samples;
      delta_samples.resize(num_samples);

      // Compute Sigma = Z^T*M_z*W_z_inv*M_z*Z
      HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::MultiVector<RealT> > M_z_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::MultiVector<RealT> > W_z_inv_M_z_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->get_z_opt());
      for(int k = 0; k < N; k++)
        {
          HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = (*Z)[k];
          z_tmp1->set(*z[k]);
          z_tmp1->axpy(-1.0,*data_interface_->get_z_opt());

          HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = (*M_z_Z)[k];
          z_prior_interface_->Apply_M_z(*z_tmp2,*z_tmp1);

          HDSA::Ptr<HDSA::Vector<RealT> > z_tmp3 = (*W_z_inv_M_z_Z)[k];
          z_prior_interface_->Apply_W_z_Inverse(*z_tmp3,*z_tmp2);
        }
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Sigma = M_z_Z->MatMat(*W_z_inv_M_z_Z);

      // Factorize Sigma = R*R^T
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(N,N);
      HDSA::Linear_Algebra::Cholesky_Factorization(*Sigma,*R);

      // Loop over samples
      for(int i = 0; i < num_samples; i++)
      {
        delta_samples[i] = HDSA::makePtr<HDSA::MultiVector<RealT> >(N,*data_interface_->get_u_opt());

        // Generate random sample i 
        HDSA::Ptr<HDSA::MultiVector<RealT> > u_rand = HDSA::makePtr<HDSA::MultiVector<RealT> >(N+1,*data_interface_->get_u_opt());
        u_prior_interface_->Sample_with_Covariance_W_u_Inverse(*u_rand);

        for(int k = 0; k < N; k++)
          {
            HDSA::Ptr<HDSA::Vector<RealT> > u_out = (*delta_samples[i])[k];
            u_out->set(*(*u_rand)[N]);
            for(int j = 0; j < k+1; j++)
              {
                u_out->axpy((*R)(j,k),*(*u_rand)[j]);
              }
          }

      }

      return delta_samples;
    }

    void Compute_Temporal_Data(void)
    {
      std::cout << "Need to implement Compute_Temporal_Data" << std::endl;
    }

    void Compute_z_pert_Data(void)
    {
      std::cout << "Need to implement Compute_z_pert_Data" << std::endl;
    }

    void Compute_Delta_z_opt_Metrics(void)
    {
      std::cout << "Need to implement Compute_Delta_z_opt_Metrics" << std::endl;
    }

    void Compute_Delta_z_pert_Metrics(void)
    {
      std::cout << "Need to implement Compute_Delta_z_pert_Metrics" << std::endl;
    }

  };

}

#endif
