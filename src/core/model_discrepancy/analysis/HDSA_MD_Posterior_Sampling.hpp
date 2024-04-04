#ifndef HDSA_MD_POSTERIOR_SAMPLING_HPP
#define HDSA_MD_POSTERIOR_SAMPLING_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Posterior_Sampling{

  private:
    HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
    HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface_;
    HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;
  
  public:
    HDSA::Ptr<HDSA::MD_Posterior_Data<RealT> > post_data;

    MD_Posterior_Sampling(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface)
    { 
      data_interface_ = data_interface;
      u_prior_interface_ = u_prior_interface;
      z_prior_interface_ = z_prior_interface;
      post_data = HDSA::makePtr<HDSA::MD_Posterior_Data<RealT> >();
    }

    ~MD_Posterior_Sampling(void)
    { }


    void Compute_Posterior_Data(RealT & alpha_d,int & num_samples)
    {
      post_data->Compute_Posterior_Data(*data_interface_,*u_prior_interface_,*z_prior_interface_,alpha_d,num_samples);
    }

    HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT> > Posterior_Discrepancy_Samples(std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > & z) const
    {
      int p = z.size();
      HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT> > delta = HDSA::makePtr<HDSA::MD_Posterior_Vectors<RealT> >(p,*z[0]);
      
      HDSA::Ptr<HDSA::MultiVector<RealT> > Zc = HDSA::makePtr<HDSA::MultiVector<RealT> >(p-1,*z[0]);
      HDSA::Ptr<HDSA::MultiVector<RealT> > W_z_inv_Zc =HDSA::makePtr<HDSA::MultiVector<RealT> >(p-1,*z[0]);
      for(int k = 0; k < p-1; k++)
	{
	  (*Zc)[k]->set(*(*post_data->Z)[k+1]);
	  (*Zc)[k]->axpy(-1.0,*data_interface_->z_opt);

	  (*W_z_inv_Zc)[k]->set(*(*post_data->W_z_inv_Z)[k+1]);
	  (*W_z_inv_Zc)[k]->axpy(-1.0,*post_data->W_z_inv_z_opt);
	}

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Zc_W_z_inv_Zc = Zc->MatMat(*W_z_inv_Zc);

      for(int k = 0; k < p; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > dz = z[k]->clone();
	  dz->set(*z[k]);
	  dz->axpy(-1.0,*post_data->W_z_inv_z_opt);

	  // Continue implementation at line 85 of the Sabl version


	}


      return delta;
    }

  };

}

#endif
