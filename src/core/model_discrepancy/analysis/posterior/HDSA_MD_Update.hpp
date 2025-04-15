#ifndef HDSA_MD_UPDATE_HPP
#define HDSA_MD_UPDATE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Update{

  private:
    HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface_;
    HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface_;
    HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface_;
    HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > opt_prob_interface_;
    HDSA::Ptr<HDSA::MD_Posterior_Sampling<RealT> > post_sampling_;
    HDSA::Ptr<HDSA::MD_Hessian_Analysis<RealT> > hessian_analysis_;
    HDSA::Ptr<HDSA::Vector<RealT> > state_grad_;
    RealT state_grad_W_u_inv_state_grad_;

  public:

    MD_Update(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > & data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > & u_prior_interface, 
	      const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > & z_prior_interface, const HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > & opt_prob_interface, 
	      const HDSA::Ptr<HDSA::MD_Posterior_Sampling<RealT> > & post_sampling, const HDSA::Ptr<HDSA::MD_Hessian_Analysis<RealT> > & hessian_analysis): 
      data_interface_(data_interface), u_prior_interface_(u_prior_interface), z_prior_interface_(z_prior_interface), opt_prob_interface_(opt_prob_interface), 
      post_sampling_(post_sampling), hessian_analysis_(hessian_analysis)
    { 
      state_grad_ = data_interface_->get_u_opt()->clone();
      opt_prob_interface_->Misfit_Gradient(*state_grad_,*data_interface_->get_u_opt(),*data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::Vector<RealT> > u_tmp = state_grad_->clone();
      u_prior_interface_->Apply_W_u_Inverse(*u_tmp,*state_grad_);
      state_grad_W_u_inv_state_grad_ = u_tmp->dot(*state_grad_);
    }

    ~MD_Update(void)
    { }


    HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT> > Posterior_Update_Samples(void) const
    {
      HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT> > posterior_samples = HDSA::makePtr<HDSA::MD_Posterior_Vectors<RealT> >(post_sampling_->post_data->num_samples,*data_interface_->get_z_opt());
      HDSA::Ptr<HDSA::Vector<RealT> > z_update_mean = Posterior_Update_Mean();
      posterior_samples->mean->set(*z_update_mean);

      //////////////////// B_theta_hat
      HDSA::Ptr<HDSA::MultiVector<RealT> > u_tmp1 =  HDSA::makePtr<HDSA::MultiVector<RealT> >(post_sampling_->post_data->num_samples,*data_interface_->get_u_opt());
      HDSA::Ptr<HDSA::MultiVector<RealT> > B_theta_hat =  HDSA::makePtr<HDSA::MultiVector<RealT> >(post_sampling_->post_data->num_samples,*data_interface_->get_z_opt());
      for(int i = 0; i < post_sampling_->post_data->N; i++)
	{
	  RealT coeff1 = post_sampling_->post_data->sum_g_vecs[i]/std::sqrt((*post_sampling_->post_data->Mu)(i,0));
	  u_tmp1->axpy(coeff1,*post_sampling_->post_data->u_i_hat[i]);

	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > coeff2 = post_sampling_->post_data->u_i_hat[i]->MatVec(*state_grad_);

	  // Compute W_z_inv_yi                                                                                                                                                                                                               
	  HDSA::Ptr<HDSA::Vector<RealT> > W_z_inv_yi = data_interface_->get_z_opt()->clone();
	  W_z_inv_yi->axpy(-post_sampling_->post_data->sum_g_vecs[i],*post_sampling_->post_data->W_z_inv_z_opt);
	  for(int j = 0; j < post_sampling_->post_data->N; j++)
	    {
	      W_z_inv_yi->axpy((*post_sampling_->post_data->g_vecs)(j,i),*(*post_sampling_->post_data->W_z_inv_Z)[j]);
	    }

	  for(int k = 0; k < post_sampling_->post_data->num_samples; k++)
	    {
	      RealT val = (*coeff2)(k,0)/std::sqrt((*post_sampling_->post_data->Mu)(i,0));
	      (*B_theta_hat)[k]->axpy(val,*W_z_inv_yi);
	    }
	}

      for(int k = 0; k < post_sampling_->post_data->num_samples; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > u_tmp2 = data_interface_->get_u_opt()->clone();
	  opt_prob_interface_->Apply_Misfit_Hessian(*u_tmp2,*(*u_tmp1)[k],*data_interface_->get_u_opt(),*data_interface_->get_z_opt());
	  
	  HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = data_interface_->get_z_opt()->clone();
	  opt_prob_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*z_tmp1,*u_tmp2,*data_interface_->get_z_opt());
	
	  (*B_theta_hat)[k]->plus(*z_tmp1);
	}

      B_theta_hat->scale(std::sqrt(post_sampling_->post_data->alpha_d));


      //////////////// B_theta_breve
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp1 = post_sampling_->post_data->z_breve->MatMat(*post_sampling_->post_data->Zc); // Dimension (N-1)x(num_post_samples)
   
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > tmp2 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(post_sampling_->post_data->N-1,post_sampling_->post_data->num_samples);
      HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*post_sampling_->post_data->Zc_W_z_inv_Zc,*tmp2,*tmp1);

      HDSA::Ptr<HDSA::MultiVector<RealT> > B_theta_breve =  HDSA::makePtr<HDSA::MultiVector<RealT> >(post_sampling_->post_data->num_samples,*data_interface_->get_z_opt());
      for(int k = 0; k < post_sampling_->post_data->num_samples; k++)
	{
	  for(int i = 0; i < post_sampling_->post_data->N-1; i++)
	    {
	      (*B_theta_breve)[k]->axpy(-(*tmp2)(i,k),*(*post_sampling_->post_data->W_z_inv_Zc)[i]);
	    }
	}
      B_theta_breve->axpy(1.0,*post_sampling_->post_data->z_breve);
      B_theta_breve->scale(std::sqrt(state_grad_W_u_inv_state_grad_));

      HDSA::Ptr<HDSA::MultiVector<RealT> > z_update_samples = posterior_samples->samples;
      for(int k = 0; k < post_sampling_->post_data->num_samples; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = data_interface_->get_z_opt()->clone();
	  z_tmp2->set(*(*B_theta_hat)[k]);
	  z_tmp2->plus(*(*B_theta_breve)[k]);
	  hessian_analysis_->Apply_RS_Hessian_Inverse(*(*z_update_samples)[k],*z_tmp2,*data_interface_->get_z_opt());
	  (*z_update_samples)[k]->scale(-1.0);
	  (*z_update_samples)[k]->plus(*z_update_mean);
	}

      return posterior_samples;
    }

    HDSA::Ptr<HDSA::Vector<RealT> > Posterior_Update_Mean(void) const
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_update_mean = data_interface_->get_z_opt()->clone();

      HDSA::Ptr<HDSA::Vector<RealT> > u_tmp1 = state_grad_->clone();
      for(int ell = 0; ell < post_sampling_->post_data->N; ell++)
	{
	  u_tmp1->plus(*(*post_sampling_->post_data->u_ell)[ell]);
	  for(int i = 0; i < post_sampling_->post_data->N; i++)
	    {
	      RealT coeff = post_sampling_->post_data->sum_g_vecs[i]*(*post_sampling_->post_data->b_i_ell)(i,ell); 
	      u_tmp1->axpy(-coeff,*(*post_sampling_->post_data->u_i_ell[i])[ell]);
	    }
	}

      HDSA::Ptr<HDSA::Vector<RealT> > u_tmp2 = u_tmp1->clone();
      opt_prob_interface_->Apply_Misfit_Hessian(*u_tmp2,*u_tmp1,*data_interface_->get_u_opt(),*data_interface_->get_z_opt());

      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = z_update_mean->clone();
      opt_prob_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*z_tmp1,*u_tmp2,*data_interface_->get_z_opt());
							
      for(int ell = 0; ell < post_sampling_->post_data->N; ell++)
	{
	  RealT coeff1 = state_grad_->dot(*(*post_sampling_->post_data->u_ell)[ell]);
	  z_tmp1->axpy(coeff1,*(*post_sampling_->post_data->W_z_inv_Z)[ell]);
	  z_tmp1->axpy(-coeff1,*post_sampling_->post_data->W_z_inv_z_opt);
	  for(int i = 0; i < post_sampling_->post_data->N; i++)
	    {
	      RealT coeff2 = state_grad_->dot(*(*post_sampling_->post_data->u_i_ell[i])[ell]);
	      coeff2 *= (*post_sampling_->post_data->b_i_ell)(i,ell); 

              z_tmp1->axpy(coeff2*post_sampling_->post_data->sum_g_vecs[i],*post_sampling_->post_data->W_z_inv_z_opt);
              for(int j = 0; j < post_sampling_->post_data->N; j++)
                {
                  z_tmp1->axpy(-coeff2*(*post_sampling_->post_data->g_vecs)(j,i),*(*post_sampling_->post_data->W_z_inv_Z)[j]);
		}

	    }
	}

      z_tmp1->scale(-1.0/post_sampling_->post_data->alpha_d);
      hessian_analysis_->Apply_RS_Hessian_Inverse(*z_update_mean,*z_tmp1,*data_interface_->get_z_opt());
      z_update_mean->plus(*data_interface_->get_z_opt());

      return z_update_mean;
    }

  };

}

#endif
