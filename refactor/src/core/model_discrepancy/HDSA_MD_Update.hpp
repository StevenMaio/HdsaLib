#ifndef HDSA_MD_UPDATE_HPP
#define HDSA_MD_UPDATE_HPP

namespace HDSA
{

template <class RealT>
class Model_Discrepancy_Update {

private:
  HDSA::Ptr<HDSA::Model_Discrepancy_Interface<RealT> > md_interface_;
  HDSA::Ptr<HDSA::Bayes_Posterior_Data<RealT> > post_data_;
  HDSA::Ptr<HDSA::Vector<RealT> > u_opt_;
  HDSA::Ptr<HDSA::Vector<RealT> > z_opt_;

public:

  Model_Discrepancy_Update(HDSA::Ptr<HDSA::Model_Discrepancy_Interface<RealT> > & md_interface)
  {   
    md_interface_ = md_interface;
    post_data_ = HDSA::makePtr<HDSA::Bayes_Posterior_Data<RealT> >();
    u_opt_ = md_interface->Load_Optimal_u();
    z_opt_ = md_interface->Load_Optimal_z();
  }

  virtual ~Model_Discrepancy_Update()
  { }

  void Compute_Posterior_Data(RealT alpha)
  {
    post_data_->Compute_Posterior_Data(*md_interface_,alpha,*u_opt_,*z_opt_);
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Posterior_Update_Mean() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_update = z_opt_->clone();

    std::vector<RealT> gi_sum = std::vector<RealT>(post_data_->N,0.0);
    for(int i = 0; i < post_data_->N; i++)
      {
	RealT val = 0.0;
	for(int k = 0; k < post_data_->N; k++)
	  {
	    val += (*post_data_->g_vecs)(k,i);
	  }
	gi_sum[i] = val;
      }
    
    HDSA::Ptr<HDSA::Vector<RealT> > B_theta_bar = z_opt_->clone();
    
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = u_opt_->clone();
    for(int ell = 0; ell < post_data_->N; ell++)
      {
	u_vec_1->plus(*(*post_data_->u_ell)[ell]);
	for(int i = 0; i < post_data_->N; i++)
	  {
	    RealT c = -1.0*gi_sum[i]*(*post_data_->b_i_ell)(i,ell);
	    u_vec_1->axpy(c,*(*post_data_->u_i_ell[i])[ell]);
	  }
      }

    HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = u_opt_->clone();
    md_interface_->Apply_Misfit_Hessian(*u_vec_2,*u_vec_1,*u_opt_,*z_opt_);
    md_interface_->Apply_Solution_Operator_z_Jacobian_Transpose(*B_theta_bar, *u_vec_2,*z_opt_);

    RealT c_run = 0.0;
    for(int ell = 0; ell < post_data_->N; ell++)
      {
	RealT c = (*post_data_->u_ell)[ell]->dot(*post_data_->state_grad);
	B_theta_bar->axpy(c,*(*post_data_->Gamma_inv_Z)[ell]);
	c_run += c;
      }
    B_theta_bar->axpy(-1.0*c_run,*post_data_->Gamma_inv_z_opt);
    
    HDSA::Ptr<HDSA::MultiVector<RealT> > Gamma_inv_w = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,*z_opt_);
    for(int i = 0; i < post_data_->N; i++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > gw = (*Gamma_inv_w)[i];
	for(int k = 0; k < post_data_->N; k++)
	  {
	    gw->axpy((*post_data_->g_vecs)(k,i),*(*post_data_->Gamma_inv_Z)[k]);
	    gw->axpy(-1.0*(*post_data_->g_vecs)(k,i),*post_data_->Gamma_inv_z_opt);
	  }
	for(int ell = 0; ell < post_data_->N; ell++)
	  {
	    RealT c = -1.0*( (*post_data_->b_i_ell)(i,ell) )*( (*post_data_->u_i_ell[i])[ell]->dot(*post_data_->state_grad) );
	    B_theta_bar->axpy(c,*gw);
	  }
      }
    
    B_theta_bar->scale(1.0/post_data_->alpha);
    md_interface_->Apply_RS_Hessian_Inverse(*z_update,*B_theta_bar,*z_opt_);
    z_update->scale(-1.0);
    z_update->plus(*z_opt_);
    return z_update;
  }


};

}

#endif


