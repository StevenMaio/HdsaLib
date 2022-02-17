#ifndef HDSA_BAYES_MODEL_ERROR_HPP
#define HDSA_BAYES_MODEL_ERROR_HPP

// This class executes the Bayesian HDSA with respect to model form error

namespace HDSA
{

  template <class RealT>
  class Bayes_Model_Error
  {
  private:
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    int num_prior_samps_;
    int num_post_samps_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z0_samps_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > U0_samps_;
    HDSA::Ptr<HDSA::Bayes_Posterior_Data<RealT> > post_data_;    

  public:
  
    Bayes_Model_Error(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
				const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
				const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    {
      num_prior_samps_ = parlist_sensitivity_->sublist("Bayes Model Error").get("Number of Prior Samples", 5);
      num_post_samps_ = parlist_sensitivity_->sublist("Bayes Model Error").get("Number of Posterior Samples", 5);
    }
    
    void Compute(void)
    {
     
      int num_procs = 1;
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,num_procs);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      HDSA::Opt_Problem_Objects_Bayes_Model_Error<RealT> &OP_Objects_subcomm_bayes_model_error = static_cast<HDSA::Opt_Problem_Objects_Bayes_Model_Error<RealT>&>(*OP_Objects_subcomm);
      HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > bayes_model_error_objects_subcomm = OP_Objects_subcomm_bayes_model_error.Get_Model_Error_Objects();
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects_subcomm->z->Clone();
      OP_Objects_subcomm->rs_obj->gradient_z(*grad_nominal,*OP_Objects_subcomm->z,*bayes_model_error_objects_subcomm->OP_Objects_->theta,true);
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());

      Compute_Prior_Samples(bayes_model_error_objects_subcomm,weight_matrices_subcomm,Sen_Op_subcomm);
      Compute_Posterior_Data(bayes_model_error_objects_subcomm,weight_matrices_subcomm,Sen_Op_subcomm);   
      Posterior_Sampling(bayes_model_error_objects_subcomm, Nom_subcomm, grad_nominal);
    }

    void Compute_Prior_Samples(const HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, 
			       const HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > & Sen_Op)
    {
      Z0_samps_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_prior_samps_,bayes_model_error_objects->OP_Objects_->z);
      U0_samps_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_prior_samps_,bayes_model_error_objects->OP_Objects_->u);
      HDSA::Ptr<HDSA::MultiVector<RealT> > z_prior_samps = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_prior_samps_,bayes_model_error_objects->OP_Objects_->z);
      HDSA::Ptr<HDSA::MultiVector<RealT> > delta_prior_samps = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_prior_samps_,bayes_model_error_objects->OP_Objects_->u);
      for(int k = 0; k < num_prior_samps_; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > z_vec = Sen_Op->Generate_Random_z_Vector();
	  HDSA::Ptr<HDSA::Vector<RealT> > zk = (*Z0_samps_)[k];	  
	  bayes_model_error_objects->Apply_Sqrt_Gamma_Mat_Inverse(zk,z_vec);

	  HDSA::Ptr<HDSA::Vector<RealT> > zpk = (*z_prior_samps)[k];
	  zpk->set(*zk);
	  zpk->plus(*bayes_model_error_objects->OP_Objects_->z);

	  HDSA::Ptr<HDSA::Vector<RealT> > u_vec = Sen_Op->Generate_Random_u_Vector();
	  HDSA::Ptr<HDSA::Vector<RealT> > uk = (*U0_samps_)[k];	  
	  bayes_model_error_objects->Apply_Sqrt_L_Mat_Inverse(uk,u_vec);

	  HDSA::Ptr<HDSA::Vector<RealT> > dpk = (*delta_prior_samps)[k];
	  dpk->set(*uk);
	  z_vec->zero();
	  bayes_model_error_objects->Apply_Gamma_Mat_Inverse(z_vec,zk);
	  RealT val = std::sqrt( 1.0 + z_vec->dot(*zk) );
	  dpk->scale(val);
	}

      std::string name = "z_prior_samples.txt";
      z_prior_samps->Write_to_File(name);
      name = "delta_prior_samples.txt";
      delta_prior_samps->Write_to_File(name);
    }

    void Compute_Posterior_Data(const HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, 
				const HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > & Sen_Op)
    {
      post_data_ = HDSA::makePtr<HDSA::Bayes_Posterior_Data<RealT> >();
      post_data_->alpha = parlist_sensitivity_->sublist("Bayes Model Error").get("alpha", 1.0);
      post_data_->Z = bayes_model_error_objects->Load_Z_Data();
      post_data_->Y = bayes_model_error_objects->Load_Y_Data();
      post_data_->N = post_data_->Z->Number_of_Vectors();
      
      post_data_->Gamma_inv_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->z);
      post_data_->Mz_inv_Gamma_inv_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->z);
      for(int k = 0; k < post_data_->N; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zk = (*post_data_->Z)[k];
	  HDSA::Ptr<HDSA::Vector<RealT> > gzk = (*post_data_->Gamma_inv_Z)[k];
	  bayes_model_error_objects->Apply_Gamma_Mat_Inverse(gzk,zk);
	  HDSA::Ptr<HDSA::Vector<RealT> > Mgzk = (*post_data_->Mz_inv_Gamma_inv_Z)[k];
	  weight_matrices->Apply_z_Weight_Mat_Inverse(Mgzk,gzk);
	}

      post_data_->G = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(post_data_->N,post_data_->N);
      for(int i = 0; i < post_data_->N; i++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zi = (*post_data_->Z)[i];
	  HDSA::Ptr<HDSA::Vector<RealT> > gzi = (*post_data_->Gamma_inv_Z)[i];
	  RealT vali = 1.0 + bayes_model_error_objects->z_star_gamma_inv_z_star_ - zi->dot(*bayes_model_error_objects->gamma_inv_z_star_);
	  for(int j = 0; j < i+1; j++)
	    {
	      HDSA::Ptr<HDSA::Vector<RealT> > zj = (*post_data_->Z)[j];
	      RealT val = vali;
	      val -= zj->dot(*bayes_model_error_objects->gamma_inv_z_star_);
	      val += zj->dot(*gzi);
	      post_data_->G->Replace_Element(i,j,val);
	    }
	}
      for(int i = 0; i < post_data_->N; i++)
	{
	  for(int j = i+1; j < post_data_->N; j++)
	    {
	      post_data_->G->Replace_Element(i,j,(*post_data_->G)(j,i));
	    }
	}

      post_data_->g_vecs = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(post_data_->N,post_data_->N);
      post_data_->Lambda = HDSA::makePtr<Std_Vector<RealT> >(post_data_->N);
      HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(post_data_->G, post_data_->g_vecs, post_data_->Lambda);

      post_data_->u_ell = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->u);      
      post_data_->u_i_ell.resize(post_data_->N);
      for(int i = 0; i < post_data_->N; i++)
	{
	  post_data_->u_i_ell[i] = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->u);
	  for(int ell = 0; ell < post_data_->N; ell++)
	    {
	      HDSA::Ptr<HDSA::Vector<RealT> > yl = (*post_data_->Y)[ell];
	      HDSA::Ptr<HDSA::Vector<RealT> > ul = (*post_data_->u_ell)[ell];
	      bayes_model_error_objects->Apply_L_Mat_Inverse(ul,yl);
	      HDSA::Ptr<HDSA::Vector<RealT> > uil = (*post_data_->u_i_ell[i])[ell];
	      bayes_model_error_objects->Apply_L_Plus_Shift_Mat_Inverse(uil,ul,(*post_data_->Lambda)(i)/post_data_->alpha);
	      uil->scale(1.0/post_data_->alpha);
	    }
	}

      post_data_->a_ell = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(post_data_->N,1);
      post_data_->b_i_ell = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(post_data_->N,post_data_->N);
      for(int ell = 0; ell < post_data_->N; ell++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zl = (*post_data_->Z)[ell];
	  RealT val_a = 1.0 - zl->dot(*bayes_model_error_objects->gamma_inv_z_star_) + bayes_model_error_objects->z_star_gamma_inv_z_star_;
	  post_data_->a_ell->Replace_Element(ell,0,val_a);
	  for(int i = 0; i < post_data_->N; i++)
	    {
	      RealT val_b = 0.0;
	      for(int k = 0; k < post_data_->N; k++)
		{
		  HDSA::Ptr<HDSA::Vector<RealT> > gzk = (*post_data_->Gamma_inv_Z)[k];
		  val_b += (*post_data_->g_vecs)(k,i)*(zl->dot(*gzk) - gzk->dot(*bayes_model_error_objects->OP_Objects_->z) + (*post_data_->a_ell)(ell,0));
		}
	      post_data_->b_i_ell->Replace_Element(i,ell,val_b);
	    }
	}
      
      post_data_->u_hat.resize(num_post_samps_);
      for(int k = 0; k < num_post_samps_; k++)
	{
	  post_data_->u_hat[k] = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->u);
	  for(int i = 0; i < post_data_->N; i++)
	    {
	      HDSA::Ptr<HDSA::Vector<RealT> > u_vec = Sen_Op->Generate_Random_u_Vector();
	      HDSA::Ptr<HDSA::Vector<RealT> > ui = (*post_data_->u_hat[k])[i];
	      bayes_model_error_objects->Apply_Sqrt_L_Plus_Shift_Mat_Inverse(ui,u_vec,(*post_data_->Lambda)(i)/post_data_->alpha);
	      ui->scale(1.0/std::sqrt(post_data_->alpha));
	    }
	}

    }

    void Posterior_Sampling(HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom,
			    HDSA::Ptr<HDSA::Vector<RealT> > & grad_nominal) const
    {
      for(int k = 0; k < post_data_->N; k++)
	{
	  HDSA::Ptr<HDSA::MultiVector<RealT> > delta_samples = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_post_samps_,bayes_model_error_objects->OP_Objects_->u);
	  HDSA::Ptr<HDSA::Vector<RealT> > delta_mean = bayes_model_error_objects->OP_Objects_->u->Clone();
	  Posterior_Discrepancy(delta_samples,delta_mean,(*post_data_->Z)[k],bayes_model_error_objects);
	  std::string name = "delta_mean_at_z" + std::to_string(k+1) + ".txt";
	  delta_mean->Write_to_File(name);
	  name = "delta_samples_at_z" + std::to_string(k+1) + ".txt";
	  delta_samples->Write_to_File(name);
	}
      
      HDSA::Ptr<HDSA::MultiVector<RealT> > z_samples = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_post_samps_,bayes_model_error_objects->OP_Objects_->z);
      HDSA::Ptr<HDSA::Vector<RealT> > z_mean = bayes_model_error_objects->OP_Objects_->z->Clone();
      Posterior_Opt_z(z_samples,z_mean,bayes_model_error_objects, Nom, grad_nominal);
      std::string name = "opt_z_mean.txt";
      z_mean->Write_to_File(name);
      name = "opt_z_samples.txt";
      z_samples->Write_to_File(name);
    }

    void Posterior_Discrepancy(const HDSA::Ptr<HDSA::MultiVector<RealT> > & delta_samples, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_mean, 
			       const HDSA::Ptr<HDSA::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects) const
    {
      HDSA::Ptr<HDSA::Vector<RealT> > dz = z->Clone();
      dz->set(*z);
      dz->axpy(-1.0,*bayes_model_error_objects->OP_Objects_->z);
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z->Clone();

      RealT val = 1.0 - dz->dot(*bayes_model_error_objects->gamma_inv_z_star_);
      std::vector<RealT> coeff_ell = post_data_->Gamma_inv_Z->dot(*dz);
      for(int ell = 0; ell < post_data_->N; ell++)
	{
	  coeff_ell[ell] += val;
	}
      std::vector<RealT> coeff_i = std::vector<RealT>(post_data_->N,0.0);
      for(int i = 0; i < post_data_->N; i++)
	{
	  RealT tmp = 0.0;
	  for(int k = 0; k < post_data_->N; k++)
	    {
	      tmp += coeff_ell[k]*(*post_data_->g_vecs)(k,i);
	    }
	  coeff_i[i] = tmp;
	}

      delta_mean->zero();
      for(int ell = 0; ell < post_data_->N; ell++)
	{
	  delta_mean->axpy(coeff_ell[ell],*(*post_data_->u_ell)[ell]);
	  for(int i = 0; i < post_data_->N; i++)
	    {
	      RealT c = -1.0*coeff_i[i]*(*post_data_->b_i_ell)(i,ell);
	      delta_mean->axpy(c,*(*post_data_->u_i_ell[i])[ell]);
	    }
	}
      delta_mean->scale(1.0/post_data_->alpha);

      for(int k = 0; k < num_post_samps_; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > delta_k = (*delta_samples)[k];
	  for(int i = 0; i < post_data_->N; i++)
	    {
	      RealT c = coeff_i[i]/std::sqrt((*post_data_->Lambda)(i));
	      delta_k->axpy(c,*(*post_data_->u_hat[k])[i]);
	    }
	  delta_k->scale(std::sqrt(post_data_->alpha));
	  delta_k->plus(*delta_mean);
	}
    }

    void Posterior_Opt_z(HDSA::Ptr<HDSA::MultiVector<RealT> > & z_samples, HDSA::Ptr<HDSA::Vector<RealT> > & z_mean,
			 HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, 
			 HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom, HDSA::Ptr<HDSA::Vector<RealT> > & grad_nominal) const
      
    {
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

      HDSA::Ptr<HDSA::Vector<RealT> > B_theta_bar = bayes_model_error_objects->OP_Objects_->z->Clone();

      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = bayes_model_error_objects->OP_Objects_->u->Clone();
      for(int ell = 0; ell < post_data_->N; ell++)
	{
	  u_vec_1->plus(*(*post_data_->u_ell)[ell]);
	  for(int i = 0; i < post_data_->N; i++)
	    {
	      RealT c = -1.0*gi_sum[i]*(*post_data_->b_i_ell)(i,ell);
	      u_vec_1->axpy(c,*(*post_data_->u_i_ell[i])[ell]);
	    }
	}
    
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = bayes_model_error_objects->OP_Objects_->u->Clone();
      bayes_model_error_objects->OP_Objects_->fs_obj->hessVec_u_u(*u_vec_2, *u_vec_1, *bayes_model_error_objects->OP_Objects_->u,*bayes_model_error_objects->OP_Objects_->z,
								  *bayes_model_error_objects->OP_Objects_->theta,false,bayes_model_error_objects->g_);
      bayes_model_error_objects->Apply_Solution_Operator_z_Jacobian_Transpose(*B_theta_bar, *u_vec_2);

      RealT c_run = 0.0;
      for(int ell = 0; ell < post_data_->N; ell++)
	{
	  RealT c = (*post_data_->u_ell)[ell]->dot(*bayes_model_error_objects->g_);
	  B_theta_bar->axpy(c,*(*post_data_->Gamma_inv_Z)[ell]);
	  c_run += c;
	}
      B_theta_bar->axpy(-1.0*c_run,*bayes_model_error_objects->gamma_inv_z_star_);

      HDSA::Ptr<HDSA::MultiVector<RealT> > Gamma_inv_w = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->z);
      for(int i = 0; i < post_data_->N; i++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > gw = (*Gamma_inv_w)[i];
	  for(int k = 0; k < post_data_->N; k++)
	    {
	      gw->axpy((*post_data_->g_vecs)(k,i),*(*post_data_->Gamma_inv_Z)[k]);
	      gw->axpy(-1.0*(*post_data_->g_vecs)(k,i),*bayes_model_error_objects->gamma_inv_z_star_);
	    }
	  for(int ell = 0; ell < post_data_->N; ell++)
	    {
	      RealT c = -1.0*( (*post_data_->b_i_ell)(i,ell) )*( (*post_data_->u_i_ell[i])[ell]->dot(*bayes_model_error_objects->g_) );
	      B_theta_bar->axpy(c,*gw);
	    }
	}

      B_theta_bar->scale(1.0/post_data_->alpha);

      Apply_Inverse_Hessian(z_mean,B_theta_bar,bayes_model_error_objects->OP_Objects_,Nom,grad_nominal);
      z_mean->scale(-1.0);
      z_mean->plus(*bayes_model_error_objects->OP_Objects_->z);
    }

    // Invert in reduced space
    void Apply_Inverse_Hessian(HDSA::Ptr<HDSA::Vector<RealT> > & x_star, const HDSA::Ptr<HDSA::Vector<RealT> > & b, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects,
			       const HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, const HDSA::Ptr<HDSA::Vector<RealT> > & grad_nominal) const
    {
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > hessian_op = HDSA::makePtr<Hessian_Operator<RealT> >(OP_Objects,grad_nominal);
      RealT tol = Nom->Get_parlist_sensitivity()->sublist("KKT Solve").get("Tolerance",1.e-5);
      std::string solver = Nom->Get_parlist_sensitivity()->sublist("KKT Solve").get("Solver","CG");
      bool verbose = Nom->Get_parlist_sensitivity()->sublist("KKT Solve").get("Verbosity",false);
      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(x_star,b,hessian_op,tol,solver,verbose);
    }

    // Overload HDSA::Linear_Operator to take matrix vector products
    template <class ScalarType>
    class Hessian_Operator : public HDSA::Linear_Operator<ScalarType>
    {
      HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > OP_Objects_;
      HDSA::Ptr<HDSA::Vector<ScalarType> > grad_;
      
    public:
      
      Hessian_Operator(const HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > & OP_Objects, const HDSA::Ptr<HDSA::Vector<ScalarType> > & grad_nominal): OP_Objects_(OP_Objects), grad_(grad_nominal)
      { }
      
      //! Dtor
      ~Hessian_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	OP_Objects_->rs_obj->hessVec_z_z(*y,*x,*OP_Objects_->z,*OP_Objects_->theta,false,grad_);  
      }
      
    };
    
  };
  
}

#endif
