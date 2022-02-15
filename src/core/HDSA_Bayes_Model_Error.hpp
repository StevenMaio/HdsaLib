#ifndef HDSA_BAYES_MODEL_ERROR_HPP
#define HDSA_BAYES_MODEL_ERROR_HPP

// This class executes the randomized GSVD solver for model error

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
    }
    
    void Compute(void)
    {
     
      int num_procs = 1;
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,num_procs);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      HDSA::Opt_Problem_Objects_Bayes_Model_Error<RealT> &OP_Objects_subcomm_bayes_model_error = static_cast<HDSA::Opt_Problem_Objects_Bayes_Model_Error<RealT>&>(*OP_Objects_subcomm);
      HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > bayes_model_error_objects = OP_Objects_subcomm_bayes_model_error.Get_Model_Error_Objects();
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects_subcomm->z->Clone();
      OP_Objects_subcomm->rs_obj->gradient_z(*grad_nominal,*OP_Objects_subcomm->z,*bayes_model_error_objects->OP_Objects_->theta,true);
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());

      Compute_Prior_Samples(bayes_model_error_objects,weight_matrices_subcomm,Sen_Op_subcomm);
      Compute_Posterior_Data(bayes_model_error_objects,weight_matrices_subcomm,Sen_Op_subcomm);   

    }

    void Compute_Prior_Samples(HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, 
			       HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > & Sen_Op)
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

    void Compute_Posterior_Data(HDSA::Ptr<HDSA::Bayes_Model_Error_Objects<RealT> > & bayes_model_error_objects, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices, 
				HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > & Sen_Op)
    {
      post_data_ = HDSA::makePtr<HDSA::Bayes_Posterior_Data<RealT> >();
      post_data_->alpha = parlist_sensitivity_->sublist("Bayes Model Error").get("alpha", 1.0);
      post_data_->Z = bayes_model_error_objects->Load_Z_Data();
      post_data_->Y = bayes_model_error_objects->Load_Y_Data();
      post_data_->N = post_data_->Z->Number_of_Vectors();
      
      post_data_->Gamma_inv_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(post_data_->N,bayes_model_error_objects->OP_Objects_->z);
      for(int k = 0; k < post_data_->N; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > zk = (*post_data_->Z)[k];
	  HDSA::Ptr<HDSA::Vector<RealT> > gzk = (*post_data_->Gamma_inv_Z)[k];
	  bayes_model_error_objects->Apply_Gamma_Mat_Inverse(gzk,zk);
	}

      // Need to conintue implementation

    }


    // Invert in reduced space
    void Apply_Inverse_Hessian(HDSA::Ptr<HDSA::Vector<RealT> > & x_star, const HDSA::Ptr<HDSA::Vector<RealT> > & b, HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects,
			       HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, HDSA::Ptr<HDSA::Vector<RealT> > & grad_nominal)
    {
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > hessian_op = HDSA::makePtr<Hessian_Operator<RealT> >(OP_Objects,Nom,grad_nominal);
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
      HDSA::Ptr<HDSA::Nominal_Data<ScalarType> > Nom_;
      HDSA::Ptr<HDSA::Vector<ScalarType> > grad_;
      ScalarType coeff_;
      
    public:
      
      Hessian_Operator(const HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > & OP_Objects, const HDSA::Ptr<HDSA::Nominal_Data<ScalarType> > & Nom, 
		       const HDSA::Ptr<HDSA::Vector<ScalarType> > & grad_nominal): OP_Objects_(OP_Objects), Nom_(Nom), grad_(grad_nominal)
      { 
	Nom_->Set_grad_nominal(grad_);
      }
      
      //! Dtor
      ~Hessian_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	OP_Objects_->rs_obj->hessVec_z_z(*y,*x,*OP_Objects_->z,*OP_Objects_->theta,false, grad_);
	Nom_->Apply_Regularization_Update(y,x);   
      }
      
    };
    
  };
  
}

#endif
