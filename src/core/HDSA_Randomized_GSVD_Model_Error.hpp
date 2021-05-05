#ifndef HDSA_RANDOMIZED_GSVD_MODEL_ERROR_HPP
#define HDSA_RANDOMIZED_GSVD_MODEL_ERROR_HPP

// This class executes the randomized GSVD solver for model error

namespace HDSA
{

  template <class RealT>
  class Randomized_GSVD_Model_Error
  {
  private:
    bool reduced_space_sen_;
    int num_sing_vals_, p_, kpp_, q_, nonzero_z_dim_, state_dim_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    HDSA::Ptr<HDSA::Vector<RealT> > S_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_;
    HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > V_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Mz_V_;

  public:
  
    Randomized_GSVD_Model_Error(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
				const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
				const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    {
      reduced_space_sen_ = parlist_sensitivity->sublist("Formulation").get("Reduced Space Sensitivities",true);
      num_sing_vals_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Singular Values",1); 
      p_ = parlist_sensitivity->sublist("GSVD Solver").get("Oversampling Factor",20);
      q_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Subspace Iterations",0);
      kpp_ = num_sing_vals_ + p_;
    }
    
    void Compute(void)
    {
      // This function executes the randomized GSVD
      std::clock_t timer_gsvd = std::clock();

      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,kpp_);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      HDSA::Opt_Problem_Objects_Model_Error<RealT> &OP_Objects_subcomm_model_error = static_cast<HDSA::Opt_Problem_Objects_Model_Error<RealT>&>(*OP_Objects_subcomm);
      HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_objects = OP_Objects_subcomm_model_error.Get_Model_Error_Objects();
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects_subcomm->z->Clone();
      OP_Objects_subcomm->rs_obj->gradient_z(*grad_nominal,*OP_Objects_subcomm->z,*model_error_objects->OP_Objects_->theta,true);
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      nonzero_z_dim_ = Nom_subcomm->Get_nonzero_z_dim();
      state_dim_ = OP_Objects_subcomm->u->dimension();
      HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm;
      if(reduced_space_sen_)
      	{
      	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
      	}
      else
      	{
	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
      	}
      // Need to perform sketch of the range of D  
      if(comm_->getRank() == 0)
      	{
      	  std::cout << " " << std::endl;
      	  std::cout << "Beginning range sketching for local sensitivity number " << sample_index_ << std::endl;
      	  std::cout << " " << std::endl;
      	}
  
      // Instantiate matrix to interface with LA routines   
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);   
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WY = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 

      comm_->barrier();
      std::clock_t timer_Loop_1 = std::clock();
      for(int k = 0; k < kpp_; k++)
      	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
      	    {
	      std::clock_t timer_range_iter_k = std::clock();
	      HDSA::Ptr<HDSA::Vector<RealT> >  u_vec_random = Sen_Op_subcomm->Generate_Random_u_Vector();
	      u_vec_random->scale(std::sqrt(1.0+std::pow(model_error_objects->Mz_star_->norm(),2.0))); 
	      HDSA::Ptr<HDSA::Vector<RealT> >  z_vec_random = Sen_Op_subcomm->Generate_Random_z_Vector();
	      z_vec_random->scale(model_error_objects->g_->norm());

	      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_1 = OP_Objects_subcomm->u->Clone();
      	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();

	      // Apply B to random vector
	      OP_Objects_subcomm->fs_obj->hessVec_u_u(*u_vec_subcomm_1, *u_vec_random, *OP_Objects_subcomm->u, 
						      *OP_Objects_subcomm->z,*model_error_objects->OP_Objects_->theta,false, model_error_objects->g_);
	      model_error_objects->Apply_Solution_Operator_z_Jacobian_Transpose(*z_vec_subcomm_1, *u_vec_subcomm_1);
	      weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_2,z_vec_random);
	      z_vec_subcomm_1->plus(*z_vec_subcomm_2);

	      // Apply H^{-1} to B*v
	      z_vec_subcomm_2->zero();
	      Apply_Inverse_Hessian(z_vec_subcomm_2, z_vec_subcomm_1,OP_Objects_subcomm,Nom_subcomm,grad_nominal);
	      Y->Write_Vector_to_Column(k,z_vec_subcomm_2);

	      // Apply weight matrix
	      z_vec_subcomm_1->zero();
	      weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_1,z_vec_subcomm_2);
	      WY->Write_Vector_to_Column(k,z_vec_subcomm_1);
	      
	      RealT Time_range_iter_k = static_cast<RealT>(std::clock()-timer_range_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
      	      if(subcomm->getRank() == 0)
      		{
      		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the range sketching loop for local sensitivity number " << sample_index_ << " in " << Time_range_iter_k << " seconds" << std::endl;
      		}
      	    }
      	}
  
      comm_->barrier();
      proc_dist->Broadcast_Matrix(Y);
      proc_dist->Broadcast_Matrix(WY);
      RealT Time_Loop_1 = static_cast<RealT>(std::clock()-timer_Loop_1)/static_cast<RealT>(CLOCKS_PER_SEC);

      // Need to orthogonalize the sketch with respect to the z weight matrix
      std::clock_t timer_ortho_z = std::clock();

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > SQ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_); 
      HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Y, WY, Q, R, SQ);

      comm_->barrier();
      RealT Time_ortho_z = static_cast<RealT>(std::clock()-timer_ortho_z)/static_cast<RealT>(CLOCKS_PER_SEC);
      if(comm_->getRank() == 0)
	{
	  std::cout << "Completed orthogonalization for local sensitivity number " << sample_index_<< " in " << Time_ortho_z << " seconds" << std::endl;
	}
      
      // Subspace iteration
      std::clock_t timer_subspace_iter = std::clock();
      
      if((comm_->getRank() == 0) & (q_ > 0) )
      	{
      	  std::cout << " " << std::endl;
      	  std::cout << "Beginning subspace iteration for local sensitivity number " << sample_index_ << std::endl;
      	  std::cout << " " << std::endl;
      	}
      
      for(int j = 0; j < q_; j++)
	{
	  HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > Y_subspace_iter_param = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);
	  HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > WY_subspace_iter_param = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);

	  Y_subspace_iter_param->Set_a(1.0);
	  Y_subspace_iter_param->Set_u(*model_error_objects->g_);
	  Y_subspace_iter_param->Set_z(*model_error_objects->Mz_star_);
	  WY_subspace_iter_param->Set_a(1.0-model_error_objects->Einv_Mz_star_->dot(*model_error_objects->Mz_star_));
	  WY_subspace_iter_param->Set_u(*model_error_objects->coeff_Linv_g_);
	  WY_subspace_iter_param->Set_z(*model_error_objects->N_min_Einv_Mz_star_);

      	  for(int k = 0; k < kpp_; k++)
      	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
      		{
		  std::clock_t timer_subspace_iter_k = std::clock();
							
		  // F_uu*S_z*H^{-1}
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_1 = OP_Objects_subcomm->u->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_2 = OP_Objects_subcomm->u->Clone();
		  SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);
		  Apply_Inverse_Hessian(z_vec_subcomm_2, z_vec_subcomm_1,OP_Objects_subcomm,Nom_subcomm,grad_nominal);
		  z_vec_subcomm_2->scale(-1.0);

		  model_error_objects->Apply_Solution_Operator_z_Jacobian(*u_vec_subcomm_1, *z_vec_subcomm_2);
		  OP_Objects_subcomm->fs_obj->hessVec_u_u(*u_vec_subcomm_2, *u_vec_subcomm_1, *OP_Objects_subcomm->u, 
							  *OP_Objects_subcomm->z,*model_error_objects->OP_Objects_->theta,false, model_error_objects->g_);

		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_1,z_vec_subcomm_2);
		  Y_subspace_iter_param->Set_uk(k,*u_vec_subcomm_2);
		  Y_subspace_iter_param->Set_zk(k,*z_vec_subcomm_1);
		  Y_subspace_iter_param->Set_bk(k,0.0);

		  // Compute WY_subspace_iter
		  model_error_objects->Apply_L_Mat_Inverse(u_vec_subcomm_1,u_vec_subcomm_2);
		  u_vec_subcomm_1->scale(model_error_objects->coeff_);
		  model_error_objects->Apply_N(z_vec_subcomm_2,z_vec_subcomm_1);
		  WY_subspace_iter_param->Set_uk(k,*u_vec_subcomm_1);
		  WY_subspace_iter_param->Set_zk(k,*z_vec_subcomm_2);
		  WY_subspace_iter_param->Set_bk(k,-1.0*model_error_objects->Einv_Mz_star_->dot(*z_vec_subcomm_1));

		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the first loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
		    }
		}
      	    }
	  
	  std::clock_t timer_ortho = std::clock();
	  comm_->barrier();
	  Y_subspace_iter_param->Broadcast_Data();
	  WY_subspace_iter_param->Broadcast_Data();
	  comm_->barrier();
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yuk_WYu = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yu_WYuk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yzk_WYz = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yz_WYzk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yu_WYu = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(1,1);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Yz_WYz = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(1,1);

	  Y_subspace_iter_param->uk->Multiply(U,WY_subspace_iter_param->uk,true,false);
	  Y_subspace_iter_param->zk->Multiply(U,WY_subspace_iter_param->zk,true,false);
	  Y_subspace_iter_param->uk->Multiply(Yuk_WYu,WY_subspace_iter_param->u,true,false);
	  WY_subspace_iter_param->uk->Multiply(Yu_WYuk,Y_subspace_iter_param->u,true,false);
	  Y_subspace_iter_param->zk->Multiply(Yzk_WYz,WY_subspace_iter_param->z,true,false);
	  WY_subspace_iter_param->zk->Multiply(Yz_WYzk,Y_subspace_iter_param->z,true,false);
	  Y_subspace_iter_param->u->Multiply(Yu_WYu,WY_subspace_iter_param->u,true,false);
	  Y_subspace_iter_param->z->Multiply(Yz_WYz,WY_subspace_iter_param->z,true,false);

	  for(int i = 0; i < kpp_; i++)
	    {
	      for(int k = 0; k < kpp_; k++)
		{
		  RealT val = 0.0;
		  val += (*U)(i,k)*( (Y_subspace_iter_param->a)*(WY_subspace_iter_param->a) + (*Yz_WYz)(0,0) );
		  val += (*Yuk_WYu)(i,0)*( (Y_subspace_iter_param->a)*(*WY_subspace_iter_param->bk)(0,k) + (*Yz_WYzk)(k,0) );
		  val += (*Yu_WYuk)(k,0)*( (WY_subspace_iter_param->a)*(*Y_subspace_iter_param->bk)(0,i) + (*Yzk_WYz)(i,0) );
		  val += (*Yu_WYu)(0,0)*( (*Y_subspace_iter_param->bk)(0,i)*(*WY_subspace_iter_param->bk)(0,k) + (*Z)(i,k) );
		  C->Replace_Element(i,k,val);
		}
	    }
	  
	  // Compute R=chol(C)
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R1 = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
	  HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(C,R1);
      
	  // Need to solve the linear systems R*x_i = e_i for i=1,2,...,n. Here e_i is the ith standard basis vector in euclidean space, hence its ith entry is 1 and others are 0
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R1inv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
	  for(int i = 0; i < kpp_; i++)
	    {
	      // Solve the linear system R*x = e_i
	      HDSA::Ptr<HDSA::Vector<RealT> > ei = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
	      ei->basis(i);
	      HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
	      HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, ei, R1);
	      R1inv->Write_Vector_to_Column(i,x);
	    }

	  HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > Q_subspace_iter_param = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);
	  HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > WQ_subspace_iter_param = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);

	  Y_subspace_iter_param->uk->Multiply(Q_subspace_iter_param->uk,R1inv);
	  Y_subspace_iter_param->zk->Multiply(Q_subspace_iter_param->zk,R1inv);
	  Y_subspace_iter_param->bk->Multiply(Q_subspace_iter_param->bk,R1inv);
	  Q_subspace_iter_param->a = Y_subspace_iter_param->a;
	  Q_subspace_iter_param->u = Y_subspace_iter_param->u;
	  Q_subspace_iter_param->z = Y_subspace_iter_param->z;

	  WY_subspace_iter_param->uk->Multiply(WQ_subspace_iter_param->uk,R1inv);
	  WY_subspace_iter_param->zk->Multiply(WQ_subspace_iter_param->zk,R1inv);
	  WY_subspace_iter_param->bk->Multiply(WQ_subspace_iter_param->bk,R1inv);
	  WQ_subspace_iter_param->a = WY_subspace_iter_param->a;
	  WQ_subspace_iter_param->u = WY_subspace_iter_param->u;
	  WQ_subspace_iter_param->z = WY_subspace_iter_param->z;

	  comm_->barrier();
	  RealT Time_ortho = static_cast<RealT>(std::clock()-timer_ortho)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << "Completed parameter orthogonalization for local sensitivity number " << sample_index_<< " in " << Time_ortho << " seconds" << std::endl;
	    }
      
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_subspace_iter_opt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WY_subspace_iter_opt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
      
	  HDSA::Ptr<HDSA::Vector<RealT> > u_vec = OP_Objects_subcomm->u->Clone();
	  HDSA::Ptr<HDSA::Vector<RealT> > z_vec = OP_Objects_subcomm->z->Clone();
	  HDSA::Ptr<HDSA::Vector<RealT> > Mz = OP_Objects_subcomm->z->Clone();
	  WQ_subspace_iter_param->u->Write_Column_to_Vector(0,u_vec);
	  WQ_subspace_iter_param->z->Write_Column_to_Vector(0,z_vec);
	  RealT c1 = WQ_subspace_iter_param->a + z_vec->dot(*model_error_objects->Mz_star_);
	  RealT c4 = u_vec->dot(*model_error_objects->g_);
	  weight_matrices_subcomm->Apply_z_Weight_Mat(Mz,z_vec);

	  for(int k = 0; k < kpp_; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
		{
		  std::clock_t timer_subspace_iter_k = std::clock();
		  
		  HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_1 = OP_Objects_subcomm->u->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_2 = OP_Objects_subcomm->u->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();
		  WQ_subspace_iter_param->uk->Write_Column_to_Vector(k,u_vec_subcomm_1);
		  WQ_subspace_iter_param->zk->Write_Column_to_Vector(k,z_vec_subcomm_1);

		  RealT c2 = model_error_objects->g_->dot(*u_vec_subcomm_1);
		  RealT c3 = (*WQ_subspace_iter_param->bk)(0,k) + model_error_objects->Mz_star_->dot(*z_vec_subcomm_1);
		  
		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_2,z_vec_subcomm_1);

		  u_vec_subcomm_1->scale(c1);
		  u_vec_subcomm_1->axpy(c3,*u_vec);

		  OP_Objects_subcomm->fs_obj->hessVec_u_u(*u_vec_subcomm_2, *u_vec_subcomm_1, *OP_Objects_subcomm->u, 
							  *OP_Objects_subcomm->z,*model_error_objects->OP_Objects_->theta,false, model_error_objects->g_);
		  model_error_objects->Apply_Solution_Operator_z_Jacobian_Transpose(*z_vec_subcomm_1, *u_vec_subcomm_2);
        
		  z_vec_subcomm_1->axpy(c2,*Mz);
		  z_vec_subcomm_1->axpy(c4,*z_vec_subcomm_2);
		  
		  Apply_Inverse_Hessian(z_vec_subcomm_2, z_vec_subcomm_1,OP_Objects_subcomm,Nom_subcomm,grad_nominal);
		  z_vec_subcomm_2->scale(-1.0);
	      
		  Y_subspace_iter_opt->Write_Vector_to_Column(k,z_vec_subcomm_2);
		  
		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_1,z_vec_subcomm_2);
		  WY_subspace_iter_opt->Write_Vector_to_Column(k,z_vec_subcomm_1);
		  
		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the second loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
		    }
		}
	    }

	  comm_->barrier();
	  proc_dist->Broadcast_Matrix(Y_subspace_iter_opt);   
	  proc_dist->Broadcast_Matrix(WY_subspace_iter_opt);
	  comm_->barrier();
	  HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Y_subspace_iter_opt, WY_subspace_iter_opt, Q, R1, SQ);    

	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Completed subspace iteration " << j+1 << " out of " << q_  << " for local sensitivity number " << sample_index_<< std::endl;
	    }
	}
      
      comm_->barrier();
      RealT Time_Loop_SubIter = static_cast<RealT>(std::clock()-timer_subspace_iter)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      // Need to project D onto the sketched subspace
      if(comm_->getRank() == 0)
      	{
	  std::cout << " " << std::endl;
      	  std::cout << "Beginning projection onto the sketched subspace for local sensitivity number " << sample_index_ << std::endl;
      	  std::cout << " " << std::endl;
      	}

      HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > Q_proj = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);
      HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > WQ_proj = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);

      WQ_proj->Set_a(1);
      WQ_proj->Set_u(*model_error_objects->g_);
      WQ_proj->Set_z(*model_error_objects->Mz_star_);
      Q_proj->Set_a(1.0-model_error_objects->Einv_Mz_star_->dot(*model_error_objects->Mz_star_));
      Q_proj->Set_u(*model_error_objects->coeff_Linv_g_);
      Q_proj->Set_z(*model_error_objects->N_min_Einv_Mz_star_);
      
      comm_->barrier();
      std::clock_t timer_Loop_2 = std::clock();
      for(int k = 0; k < kpp_; k++)
      	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
      	    { 
	      std::clock_t timer_proj_iter_k = std::clock();
	      // F_uu*S_z*H^{-1}
	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_1 = OP_Objects_subcomm->u->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_subcomm_2 = OP_Objects_subcomm->u->Clone();
	      SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);

	      Apply_Inverse_Hessian(z_vec_subcomm_2, z_vec_subcomm_1,OP_Objects_subcomm,Nom_subcomm,grad_nominal);
	      z_vec_subcomm_2->scale(-1.0);

	      model_error_objects->Apply_Solution_Operator_z_Jacobian(*u_vec_subcomm_1, *z_vec_subcomm_2);
	      OP_Objects_subcomm->fs_obj->hessVec_u_u(*u_vec_subcomm_2, *u_vec_subcomm_1, *OP_Objects_subcomm->u, 
						      *OP_Objects_subcomm->z,*model_error_objects->OP_Objects_->theta,false, model_error_objects->g_);

	      weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_1,z_vec_subcomm_2);
	      WQ_proj->Set_uk(k,*u_vec_subcomm_2);
	      WQ_proj->Set_zk(k,*z_vec_subcomm_1);
	      WQ_proj->Set_bk(k,0.0);

	      // Compute WQ_proj
	      model_error_objects->Apply_L_Mat_Inverse(u_vec_subcomm_1,u_vec_subcomm_2);
	      u_vec_subcomm_1->scale(model_error_objects->coeff_);
	      model_error_objects->Apply_N(z_vec_subcomm_2,z_vec_subcomm_1);
	      Q_proj->Set_uk(k,*u_vec_subcomm_1);
	      Q_proj->Set_zk(k,*z_vec_subcomm_2);
	      Q_proj->Set_bk(k,-1.0*model_error_objects->Einv_Mz_star_->dot(*z_vec_subcomm_1));

	      RealT Time_proj_iter_k = static_cast<RealT>(std::clock()-timer_proj_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
      	      if(subcomm->getRank() == 0)
      		{
      		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the projection loop for local sensitivity number " << sample_index_ << " in " << Time_proj_iter_k << " seconds" << std::endl;
      		}
      	    }
      	}
      
      comm_->barrier();
      Q_proj->Broadcast_Data();
      WQ_proj->Broadcast_Data();
      comm_->barrier();

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Quk_WQu = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qu_WQuk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qzk_WQz = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qz_WQzk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,1);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qu_WQu = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(1,1);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Qz_WQz = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(1,1);

      Q_proj->uk->Multiply(U_proj,WQ_proj->uk,true,false);
      Q_proj->zk->Multiply(Z_proj,WQ_proj->zk,true,false);
      Q_proj->uk->Multiply(Quk_WQu,WQ_proj->u,true,false);
      WQ_proj->uk->Multiply(Qu_WQuk,Q_proj->u,true,false);
      Q_proj->zk->Multiply(Qzk_WQz,WQ_proj->z,true,false);
      WQ_proj->zk->Multiply(Qz_WQzk,Q_proj->z,true,false);
      Q_proj->u->Multiply(Qu_WQu,WQ_proj->u,true,false);
      Q_proj->z->Multiply(Qz_WQz,WQ_proj->z,true,false);

      for(int i = 0; i < kpp_; i++)
	{
	  for(int k = 0; k < kpp_; k++)
	    {
	      RealT val = 0.0;
	      val += (*U_proj)(i,k)*( (Q_proj->a)*(WQ_proj->a) + (*Qz_WQz)(0,0) );
	      val += (*Quk_WQu)(i,0)*( (Q_proj->a)*(*WQ_proj->bk)(0,k) + (*Qz_WQzk)(k,0) );
	      val += (*Qu_WQuk)(k,0)*( (WQ_proj->a)*(*Q_proj->bk)(0,i) + (*Qzk_WQz)(i,0) );
	      val += (*Qu_WQu)(0,0)*( (*Q_proj->bk)(0,i)*(*WQ_proj->bk)(0,k) + (*Z_proj)(i,k) );
	      C_proj->Replace_Element(i,k,val);
	    }
	}
      
      // Compute R=chol(C)
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(C_proj,R_proj);

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_proj_inv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      // Need to solve the linear systems R*x_i = e_i for i=1,2,...,n. Here e_i is the ith standard basis vector in euclidean space, hence its ith entry is 1 and others are 0
      for(int i = 0; i < kpp_; i++)
	{
	  // Solve the linear system R*x = e_i
	  HDSA::Ptr<HDSA::Vector<RealT> > ei = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
	  ei->basis(i);
	  HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
	  HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, ei, R_proj);
	  R_proj_inv->Write_Vector_to_Column(i,x);
	}

      HDSA::Ptr<HDSA::Model_Error_Kronecker_Matrix<RealT> > Q_param = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);
      Q_proj->uk->Multiply(Q_param->uk,R_proj_inv);
      Q_proj->zk->Multiply(Q_param->zk,R_proj_inv);
      Q_proj->bk->Multiply(Q_param->bk,R_proj_inv);
      Q_param->a = Q_proj->a;
      Q_param->u = Q_proj->u;
      Q_param->z = Q_proj->z;

      RealT Time_Loop_2 = static_cast<RealT>(std::clock()-timer_Loop_2)/static_cast<RealT>(CLOCKS_PER_SEC);
  
      // Need to compute the SVD of R_proj
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      S_ = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
      HDSA::Linear_Algebra::SVD<RealT>(R_proj, U, VT, S_); 
      U_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
      Q->Multiply(U_,VT,false,true);
      V_ = HDSA::makePtr<HDSA::Model_Error_Kronecker_Matrix<RealT> >(state_dim_,nonzero_z_dim_,kpp_,proc_dist);
      Q_param->uk->Multiply(V_->uk,U);
      Q_param->zk->Multiply(V_->zk,U);
      Q_param->bk->Multiply(V_->bk,U);
      V_->a = Q_param->a;
      V_->u = Q_param->u;
      V_->z = Q_param->z;

      Mz_V_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_+1);
      HDSA::Ptr<HDSA::Vector<RealT> > z1 = OP_Objects_subcomm->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z2 = OP_Objects_subcomm->z->Clone();
      for(int k = 0; k < kpp_; k++)
	{
	  V_->zk->Write_Column_to_Vector(k,z1);
	  weight_matrices_subcomm->Apply_z_Weight_Mat(z2,z1);
	  Mz_V_->Write_Vector_to_Column(k,z2);
	}
      V_->z->Write_Column_to_Vector(0,z1);
      weight_matrices_subcomm->Apply_z_Weight_Mat(z2,z1);
      Mz_V_->Write_Vector_to_Column(kpp_,z2);

      // The columns of U contain the approximate opt singular vectors
      // The columns of V contain the approximate parameter singular vectors (decomposed in kronecker form)
        
      RealT Time_gsvd = static_cast<RealT>(std::clock()-timer_gsvd)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Total time for GSVD computation: " << Time_gsvd << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by the first loop: " << Time_Loop_1 << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by the subspace iteration: " << Time_Loop_SubIter << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by the second loop: " << Time_Loop_2 << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by everything else: " << Time_gsvd-Time_Loop_1-Time_Loop_2-Time_Loop_SubIter << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
      
      Write_Solution(OP_Objects_subcomm->z);
    }

    void Write_Solution(HDSA::Ptr<HDSA::Vector<RealT> > & z)
    {
      std::clock_t timer_write = std::clock();
      if(comm_->getRank() == 0)
	{
	  // Write to text files
	  std::string name;
	  std::ofstream fout;

	  name = "z_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < z->dimension(); i++)
	    {
	      if(z->Is_entry_zero(i))
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << 0.0 << std::setw(20);
		    }
		}
	      else
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << (*U_)(z->Get_map_full_to_reduced(i),k) << std::setw(20);
		    }
		}
	      fout << " " << std::endl;
	    }
	  fout.close();
	
	  name = "theta_Singular_Vector_zk_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < z->dimension(); i++)
	    {
	      if(z->Is_entry_zero(i))
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << 0.0 << std::setw(20);
		    }
		}
	      else
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << (*V_->zk)(z->Get_map_full_to_reduced(i),k) << std::setw(20);
		    }
		}
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_z_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < z->dimension(); i++)
	    {
	      if(z->Is_entry_zero(i))
		{
		  fout << 0.0 << std::setw(20);
		}
	      else
		{
		  fout << (*V_->z)(z->Get_map_full_to_reduced(i),0) << std::setw(20);
		}
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_Mzk_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < z->dimension(); i++)
	    {
	      if(z->Is_entry_zero(i))
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << 0.0 << std::setw(20);
		    }
		}
	      else
		{
		  for(int k = 0; k < num_sing_vals_; k++)
		    {   
		      fout << (*Mz_V_)(z->Get_map_full_to_reduced(i),k) << std::setw(20);
		    }
		}
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_Mz_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < z->dimension(); i++)
	    {
	      if(z->Is_entry_zero(i))
		{
		  fout << 0.0 << std::setw(20);
		}
	      else
		{
		  fout << (*Mz_V_)(z->Get_map_full_to_reduced(i),kpp_) << std::setw(20);
		}
	      fout << " " << std::endl;
	    }
	  fout.close();	  
	  name = "theta_Singular_Vector_uk_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < state_dim_; i++)
	    {
	      for(int k = 0; k < num_sing_vals_; k++)
		{   
		  fout << (*V_->uk)(i,k) << std::setw(20);
		}
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_u_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < state_dim_; i++)
	    {
	      fout << (*V_->u)(i,0) << std::setw(20);
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_bk_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {   
	      fout << (*V_->bk)(0,k) << std::setw(20);
	      fout << " " << std::endl;
	    }
	  fout.close();
	  name = "theta_Singular_Vector_a_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  fout << V_->a << std::endl;
	  fout.close();

	  name = "Singular_Values_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      fout << (*S_)(k) << std::setw(20);
	    }
	  fout.close();
	     
	}
      comm_->barrier();
      RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);        

      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}

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
