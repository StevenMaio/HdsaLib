#ifndef HDSA_RANDOMIZED_HESSIAN_EVP_HPP
#define HDSA_RANDOMIZED_HESSIAN_EVP_HPP

namespace HDSA
{

  template <class RealT>
  class Randomized_Hessian_EVP{
  private:
    std::string eig_alg_;
    int oversampling_factor_;
    int max_iter_;
    RealT eig_ratio_tol_;
    RealT residual_bound_;
    RealT residual_bound_tol_;
    int initial_samples_;
    int nonzero_z_dim_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    int vecs_so_far_;
    HDSA::Ptr<HDSA::Vector<RealT> > S_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_;
 
  public:
  
    Randomized_Hessian_EVP(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
			   const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
			   const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), 
      weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    {
      eig_alg_ = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Eigenvalue Algorithm", "Nystrom");
      int target_rank = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Initial Target Rank", 1);
      oversampling_factor_ = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Oversampling Factor", 10);
      max_iter_ = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Maximum Number of Iterations", 2);
      eig_ratio_tol_ = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Eigenvalue Ratio Tolerance", 1.e-3);
      residual_bound_ = 1.e4;
      residual_bound_tol_ = parlist_sensitivity->sublist("Randomized Hessian EVP").get("Residual Bound Tolerance", 1.e-1);
      initial_samples_ = target_rank + oversampling_factor_;
    }
    
    void Compute(void)
    {
      // This function executes the algorithm
      std::clock_t timer_alg = std::clock();
      
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,initial_samples_);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
  
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      nonzero_z_dim_ = Nom_subcomm->Get_nonzero_z_dim();

      // There are "numSubComm" subcommunicators. They will compute "initial_samples_" matvecs on the first pass, and then augment the data with "numSubComm" matvecs for each subsequent iteration
      // The maximum number of matrix vector products will be = initial_samples_ + (max_iter_-1)*numSubComm
      int max_matvecs = initial_samples_ + (max_iter_-1)*proc_dist->Get_numSubComm();
      // Preallocate arrays
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,max_matvecs);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > AQ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,max_matvecs);  
      int iteration_count = 1;
      vecs_so_far_ = 0;
      RealT current_eig_ratio = eig_ratio_tol_ + 1.0;
      
      // Precompute nominal gradient
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects_subcomm->z->Clone();
      OP_Objects_subcomm->rs_obj->gradient_z(*grad_nominal,*OP_Objects_subcomm->z, *OP_Objects_subcomm->theta);
      Nom_subcomm->Set_grad_nominal(grad_nominal);
      
      while( ( (current_eig_ratio > eig_ratio_tol_) || (residual_bound_ > residual_bound_tol_) ) && (iteration_count <= max_iter_) )
	{
	  int num_samples = proc_dist->Get_numSubComm();
	  if(iteration_count==1)
	    {
	      num_samples = initial_samples_;
	    }
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,num_samples);
	  
	  // Need to perform sketch of the range of K 
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Beginning range finder sketch on iteration number " << iteration_count << " and sample number " << sample_index_ << std::endl;
	      std::cout << " " << std::endl;
	    }      
	  
	  comm_->barrier();
	  std::default_random_engine generator;
	  bool time_seed = parlist_sensitivity_->sublist("Formulation").get("System Time Seed",false);
	  int seed;
	  if(time_seed)
	    {
	      seed = time(NULL)*(1+proc_dist->Get_Comm_Split_Ranks()[0]);
	    }
	  else
	    {
	      seed = 477*(1+proc_dist->Get_Comm_Split_Ranks()[0]);
	    }
	  generator.seed(seed);
	  std::normal_distribution<RealT> distribution = std::normal_distribution<RealT>(0.0,1.0);
	  
	  comm_->barrier();
	  for(int k = 0; k < num_samples; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
		{
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_random = OP_Objects_subcomm->z->Clone();
		  // Populate vectors with standard normal samples
		  for(int l = 0; l < nonzero_z_dim_; l++)
		    {
		      z_vec_random->Replace_Element(z_vec_random->Get_map_reduced_to_full(l),distribution(generator));
		    }
		  
		  OP_Objects_subcomm->rs_obj->hessVec_z_z(*z_vec_subcomm,*z_vec_random,*OP_Objects_subcomm->z,*OP_Objects_subcomm->theta,false,grad_nominal);
		  Nom_subcomm->Apply_Regularization_Update(z_vec_subcomm,z_vec_random);
		  Y->Write_Vector_to_Column(k,z_vec_subcomm);
		  std::cout << "Completed matrix vector product " << k+1 << " out of " << num_samples << std::endl;
		}
	    }
	  
	  comm_->barrier();   
	  for(int k = 0; k < num_samples; k++)
	    {
	      char *buff = (char*)Y->Get_Element_Ptr(0,k);
	      comm_->broadcast(proc_dist->Get_Procs_Loop_Distribution()[k][0],8*nonzero_z_dim_,buff);
	    }
	  comm_->barrier();
	  
	  // Need to orthogonalize columns of Y against Q(:,1:vecs_so_far_) and write into Q(:,(vecs_so_far_+1):(vecs_so_far_+num_samples))
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > QI = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_+num_samples);
	  for(int k = 0; k < vecs_so_far_; k++)
	    {
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  QI->Replace_Element(i,k,(*Q)(i,k));
		}
	    }
	  
	  for(int k = 0; k < num_samples; k++)
	    {
	      // Normalized Y for numerical stability
	      RealT normalize = 0.0;
	      for(int j = 0; j < nonzero_z_dim_; j++)
		{
		  normalize += (*Y)(j,k)*(*Y)(j,k);
		}
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  QI->Replace_Element(i,vecs_so_far_+k,(*Y)(i,k)/std::sqrt(normalize));
		}
	    }
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_+num_samples);
	  HDSA::Linear_Algebra::QR_Factorization<RealT>(QI, Z);
	  for(int k = 0; k < num_samples; k++)
	    {
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  Q->Replace_Element(i,vecs_so_far_+k,(*Z)(i,vecs_so_far_+k));
		}
	    }
	  
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Beginning eigenvalue estimation on iteration number " << iteration_count << " and sample number " << sample_index_ << std::endl;
	      std::cout << " " << std::endl;
	    } 
	  
	  // Compute A*Q for new columns of Q
	  comm_->barrier();
	  for(int k = 0; k < num_samples; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
		{
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_in = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_out = OP_Objects_subcomm->z->Clone();
				  
		  Q->Write_Column_to_Vector(vecs_so_far_+k,z_vec_subcomm_in);
		  OP_Objects_subcomm->rs_obj->hessVec_z_z(*z_vec_subcomm_out,*z_vec_subcomm_in,*OP_Objects_subcomm->z,*OP_Objects_subcomm->theta,false,grad_nominal);
		  Nom_subcomm->Apply_Regularization_Update(z_vec_subcomm_out,z_vec_subcomm_in);
		  AQ->Write_Vector_to_Column(vecs_so_far_+k,z_vec_subcomm_out);		  
		  std::cout << "Completed matrix vector product " << k+1 << " out of " << num_samples << std::endl;
		}
	    }
	  
	  // Gather data in AQ
	  comm_->barrier();   
	  for(int k = 0; k < num_samples; k++)
	    {
	      char *buff = (char*)AQ->Get_Element_Ptr(0,vecs_so_far_+k);
	      comm_->broadcast(proc_dist->Get_Procs_Loop_Distribution()[k][0],8*nonzero_z_dim_,buff);
	    }
	  comm_->barrier();
	  
	  // Compute B
	  vecs_so_far_ += num_samples;
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	  for(int i = 0; i < vecs_so_far_; i++)
	    {
	      for(int j = 0; j < vecs_so_far_; j++)
	  	{
	  	  RealT val = 0.0;
	  	  for(int k = 0; k < nonzero_z_dim_; k++)
	  	    {
	  	      val += (*Q)(k,i)*(*AQ)(k,j);
	  	    }
	  	  B->Replace_Element(i,j,val);
	  	}
	    }
	  
	  HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(vecs_so_far_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_);
	  if(eig_alg_ == "Nystrom")
	    {
	      // Compute C = Chol(B)
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(B,C);	//B = C^T*C      

	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > F = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_);
	      for(int j = 0; j < vecs_so_far_; j++)
		{
		  // F(:,j) = AQ*C^{-1}*e_j
		  
		  // Compute C^{-1}*e_j
		  std::vector<RealT> c_inv = std::vector<RealT>(vecs_so_far_,0.0);
		  for(int k = vecs_so_far_-1; k >= 0; k--)
		    {
		      RealT val = 0.0;
		      if(k == j)
			{
			  val = 1.0;
			}
		      for(int l = vecs_so_far_-1; l > k; l--)
			{
			  val -= (*C)(k,l)*c_inv[l];
			}
		      c_inv[k] = val/(*C)(k,k);
		    }
		  
		  // F(:,j) = AQ*c_inv
		  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      RealT val = 0.0;
		      for(int k = 0; k < vecs_so_far_; k++)
			{
			  val += c_inv[k]*(*AQ)(i,k);
			}
		      F->Replace_Element(i,j,val);
		    }
		}
	      
	      // Need to compute the SVD of F
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::SVD<RealT>(F, U, VT, S); 

	      // Square singular values to get eigenvalues
	      for(int k = 0; k < vecs_so_far_; k++)
		{
		  S->Replace_Element(k,std::pow((*S)(k),2));
		}	      
	    }
	  else if(eig_alg_ == "Double Pass")
	    {
	      // Need to compute eigenvalue decomposition of B
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(B, V, S);
	      Q->Multiply(U,V); // Need to compute U = Q*V to get eigenvectors
	    }
	  
	  // Need to determine convergence and write solutions to text file
	  current_eig_ratio = (*S)(vecs_so_far_-1-oversampling_factor_)/(*S)(0);
	  
	  int k = vecs_so_far_-1-oversampling_factor_;
	  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_k = OP_Objects_subcomm->z->Clone();
	  HDSA::Ptr<HDSA::Vector<RealT> > z_matvec = OP_Objects_subcomm->z->Clone();

	  U->Write_Column_to_Vector(k,z_vec_k);

	  z_vec_k->scale(1.0/z_vec_k->norm());
	  OP_Objects_subcomm->rs_obj->hessVec_z_z(*z_matvec,*z_vec_k,*OP_Objects_subcomm->z,*OP_Objects_subcomm->theta,false,grad_nominal);
	  z_matvec->axpy(-(*S)(k),*z_vec_k);
	  RealT res_k = z_matvec->norm();
	  
	  RealT diff_up = std::abs((*S)(k)-(*S)(k-1));
	  RealT diff_down = std::abs((*S)(k)-(*S)(k+1)); 
	  RealT delta = std::min(diff_up,diff_down);
	  residual_bound_ = std::pow(res_k,2)/(delta*(*S)(k));
	  
	  if(comm_->getRank()==0)
	    {
	      std::cout << "Iteration number: " << iteration_count << " has computed " << vecs_so_far_ << " vectors and " << vecs_so_far_-oversampling_factor_ << " eigenvalues" << std::endl;
	      std::cout << "Smallest eigenvalue/largest eigenvalue = " << current_eig_ratio << std::endl;
	      std::cout << "Estimated residual bound for smallest eigenvalue = " << residual_bound_ << std::endl;
	      std::cout << "The eigenvalues computed thus far are: ";
	      for(int k = 0; k < vecs_so_far_-oversampling_factor_-1; k++)
		{
		  std::cout << (*S)(k) << " , ";
		}
	      std::cout << (*S)(vecs_so_far_-oversampling_factor_-1) << "  " << std::endl;;
	    }
	  
	  iteration_count += 1;
	  if( (current_eig_ratio < eig_ratio_tol_) & (residual_bound_ > residual_bound_tol_) )
	    {
	      oversampling_factor_ += proc_dist->Get_numSubComm();
	    }

	  if( ( (current_eig_ratio < eig_ratio_tol_) && (residual_bound_ < residual_bound_tol_) )  || (iteration_count > max_iter_) )
	    {
	      S_ = S;
	      U_ = U;
	    }
	  
	}
      
      RealT Time_alg = static_cast<RealT>(std::clock()-timer_alg)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Total time for eigenvalue computation: " << Time_alg << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
      
      Write_Solution();
    }

    void Write_Solution(void)
    {
      std::clock_t timer_SI = std::clock();
      // Estimate sensitivity indices
      int theta_dim = theta_->dimension();
      std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim,0.0);

      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
      OP_Objects->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = weight_matrices_factory_->Construct_Weight_Matrices(theta_,comm_); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
      // Precompute nominal gradient
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects->z->Clone();
      OP_Objects->rs_obj->gradient_z(*grad_nominal,*OP_Objects->z, *OP_Objects->theta);
     
      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec = OP_Objects->theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_matvec = OP_Objects->theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec = OP_Objects->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_k = OP_Objects->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_j = OP_Objects->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_out = OP_Objects->z->Clone();
      int r = vecs_so_far_-oversampling_factor_;
      RealT val = 0.0;
      RealT sk = 0.0;
      RealT sj = 0.0;
         
      for(int i = 0; i < theta_dim; i++)
	{
	  if(comm_->getRank() == 0)
	    {
	      std::cout << "Computing sensitivity index " << i+1 << " out of " << theta_dim << std::endl;
	    }
	  theta_vec->basis(i);
	  OP_Objects->rs_obj->hessVec_z_theta(*z_vec, *theta_vec, *OP_Objects->z, *OP_Objects->theta,true,grad_nominal);
	  val = 0.0;
	  for(int k = 0; k < r; k++)
	    {
	      U_->Write_Column_to_Vector(k,z_vec_k);
	      z_vec_k->scale(1.0/z_vec_k->norm());
	      sk = (z_vec_k->dot(*z_vec))/(*S_)(k);
	      weight_matrices->Apply_z_Weight_Mat(z_vec_out,z_vec_k);
	      for(int j = 0; j < r; j++)
		{
		  U_->Write_Column_to_Vector(j,z_vec_j);
		  z_vec_j->scale(1.0/z_vec_j->norm());
		  sj = (z_vec_j->dot(*z_vec))/(*S_)(j);
		  val += sk*sj*(z_vec_j->dot(*z_vec_out));
		}
	    }
	  theta_matvec->zero();
	  weight_matrices->Apply_theta_Weight_Mat(theta_matvec,theta_vec);
	  sensitivity_indices[i] = std::sqrt(val)/std::sqrt(theta_matvec->dot(*theta_vec));
	}
      RealT Time_SI = static_cast<RealT>(std::clock()-timer_SI)/static_cast<RealT>(CLOCKS_PER_SEC);      

      std::clock_t timer_write = std::clock();
      // Write to text files
      std::string name;
      std::ofstream fout;
      name = "Sensitivity_Indices_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < theta_dim; k++)
	{
	  fout << sensitivity_indices[k] << std::setw(20);
	}
      fout.close(); 

     // Write eigenvalues/vectors to text files
      name = "Hessian_evals_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < vecs_so_far_-oversampling_factor_; k++)
	{
	  fout << (*S_)(k) << std::setw(20);
	}
      fout.close();

      name = "Hessian_evecs_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int i  = 0; i < z_vec->dimension(); i++)
	{
	  if(z_vec->Is_entry_zero(i))
	    {
	      for(int k = 0; k < vecs_so_far_-oversampling_factor_; k++)
		{   
		  fout << 0.0 << std::setw(20);
		}
	    }
	  else
	    {
	      for(int k = 0; k < vecs_so_far_-oversampling_factor_; k++)
		{   
		  fout << (*U_)(z_vec->Get_map_full_to_reduced(i),k) << std::setw(20);
		}
	    }
	  fout << " " << std::endl;
	}
      fout.close();

      RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);     
      
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Time estimating sensitivity indices from hessian decomposition: " << Time_SI << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
   

    }


 
  };

}

#endif 
