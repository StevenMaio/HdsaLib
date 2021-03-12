#ifndef HDSA_RANDOMIZED_LIKELIHOOD_INFORMED_SUBSPACE_EVP_HPP
#define HDSA_RANDOMIZED_LIKELIHOOD_INFORMED_SUBSPACE_EVP_HPP

namespace HDSA
{

  template <class RealT>
  class Randomized_Likelihood_Informed_Subspace_EVP{
  private:
    std::string eig_alg_;
    int oversampling_factor_;
    int max_iter_;
    RealT eig_min_;
    int initial_samples_;
    int nonzero_z_dim_;
    RealT reg_hess_inv_tol_;
    bool reg_hess_inv_verbose_;
    std::string reg_hess_inv_solver_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    int vecs_so_far_;
    int numProcs_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_; // Matrix of eigenvectors
    HDSA::Ptr<HDSA::Vector<RealT> > S_; // vector of eigenvalues
 
  public:
  
    Randomized_Likelihood_Informed_Subspace_EVP(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
						const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
						const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), 
      weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    {
      eig_alg_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Eigenvalue Algorithm", "Nystrom");
      int target_rank = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Initial Target Rank", 1);
      oversampling_factor_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Oversampling Factor", 10);
      max_iter_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Maximum Number of Iterations", 2);
      eig_min_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Minimum Eigenvalue Tolerance", 1.0);
      reg_hess_inv_tol_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Regularization Hessian Inversion Tolerance", 1.e-12);
      reg_hess_inv_verbose_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Regularization Hessian Inversion Verbosity", false);
      reg_hess_inv_solver_ = parlist_sensitivity->sublist("Randomized Likelihood Informed Subspace EVP").get("Regularization Hessian Inversion Solver", "CG");
      initial_samples_ = target_rank + oversampling_factor_;
      numProcs_ = comm->getSize();
      if( (initial_samples_%numProcs_) !=0 )
	{
	  std::cout << "Number of samples does not evenly divide the number of processors" << std::endl;
	}
      if( numProcs_ > initial_samples_ )
	{
	  std::cout << "The code does not support using more processors than samples" << std::endl;
	}
    }
    
    void Compute(void)
    {
      // This function executes the algorithm
      std::clock_t timer_alg = std::clock();

      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Starting randomized LIS eigensolver set up" << std::endl;
	  std::cout << " " << std::endl;
	} 
      std::clock_t timer_setup = std::clock();
      
      int r = std::min(initial_samples_,numProcs_);
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,r);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
  
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
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
      RealT current_eig_min = eig_min_ + 1.0;
      
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

      // Precompute nominal gradient
      HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal = OP_Objects_subcomm->z->Clone();
      OP_Objects_subcomm->rs_obj->gradient_z(*grad_nominal,*OP_Objects_subcomm->z, *OP_Objects_subcomm->theta);
      Nom_subcomm->Set_grad_nominal(grad_nominal);

      HDSA::Ptr<HDSA::Linear_Operator<RealT> > reg_op = HDSA::makePtr<Reg_Operator<RealT> >(OP_Objects_subcomm,Nom_subcomm,grad_nominal);

      comm_->barrier();
      RealT Time_setup = static_cast<RealT>(std::clock()-timer_setup)/static_cast<RealT>(CLOCKS_PER_SEC);
         if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Completed randomized LIS eigensolver set up in " << Time_setup << " seconds" << std::endl;
	  std::cout << " " << std::endl;
	} 

      while( (current_eig_min > eig_min_) && (iteration_count <= max_iter_) )
	{
	  int num_samples = proc_dist->Get_numSubComm();
	  int num_loops = 1;
	  if(iteration_count==1)
	    {
	      num_samples = initial_samples_;
	      num_loops = num_samples/numProcs_;
	    }
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,num_samples);
	  // Need to perform sketch of the range of Y
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Beginning range finder sketch on iteration number " << iteration_count << " and sample number " << sample_index_ << std::endl;
	      std::cout << " " << std::endl;
	    }      
	  	 
	  comm_->barrier();
	  for(int i = 0; i < num_loops; i++)
	    {
	      for(int j = 0; j < numProcs_; j++)
		{
		  int k = i*numProcs_ + j;
		  if(proc_dist->Does_Processor_Own_Vector(j))
		    {
		      std::clock_t timer_matvec = std::clock();
		      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
		      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_random = OP_Objects_subcomm->z->Clone();
		      // Populate vectors with standard normal samples
		      for(int l = 0; l < nonzero_z_dim_; l++)
			{
			  z_vec_random->Replace_Element(z_vec_random->Get_map_reduced_to_full(l),distribution(generator));
			}
		      OP_Objects_subcomm->rs_obj->Misfit_hessVec_z_z(*z_vec_subcomm,*z_vec_random,*OP_Objects_subcomm->z,*OP_Objects_subcomm->theta,false,grad_nominal);
		      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(z_vec_random,z_vec_subcomm,reg_op,reg_hess_inv_tol_,reg_hess_inv_solver_,reg_hess_inv_verbose_);
		      Y->Write_Vector_to_Column(k,z_vec_random);
		      RealT Time_matvec = static_cast<RealT>(std::clock()-timer_matvec)/static_cast<RealT>(CLOCKS_PER_SEC);
		      std::cout << "Completed matrix vector product " << k+1 << " out of " << num_samples << " in " << Time_matvec << " seconds" << std::endl;
		    }
		}
	    }

	  comm_->barrier();  
	  std::clock_t timer_broadcast1 = std::clock();
	  for(int i = 0; i < num_loops; i++)
	    {
	      for(int j = 0; j < numProcs_; j++)
		{
		  int k = i*numProcs_ + j; 
		  char *buff = (char*)Y->Get_Element_Ptr(0,k);
		  comm_->broadcast(proc_dist->Get_Procs_Loop_Distribution()[j][0],8*nonzero_z_dim_,buff);
		}
	    }
	  comm_->barrier();
	  RealT Time_broadcast1 = static_cast<RealT>(std::clock()-timer_broadcast1)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Completed the first data broadcast in " << Time_broadcast1 << " seconds" << std::endl;
	      std::cout << " " << std::endl;
	    }  
	  
	  comm_->barrier();  
	  std::clock_t timer_denseLA1 = std::clock(); 
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
		  normalize += std::pow((*Y)(j,k),2);
		}
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  QI->Replace_Element(i,vecs_so_far_+k,(*Y)(i,k)/std::sqrt(normalize));
		}
	    }
 
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_+num_samples);
	  HDSA::Linear_Algebra::CholQR<RealT>(QI, Z, reg_op, OP_Objects_subcomm->z);
	
	  for(int k = 0; k < num_samples; k++)
	    {
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  Q->Replace_Element(i,vecs_so_far_+k,(*Z)(i,vecs_so_far_+k));
		}
	    }

	  comm_->barrier();
	  RealT Time_denseLA1 = static_cast<RealT>(std::clock()-timer_denseLA1)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Completed the first dense linear algebra in " << Time_denseLA1 << " seconds" << std::endl;
	      std::cout << " " << std::endl;
	    }  
	  
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Beginning eigenvalue estimation on iteration number " << iteration_count << " and sample number " << sample_index_ << std::endl;
	      std::cout << " " << std::endl;
	    } 
	  
	  // Compute A*Q for new columns of Q
	  comm_->barrier();
	  for(int i = 0; i < num_loops; i++)
	    {
	      for(int j = 0; j < numProcs_; j++)
		{
		  int k = i*numProcs_ + j;
		  if(proc_dist->Does_Processor_Own_Vector(j))
		    {
		      std::clock_t timer_matvec = std::clock();
		      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_in = OP_Objects_subcomm->z->Clone();
		      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_out = OP_Objects_subcomm->z->Clone();
		      
		      for(int i = 0; i < nonzero_z_dim_; i++)
			{
			  z_vec_subcomm_in->Replace_Element(z_vec_subcomm_in->Get_map_reduced_to_full(i),(*Q)(i,vecs_so_far_+k));
			}  

		      OP_Objects_subcomm->rs_obj->Misfit_hessVec_z_z(*z_vec_subcomm_out,*z_vec_subcomm_in,*OP_Objects_subcomm->z,*OP_Objects_subcomm->theta,false,grad_nominal);
		      AQ->Write_Vector_to_Column(vecs_so_far_+k,z_vec_subcomm_out);
		      RealT Time_matvec = static_cast<RealT>(std::clock()-timer_matvec)/static_cast<RealT>(CLOCKS_PER_SEC);
		      std::cout << "Completed matrix vector product " << k+1 << " out of " << num_samples << " in " << Time_matvec << " seconds" <<  std::endl;
		    }
		}
	    }
	  
	  // Gather data in AQ
	  comm_->barrier();   
	  std::clock_t timer_broadcast2 = std::clock(); 
	  for(int i = 0; i < num_loops; i++)
	    {
	      for(int j = 0; j < numProcs_; j++)
		{
		  int k = i*numProcs_ + j; 
		  char *buff = (char*)AQ->Get_Element_Ptr(0,vecs_so_far_+k);
		  comm_->broadcast(proc_dist->Get_Procs_Loop_Distribution()[j][0],8*nonzero_z_dim_,buff);
		}
	    }
	  comm_->barrier();
	  RealT Time_broadcast2 = static_cast<RealT>(std::clock()-timer_broadcast2)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Completed the second data broadcast in " << Time_broadcast2 << " seconds" << std::endl;
	      std::cout << " " << std::endl;
	    }  
	  
	  comm_->barrier();  
	  std::clock_t timer_denseLA2 = std::clock(); 
	  // Compute T
	  vecs_so_far_ += num_samples;
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > T = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	  for(int i = 0; i < vecs_so_far_; i++)
	    {
	      for(int j = 0; j < i; j++)
		{
		  T->Replace_Element(i,j,(*T)(j,i));
		}
	      for(int j = i; j < vecs_so_far_; j++)
		{
		  RealT val = 0.0;
		  for(int k = 0; k < nonzero_z_dim_; k++)
		    {
		      val += (*Q)(k,i)*(*AQ)(k,j);
		    }
		  T->Replace_Element(i,j,val);
		}
	    }
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_); // Matrix of eigenvectors
	  HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(vecs_so_far_); // Vector of eigenvalues
	  if(eig_alg_ == "Nystrom")
	    {
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      // Compute Cholesky factorization of T, T=C^T*C
	      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(T, C);
	      // Compute M = AQ*C^{-1}
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > M = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_);
	      for(int k = 0; k < vecs_so_far_; k++)
		{
		  HDSA::Ptr<HDSA::Vector<RealT> > ek = HDSA::makePtr<Std_Vector<RealT> >(vecs_so_far_);
		  ek->basis(k);
		  HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(vecs_so_far_);
		  HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>( x, ek, C);

		  HDSA::Ptr<HDSA::Vector<RealT> > AQx = HDSA::makePtr<Std_Vector<RealT> >(nonzero_z_dim_);
		  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      RealT val = 0.0;
		      for(int l = 0; l < vecs_so_far_; l++)
			{
			  val += (*AQ)(i,l)*(*x)(l);
			}
		      AQx->Replace_Element(i,val);
		    }
		  M->Write_Vector_to_Column(k, AQx);
		}
	      // Compute CholQR of M with reg_op inverse inner product
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > QM = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_);
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > QM_hat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,vecs_so_far_);
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > RM = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::CholQR<RealT>(M, QM, reg_op, OP_Objects_subcomm->z, QM_hat, RM, true);
	      // Compute SVD of RM
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > UR = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VRT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::SVD<RealT>(RM,UR,VRT,S);
	      // Square singular values and multiply QM*UR to get eigenvectors
	      for(int k = 0; k < vecs_so_far_; k++)
		{
		  S->Replace_Element(k,std::pow((*S)(k),2));
		  QM_hat->Multiply(U,UR);
		}
	    }
	  else if(eig_alg_ == "Double Pass")
	    {
	      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(vecs_so_far_,vecs_so_far_);
	      HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(T,V,S); 
  	      
	      // Need to compute U = Q*V to get eigenvectors
	      for(int i = 0; i < nonzero_z_dim_; i++)
		{
		  for(int j = 0; j < vecs_so_far_; j++)
		    {
		      RealT val = 0.0;
		      for(int k = 0; k < vecs_so_far_; k++)
			{
			  val += (*Q)(i,k)*(*V)(k,j);
			}
		      U->Replace_Element(i,j,val);
		    }
		}
	    }
	  
	  comm_->barrier();
	  RealT Time_denseLA2 = static_cast<RealT>(std::clock()-timer_denseLA2)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << " " << std::endl;
	      std::cout << "Completed the second dense linear algebra in " << Time_denseLA2 << " seconds" << std::endl;
	      std::cout << " " << std::endl;
	    } 

	  // Need to determine convergence and write solutions to text file
	  current_eig_min = (*S)(vecs_so_far_-1-oversampling_factor_);
	  
	  if(comm_->getRank()==0)
	    {
	      std::cout << "Iteration number: " << iteration_count << " has computed " << vecs_so_far_ << " vectors and " << vecs_so_far_-oversampling_factor_ << " eigenvalues" << std::endl;
	      std::cout << "Smallest eigenvalue = " << current_eig_min << std::endl;
	      std::cout << "The eigenvalues computed thus far are: ";
	      for(int k = 0; k < vecs_so_far_-oversampling_factor_-1; k++)
		{
		  std::cout << (*S)(k) << " , ";
		}
	      std::cout << (*S)(vecs_so_far_-oversampling_factor_-1) << "  " << std::endl;;
	    }
	  
	  iteration_count += 1;

	  if( (current_eig_min < eig_min_)  || (iteration_count > max_iter_) )
	    {
	      int n = S->dimension();
	      S_ = HDSA::makePtr<Std_Vector<RealT> >(n);
	      S_->set(*S);
	      U_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,n);
	      for(int j = 0; j < n; j++)
		{ 
		  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      U_->Replace_Element(i,j,(*U)(i,j));
		    }
		}
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

    void Construct_B_Parallel(void)
    {
      // This function executes the algorithm
      std::clock_t timer_alg = std::clock();
      
      int theta_dim = theta_->dimension();
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,theta_dim);
      HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
  
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      nonzero_z_dim_ = Nom_subcomm->Get_nonzero_z_dim();
      HDSA::Ptr<HDSA::Sensitivity_Operators_RS<RealT> > Sen_Op_subcomm 
	= HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,theta_dim);
      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec = theta_->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec = OP_Objects_subcomm->z->Clone();
      for(int k = 0; k < theta_dim; k++)
	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
	    {
	      std::clock_t timer_matvec = std::clock();
	      theta_vec->basis(k);
	      Sen_Op_subcomm->Apply_B(z_vec,theta_vec);
	      B->Write_Vector_to_Column(k,z_vec);
	      RealT Timer_matvec = static_cast<RealT>(std::clock()-timer_matvec)/static_cast<RealT>(CLOCKS_PER_SEC);
	      if(subcomm->getRank() == 0)
      		{
      		  std::cout << "Computed column " << k+1 << " out of " << theta_dim << " of B in " << Timer_matvec << " seconds" << std::endl;
      		}
	    }
	}
      comm_->barrier();
      proc_dist->Broadcast_Matrix(B);
      comm_->barrier();
      RealT Timer_alg = static_cast<RealT>(std::clock()-timer_alg)/static_cast<RealT>(CLOCKS_PER_SEC);
      if(comm_->getRank() == 0)
	{
	  std::cout << "Completed construction of B in " << Timer_alg << " seconds" << std::endl;
	}

      // Write solutions to text files
      if(comm_->getRank() == 0)
	{
	  RealT timer_write = std::clock();
	  std::string name;
	  std::ofstream fout;
	  name = "B.txt";
	  fout.open(name);
	  for(int i = 0; i < nonzero_z_dim_; i++)
	    {
	      for(int j = 0; j < theta_dim; j++)
		{
		  fout << std::setprecision(16) << (*B)(i,j) << "  ";
		}
	      fout << "  " << std::endl;
	    }
	  fout.close();
	  RealT Timer_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);
	  std::cout << "Wrote B to a text file in " << Timer_write << " seconds" << std::endl;
	}

    }

    void Write_Solution(void)
    {
      std::clock_t timer_normalize = std::clock();
      
      int r = vecs_so_far_-oversampling_factor_;
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
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > reg_op = HDSA::makePtr<Reg_Operator<RealT> >(OP_Objects,Nom,grad_nominal);
 
      std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > z_vecs;
      z_vecs.resize(r);
      std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > z_matvecs;
      z_matvecs.resize(r);
      for(int k = 0; k < r; k++)
	{
	  z_vecs[k] = OP_Objects->z->Clone();
	  z_matvecs[k] = OP_Objects->z->Clone();
	  U_->Write_Column_to_Vector(k,z_vecs[k]);
	  reg_op->matvec(z_matvecs[k],z_vecs[k]);
	  z_vecs[k]->scale(1.0/std::sqrt(z_vecs[k]->dot(*z_matvecs[k])));
	  z_matvecs[k]->zero();
	  weight_matrices->Apply_z_Weight_Mat(z_matvecs[k],z_vecs[k]);
	}

      RealT Time_normalize = static_cast<RealT>(std::clock()-timer_normalize)/static_cast<RealT>(CLOCKS_PER_SEC); 
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Time normalizing: " << Time_normalize << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
     
      std::clock_t timer_write = std::clock();

      // Write to text files
      std::string name;
      std::ofstream fout;
	  
      // Write eigenvalues/vectors to text files
      name = "LIS_evals_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < r; k++)
	{
	  fout << (*S_)(k) << std::setw(20);
	}
      fout.close();
	  
      name = "LIS_evecs_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int i = 0; i < z_vecs[0]->dimension(); i++)
	{
	  for(int j = 0; j < r; j++)
	    {
	      fout << (*z_vecs[j])(i) << std::setw(20);
	    }
	  fout << "  " << std::endl;
	}
      fout.close();
      
      name = "eigenvector_weightmat_products_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int i = 0; i < r; i++)
	{
	  for(int j = 0; j < r; j++)
	    {
	      fout << z_vecs[j]->dot(*z_matvecs[i]) << std::setw(20);
	    }
	  fout << "  " << std::endl;
	}
      fout.close();   
	  
      RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);  

      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
   
    }

    // Overload HDSA::Linear_Operator to take matrix vector products for the regularization hessian solve
    template <class ScalarType>
    class Reg_Operator : public HDSA::Linear_Operator<ScalarType>
    {
      HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > OP_Objects_;
      HDSA::Ptr<HDSA::Nominal_Data<ScalarType> > Nom_;
      HDSA::Ptr<HDSA::Vector<ScalarType> > grad_;
      ScalarType coeff_;
    
    public:
      
      Reg_Operator(const HDSA::Ptr<HDSA::Opt_Problem_Objects<ScalarType> > & OP_Objects, const HDSA::Ptr<HDSA::Nominal_Data<ScalarType> > & Nom, 
		   const HDSA::Ptr<HDSA::Vector<ScalarType> > & grad_nominal): OP_Objects_(OP_Objects), Nom_(Nom), grad_(grad_nominal)
      { 
	Nom_->Set_grad_nominal(grad_);
      }
    
      //! Dtor
      ~Reg_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const 
      {
	OP_Objects_->rs_obj->Regularization_hessVec_z_z(*y,*x,*OP_Objects_->z,*OP_Objects_->theta,false, grad_);
	Nom_->Apply_Regularization_Update(y,x);
      }
      
    };
 
  };

}

#endif 
