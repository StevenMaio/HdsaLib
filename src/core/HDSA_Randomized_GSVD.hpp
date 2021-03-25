// #ifndef HDSA_RANDOMIZED_GSVD_HPP
// #define HDSA_RANDOMIZED_GSVD_HPP

// // This class executes the randomized GSVD solver

// namespace HDSA
// {

//   template <class RealT>
//   class Randomized_GSVD
//   {
//   private:
//     bool reduced_space_sen_;
//     int num_sing_vals_, p_, kpp_, q_, nonzero_z_dim_, theta_dim_;
//     HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
//     HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
//     HDSA::Ptr<const HDSA::Comm<int> > comm_;
//     HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
//     HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
//     int sample_index_;
//     HDSA::Ptr<HDSA::Vector<RealT> > S_;
//     HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_;
//     HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V_;

//   public:
  
//     Randomized_GSVD(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
// 		    const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
// 		    const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
//       theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
//     {
//       reduced_space_sen_ = parlist_sensitivity->sublist("Formulation").get("Reduced Space Sensitivities",true);
//       num_sing_vals_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Singular Values",1); 
//       p_ = parlist_sensitivity->sublist("GSVD Solver").get("Oversampling Factor",20);
//       q_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Subspace Iterations",0);
//       kpp_ = num_sing_vals_ + p_;
//     }
    
//     void Compute(void)
//     {
//       // This function executes the randomized GSVD
//       std::clock_t timer_gsvd = std::clock();

//       HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,kpp_);
//       HDSA::Ptr<HDSA::Comm<int> > subcomm = comm_->createSubcommunicator(proc_dist->Get_Comm_Split_Ranks());
//       HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_subcomm = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,subcomm);
//       OP_Objects_subcomm->Load_Optimal_Solution();
//       bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
//       if(enforce_z_zeros_)
// 	{
// 	  OP_Objects_subcomm->z->Enforce_Zeros();
// 	}
//       theta_dim_ = OP_Objects_subcomm->theta->dimension();
//       HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
//       HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
//       nonzero_z_dim_ = Nom_subcomm->Get_nonzero_z_dim();
//       HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm;
//       if(reduced_space_sen_)
//       	{
//       	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
//       	}
//       else
//       	{
// 	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
//       	}
//       // Need to perform sketch of the range of D  
//       if(comm_->getRank() == 0)
//       	{
//       	  std::cout << " " << std::endl;
//       	  std::cout << "Beginning range sketching for local sensitivity number " << sample_index_ << std::endl;
//       	  std::cout << " " << std::endl;
//       	}
  
//       // Instantiate matrix to interface with LA routines   
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);   
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WY = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 

//       comm_->barrier();
//       std::clock_t timer_Loop_1 = std::clock();
//       for(int k = 0; k < kpp_; k++)
//       	{
// 	  if(proc_dist->Does_Processor_Own_Vector(k))
//       	    {
// 	      std::clock_t timer_range_iter_k = std::clock();
// 	      HDSA::Ptr<HDSA::Vector<RealT> >  theta_vec_random = Sen_Op_subcomm->Generate_Random_theta_Vector();
//       	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
//       	      Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm,theta_vec_random);
// 	      Y->Write_Vector_to_Column(k,z_vec_subcomm);

// 	      HDSA::Ptr<HDSA::Vector<RealT> > W_z_vec_subcomm = OP_Objects_subcomm->z->Clone();
// 	      weight_matrices_subcomm->Apply_z_Weight_Mat(W_z_vec_subcomm,z_vec_subcomm);
// 	      WY->Write_Vector_to_Column(k,W_z_vec_subcomm);
	      
// 	      RealT Time_range_iter_k = static_cast<RealT>(std::clock()-timer_range_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
//       	      if(subcomm->getRank() == 0)
//       		{
//       		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the range sketching loop for local sensitivity number " << sample_index_ << " in " << Time_range_iter_k << " seconds" << std::endl;
//       		}
//       	    }
//       	}
  
//       comm_->barrier();
//       proc_dist->Broadcast_Matrix(Y);
//       proc_dist->Broadcast_Matrix(WY);
//       RealT Time_Loop_1 = static_cast<RealT>(std::clock()-timer_Loop_1)/static_cast<RealT>(CLOCKS_PER_SEC);
      
//       // Need to orthogonalize the sketch with respect to the z weight matrix
//       std::clock_t timer_ortho_z = std::clock();

//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > SQ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_); 
//       HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Y, WY, Q, R, SQ);

//       comm_->barrier();
//       RealT Time_ortho_z = static_cast<RealT>(std::clock()-timer_ortho_z)/static_cast<RealT>(CLOCKS_PER_SEC);
//   	  if(comm_->getRank() == 0)
// 	    {
// 	      std::cout << "Completed orthogonalization for local sensitivity number " << sample_index_<< " in " << Time_ortho_z << " seconds" << std::endl;
// 	    }
      
//       // Subspace iteration
//       std::clock_t timer_subspace_iter = std::clock();
      
//       if((comm_->getRank() == 0) & (q_ > 0) )
//       	{
//       	  std::cout << " " << std::endl;
//       	  std::cout << "Beginning subspace iteration for local sensitivity number " << sample_index_ << std::endl;
//       	  std::cout << " " << std::endl;
//       	}

//       for(int j = 0; j < q_; j++)
//       	{
// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_); 	  
// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WY_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);

//       	  for(int k = 0; k < kpp_; k++)
//       	    {
// 	      if(proc_dist->Does_Processor_Own_Vector(k))
//       		{
// 		  std::clock_t timer_subspace_iter_k = std::clock();
// 		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
// 		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  
// 		  SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);
// 		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm, z_vec_subcomm_1);		  

// 		  Y_subspace_iter_param->Write_Vector_to_Column(k,theta_vec_subcomm);

// 		  HDSA::Ptr<HDSA::Vector<RealT> > W_theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
// 		  weight_matrices_subcomm->Apply_theta_Weight_Mat_Inverse(W_theta_vec_subcomm,theta_vec_subcomm);
// 		  WY_subspace_iter_param->Write_Vector_to_Column(k,W_theta_vec_subcomm);

// 		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
// 		  if(subcomm->getRank() == 0)
// 		    {
// 		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the first loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
// 		    }
// 		}
//       	    }
// 	  comm_->barrier();
// 	  proc_dist->Broadcast_Matrix(Y_subspace_iter_param);
// 	  proc_dist->Broadcast_Matrix(WY_subspace_iter_param);
	  
// 	  std::clock_t timer_ortho = std::clock();

// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_Q_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
// 	  HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Y_subspace_iter_param, WY_subspace_iter_param, Q_subspace_iter_param, R,  W_Q_subspace_iter_param);

//       	  comm_->barrier();
// 	  RealT Time_ortho = static_cast<RealT>(std::clock()-timer_ortho)/static_cast<RealT>(CLOCKS_PER_SEC);
// 	  if(comm_->getRank() == 0)
// 	    {
// 	      std::cout << "Completed parameter orthogonalization for local sensitivity number " << sample_index_<< " in " << Time_ortho << " seconds" << std::endl;
// 	    }
	  
// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_subspace_iter_opt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
// 	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WY_subspace_iter_opt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
//       	  for(int k = 0; k < kpp_; k++)
//       	    {
// 	      if(proc_dist->Does_Processor_Own_Vector(k))
//       		{
// 		  std::clock_t timer_subspace_iter_k = std::clock();
// 		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
//       		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
 
// 		  W_Q_subspace_iter_param->Write_Column_to_Vector(k,theta_vec_subcomm);
//       		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm, theta_vec_subcomm);
// 		  Y_subspace_iter_opt->Write_Vector_to_Column(k,z_vec_subcomm);

// 		  HDSA::Ptr<HDSA::Vector<RealT> > W_z_vec_subcomm = OP_Objects_subcomm->z->Clone();
// 		  weight_matrices_subcomm->Apply_z_Weight_Mat(W_z_vec_subcomm,z_vec_subcomm);
// 		  WY_subspace_iter_opt->Write_Vector_to_Column(k,W_z_vec_subcomm);

// 		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
// 		  if(subcomm->getRank() == 0)
// 		    {
// 		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the second loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
// 		    }
//       		}
//       	    }
//       	  comm_->barrier();
// 	  proc_dist->Broadcast_Matrix(Y_subspace_iter_opt);   
// 	  proc_dist->Broadcast_Matrix(WY_subspace_iter_opt);

// 	  HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Y_subspace_iter_opt, WY_subspace_iter_opt, Q, R, SQ);
//       	  comm_->barrier();
	  
//       	  if(comm_->getRank() == 0)
//       	    {
// 	      std::cout << " " << std::endl;
//       	      std::cout << "Completed subspace iteration " << j+1 << " out of " << q_  << " for local sensitivity number " << sample_index_<< std::endl;
//       	    }
	  
//       	}
  
//       comm_->barrier();
//       RealT Time_Loop_SubIter = static_cast<RealT>(std::clock()-timer_subspace_iter)/static_cast<RealT>(CLOCKS_PER_SEC);
      
//       // Need to project D onto the sketched subspace
//       if(comm_->getRank() == 0)
//       	{
// 	  std::cout << " " << std::endl;
//       	  std::cout << "Beginning projection onto the sketched subspace for local sensitivity number " << sample_index_ << std::endl;
//       	  std::cout << " " << std::endl;
//       	}
      
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);     
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > WQ_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_); 

//       comm_->barrier();
//       std::clock_t timer_Loop_2 = std::clock();
//       for(int k = 0; k < kpp_; k++)
//       	{
// 	  if(proc_dist->Does_Processor_Own_Vector(k))
//       	    { 
// 	      std::clock_t timer_proj_iter_k = std::clock();
// 	      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm_1 = OP_Objects_subcomm->theta->Clone();
// 	      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm_2 = OP_Objects_subcomm->theta->Clone();
// 	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();

// 	      SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);
// 	      Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm_1, z_vec_subcomm_1);
//       	      weight_matrices_subcomm->Apply_theta_Weight_Mat_Inverse(theta_vec_subcomm_2,theta_vec_subcomm_1);
	      
// 	      Q_proj->Write_Vector_to_Column(k,theta_vec_subcomm_2);
// 	      WQ_proj->Write_Vector_to_Column(k,theta_vec_subcomm_1);
	      
// 	      RealT Time_proj_iter_k = static_cast<RealT>(std::clock()-timer_proj_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
//       	      if(subcomm->getRank() == 0)
//       		{
//       		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the projection loop for local sensitivity number " << sample_index_ << " in " << Time_proj_iter_k << " seconds" << std::endl;
//       		}
//       	    }
//       	}
      
//       comm_->barrier();
//       proc_dist->Broadcast_Matrix(Q_proj);
//       proc_dist->Broadcast_Matrix(WQ_proj);
//       RealT Time_Loop_2 = static_cast<RealT>(std::clock()-timer_Loop_2)/static_cast<RealT>(CLOCKS_PER_SEC);
  
//       // Need to orthogonalize the projection with respect to the parameter mass matrix
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);  
//       HDSA::Linear_Algebra::CholQR_Pre_W<RealT>(Q_proj, WQ_proj, Q_param, R);


//       // Need to compute the SVD of R
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
//       HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
//       S_ = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
//       HDSA::Linear_Algebra::SVD<RealT>(R, U, VT, S_); 
//       U_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
//       V_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
//       Q->Multiply(U_,VT,false,true);
//       Q_param->Multiply(V_,U);

//       // The columns of U contain the approximate opt singular vectors
//       // The columns of V contain the approximate parameter singular vectors
        
//       RealT Time_gsvd = static_cast<RealT>(std::clock()-timer_gsvd)/static_cast<RealT>(CLOCKS_PER_SEC);
      
//       if(comm_->getRank() == 0)
// 	{
// 	  std::cout << " " << std::endl;
// 	  std::cout << "Total time for GSVD computation: " << Time_gsvd << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << "Time contributed by the first loop: " << Time_Loop_1 << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << "Time contributed by the subspace iteration: " << Time_Loop_SubIter << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << "Time contributed by the second loop: " << Time_Loop_2 << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << "Time contributed by everything else: " << Time_gsvd-Time_Loop_1-Time_Loop_2-Time_Loop_SubIter << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << " " << std::endl;
// 	}

//       if( OP_Objects_subcomm->theta->dimension() == theta_->dimension() )
// 	{
// 	  Write_Solution(OP_Objects_subcomm->z,weight_matrices_subcomm);
// 	}
//       else
// 	{
// 	  Write_Solution_Model_Error(OP_Objects_subcomm->z);
// 	}
//     }

//     void Write_Solution(HDSA::Ptr<HDSA::Vector<RealT> > & z, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices)
//     {
//       std::clock_t timer_SI = std::clock();
//       // Estimate sensitivity indices
//       std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_vecs;
//       theta_vecs.resize(num_sing_vals_);
//       HDSA::Ptr<HDSA::Vector<RealT> > theta = theta_->Clone();
//       for(int k = 0; k < num_sing_vals_; k++)
// 	{
// 	  theta->zero();
// 	  for(int i = 0; i < theta_dim_; i++)
// 	    {
// 	      theta->Replace_Element(i,(*V_)(i,k));
// 	    }
// 	  theta_vecs[k] = theta_->Clone();
// 	  weight_matrices->Apply_theta_Weight_Mat(theta_vecs[k],theta);
// 	}
      
//       HDSA::Ptr<HDSA::Vector<RealT> > theta_matvec = theta_->Clone();
//       std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim_,0.0);

//       for(int i = 0; i < theta_dim_; i++)
// 	{
// 	  theta->basis(i);
// 	  weight_matrices->Apply_theta_Weight_Mat(theta_matvec,theta);
// 	  for(int k = 0; k < num_sing_vals_; k++)
// 	    {
// 	      sensitivity_indices[i] += std::pow((*S_)(k),2)*std::pow((*theta_vecs[k])(i),2);
// 	    }
// 	  sensitivity_indices[i] = std::sqrt(sensitivity_indices[i])/std::sqrt((*theta_matvec)(i));
// 	}
//       RealT Time_SI = static_cast<RealT>(std::clock()-timer_SI)/static_cast<RealT>(CLOCKS_PER_SEC);      

//       std::clock_t timer_write = std::clock();
//       // Write to text files
//       std::string name;
//       std::ofstream fout;
//       name = "Sensitivity_Indices_" + std::to_string(sample_index_) + ".txt";
//       fout.open(name);
//       for(int k = 0; k < theta_dim_; k++)
// 	{
// 	  fout << sensitivity_indices[k] << std::setw(20);
// 	}
//       fout.close();

//       name = "z_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
//       fout.open(name);
//       for(int i  = 0; i < z->dimension(); i++)
// 	{
// 	  if(z->Is_entry_zero(i))
// 	    {
// 	      for(int k = 0; k < num_sing_vals_; k++)
// 		{   
// 		  fout << 0.0 << std::setw(20);
// 		}
// 	    }
// 	  else
// 	    {
// 	      for(int k = 0; k < num_sing_vals_; k++)
// 		{   
// 		  fout << (*U_)(z->Get_map_full_to_reduced(i),k) << std::setw(20);
// 		}
// 	    }
// 	  fout << " " << std::endl;
// 	}
//       fout.close();

//       name = "theta_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
//       fout.open(name);
//       for(int i  = 0; i < theta_dim_; i++)
// 	{
// 	  for(int k = 0; k < num_sing_vals_; k++)
// 	    {    
// 	      fout << (*V_)(i,k) << std::setw(20);
// 	    }
// 	  fout << " " << std::endl;
// 	}
//       fout.close();
	    
//       name = "Singular_Values_" + std::to_string(sample_index_) + ".txt";
//       fout.open(name);
//       for(int k = 0; k < num_sing_vals_; k++)
// 	{
// 	  fout << (*S_)(k) << std::setw(20);
// 	}
//       fout.close();

//       name = "theta_Sample_" + std::to_string(sample_index_) + ".txt";
//       fout.open(name);
//       for(int k = 0; k < theta_dim_; k++)
// 	{
// 	  fout << (*theta_)(k) << std::setw(20);
// 	}
//       fout.close();
    
//       RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);     

//       if(comm_->getRank() == 0)
// 	{
// 	  std::cout << " " << std::endl;
// 	  std::cout << "Time estimating sensitivity indices from GSVD: " << Time_SI << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << " " << std::endl;
// 	}

//     }

//     void Write_Solution_Model_Error(HDSA::Ptr<HDSA::Vector<RealT> > & z)
//     {
//       std::clock_t timer_write = std::clock();
//       if(comm_->getRank() == 0)
// 	{
// 	  // Write to text files
// 	  std::string name;
// 	  std::ofstream fout;

// 	  name = "z_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
// 	  fout.open(name);
// 	  for(int i  = 0; i < z->dimension(); i++)
// 	    {
// 	      if(z->Is_entry_zero(i))
// 		{
// 		  for(int k = 0; k < num_sing_vals_; k++)
// 		    {   
// 		      fout << 0.0 << std::setw(20);
// 		    }
// 		}
// 	      else
// 		{
// 		  for(int k = 0; k < num_sing_vals_; k++)
// 		    {   
// 		      fout << (*U_)(z->Get_map_full_to_reduced(i),k) << std::setw(20);
// 		    }
// 		}
// 	      fout << " " << std::endl;
// 	    }
// 	  fout.close();
	  
// 	  name = "theta_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
// 	  fout.open(name);
// 	  for(int i  = 0; i < theta_dim_; i++)
// 	    {
// 	      for(int k = 0; k < num_sing_vals_; k++)
// 		{    
// 		  fout << (*V_)(i,k) << std::setw(20);
// 		}
// 	      fout << " " << std::endl;
// 	    }
// 	  fout.close();
	  
// 	  name = "Singular_Values_" + std::to_string(sample_index_) + ".txt";
// 	  fout.open(name);
// 	  for(int k = 0; k < num_sing_vals_; k++)
// 	    {
// 	      fout << (*S_)(k) << std::setw(20);
// 	    }
// 	  fout.close();
	     
// 	}
//       comm_->barrier();
//       RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);        

//       if(comm_->getRank() == 0)
// 	{
// 	  std::cout << " " << std::endl;
// 	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
// 	  std::cout << " " << std::endl;
// 	}

//     }

    
//   };
  
// }

// #endif



#ifndef HDSA_RANDOMIZED_GSVD_HPP
#define HDSA_RANDOMIZED_GSVD_HPP

// This class executes the randomized GSVD solver

namespace HDSA
{

  template <class RealT>
  class Randomized_GSVD
  {
  private:
    bool reduced_space_sen_;
    int num_sing_vals_, p_, kpp_, q_, nonzero_z_dim_, theta_dim_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    HDSA::Ptr<HDSA::Vector<RealT> > S_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V_;

  public:
  
    Randomized_GSVD(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
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
      OP_Objects_subcomm->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects_subcomm->z->Enforce_Zeros();
	}
      theta_dim_ = OP_Objects_subcomm->theta->dimension();
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_subcomm = weight_matrices_factory_->Construct_Weight_Matrices(theta_,subcomm); 
      HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_subcomm = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects_subcomm);
      nonzero_z_dim_ = Nom_subcomm->Get_nonzero_z_dim();
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

      comm_->barrier();
      std::clock_t timer_Loop_1 = std::clock();
      for(int k = 0; k < kpp_; k++)
      	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
      	    {
	      std::clock_t timer_range_iter_k = std::clock();
	      HDSA::Ptr<HDSA::Vector<RealT> >  theta_vec_random = Sen_Op_subcomm->Generate_Random_theta_Vector();
      	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
      	      Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm,theta_vec_random);
	      
	      Y->Write_Vector_to_Column(k,z_vec_subcomm);

	      RealT Time_range_iter_k = static_cast<RealT>(std::clock()-timer_range_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
      	      if(subcomm->getRank() == 0)
      		{
      		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the range sketching loop for local sensitivity number " << sample_index_ << " in " << Time_range_iter_k << " seconds" << std::endl;
      		}
      	    }
      	}
  
      comm_->barrier();
      proc_dist->Broadcast_Matrix(Y);
      RealT Time_Loop_1 = static_cast<RealT>(std::clock()-timer_Loop_1)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      // Need to orthogonalize the sketch with respect to the z weight matrix
      std::clock_t timer_ortho_z = std::clock();
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > SQ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_); 
      std::string type = "z";
      weight_matrices_subcomm->CholQR(Q,Y,type,OP_Objects_subcomm,SQ);
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
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_); 	  

      	  for(int k = 0; k < kpp_; k++)
      	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
      		{
		  std::clock_t timer_subspace_iter_k = std::clock();
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  
		  SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);
		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm, z_vec_subcomm_1);		  

		  Y_subspace_iter_param->Write_Vector_to_Column(k,theta_vec_subcomm);
		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the first loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
		    }
		}
      	    }
	  comm_->barrier();
	  proc_dist->Broadcast_Matrix(Y_subspace_iter_param);
	  
	  std::clock_t timer_ortho = std::clock();
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > W_Q_subspace_iter_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
	  type = "theta inverse";
	  weight_matrices_subcomm->CholQR(Q_subspace_iter_param,Y_subspace_iter_param,type,OP_Objects_subcomm,W_Q_subspace_iter_param);
      	  comm_->barrier();
	  RealT Time_ortho = static_cast<RealT>(std::clock()-timer_ortho)/static_cast<RealT>(CLOCKS_PER_SEC);
	  if(comm_->getRank() == 0)
	    {
	      std::cout << "Completed parameter orthogonalization for local sensitivity number " << sample_index_<< " in " << Time_ortho << " seconds" << std::endl;
	    }
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_subspace_iter_opt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
      	  for(int k = 0; k < kpp_; k++)
      	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
      		{
		  std::clock_t timer_subspace_iter_k = std::clock();
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
      		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();
 
		  W_Q_subspace_iter_param->Write_Column_to_Vector(k,theta_vec_subcomm);
      		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm, theta_vec_subcomm);
		  Y_subspace_iter_opt->Write_Vector_to_Column(k,z_vec_subcomm);

		  RealT Time_subspace_iter_k = static_cast<RealT>(std::clock()-timer_subspace_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the second loop of subspace iteration " << j+1  << " for local sensitivity number " << sample_index_<< " in " << Time_subspace_iter_k << " seconds" << std::endl;
		    }
      		}
      	    }
      	  comm_->barrier();
	  proc_dist->Broadcast_Matrix(Y_subspace_iter_opt);   

      	  type = "z";
	  weight_matrices_subcomm->CholQR(Q,Y_subspace_iter_opt,type,OP_Objects_subcomm,SQ);
      	  comm_->barrier();
	  
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
      
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);     
 
      comm_->barrier();
      std::clock_t timer_Loop_2 = std::clock();
      for(int k = 0; k < kpp_; k++)
      	{
	  if(proc_dist->Does_Processor_Own_Vector(k))
      	    { 
	      std::clock_t timer_proj_iter_k = std::clock();
	      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm_1 = OP_Objects_subcomm->theta->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm_2 = OP_Objects_subcomm->theta->Clone();
	      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();

	      SQ->Write_Column_to_Vector(k,z_vec_subcomm_1);
	      Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm_1, z_vec_subcomm_1);
      	      weight_matrices_subcomm->Apply_theta_Weight_Mat_Inverse(theta_vec_subcomm_2,theta_vec_subcomm_1);
	      
	      Q_proj->Write_Vector_to_Column(k,theta_vec_subcomm_2);
	      RealT Time_proj_iter_k = static_cast<RealT>(std::clock()-timer_proj_iter_k)/static_cast<RealT>(CLOCKS_PER_SEC);
      	      if(subcomm->getRank() == 0)
      		{
      		  std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the projection loop for local sensitivity number " << sample_index_ << " in " << Time_proj_iter_k << " seconds" << std::endl;
      		}
      	    }
      	}
      
      comm_->barrier();
      proc_dist->Broadcast_Matrix(Q_proj);
      RealT Time_Loop_2 = static_cast<RealT>(std::clock()-timer_Loop_2)/static_cast<RealT>(CLOCKS_PER_SEC);
  
      // Need to orthogonalize the projection with respect to the parameter mass matrix
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_param = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);  
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_Trans = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);  
      type = "theta";
      weight_matrices_subcomm->CholQR(Q_param,Q_proj,type,OP_Objects_subcomm, HDSA::nullPtr, R_Trans);

      // Need to compute the SVD of R^T
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      S_ = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
      HDSA::Linear_Algebra::SVD<RealT>(R_Trans, U, VT, S_); 
      U_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(nonzero_z_dim_,kpp_);
      Q->Multiply(U_,U);
      V_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(theta_dim_,kpp_);
      Q_param->Multiply(V_,VT,false,true);

      // The columns of U contain the approximate opt singular vectors
      // The columns of V contain the approximate parameter singular vectors
        
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

      if( OP_Objects_subcomm->theta->dimension() == theta_->dimension() )
	{
	  Write_Solution(OP_Objects_subcomm->z,weight_matrices_subcomm);
	}
      else
	{
	  Write_Solution_Model_Error(OP_Objects_subcomm->z);
	}
    }

    void Write_Solution(HDSA::Ptr<HDSA::Vector<RealT> > & z, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices)
    {
      std::clock_t timer_SI = std::clock();
      // Estimate sensitivity indices
      std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_vecs;
      theta_vecs.resize(num_sing_vals_);
      HDSA::Ptr<HDSA::Vector<RealT> > theta = theta_->Clone();
      for(int k = 0; k < num_sing_vals_; k++)
	{
	  theta->zero();
	  for(int i = 0; i < theta_dim_; i++)
	    {
	      theta->Replace_Element(i,(*V_)(i,k));
	    }
	  theta_vecs[k] = theta_->Clone();
	  weight_matrices->Apply_theta_Weight_Mat(theta_vecs[k],theta);
	}
      
      HDSA::Ptr<HDSA::Vector<RealT> > theta_matvec = theta_->Clone();
      std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim_,0.0);

      for(int i = 0; i < theta_dim_; i++)
	{
	  theta->basis(i);
	  weight_matrices->Apply_theta_Weight_Mat(theta_matvec,theta);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      sensitivity_indices[i] += std::pow((*S_)(k),2)*std::pow((*theta_vecs[k])(i),2);
	    }
	  sensitivity_indices[i] = std::sqrt(sensitivity_indices[i])/std::sqrt((*theta_matvec)(i));
	}
      RealT Time_SI = static_cast<RealT>(std::clock()-timer_SI)/static_cast<RealT>(CLOCKS_PER_SEC);      

      std::clock_t timer_write = std::clock();
      // Write to text files
      std::string name;
      std::ofstream fout;
      name = "Sensitivity_Indices_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < theta_dim_; k++)
	{
	  fout << sensitivity_indices[k] << std::setw(20);
	}
      fout.close();

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

      name = "theta_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int i  = 0; i < theta_dim_; i++)
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {    
	      fout << (*V_)(i,k) << std::setw(20);
	    }
	  fout << " " << std::endl;
	}
      fout.close();
	    
      name = "Singular_Values_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < num_sing_vals_; k++)
	{
	  fout << (*S_)(k) << std::setw(20);
	}
      fout.close();

      name = "theta_Sample_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < theta_dim_; k++)
	{
	  fout << (*theta_)(k) << std::setw(20);
	}
      fout.close();
    
      RealT Time_write = static_cast<RealT>(std::clock()-timer_write)/static_cast<RealT>(CLOCKS_PER_SEC);     

      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Time estimating sensitivity indices from GSVD: " << Time_SI << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << "Time writing to file: " << Time_write << " seconds for local sensitivity number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}

    }

    void Write_Solution_Model_Error(HDSA::Ptr<HDSA::Vector<RealT> > & z)
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
	  
	  name = "theta_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
	  fout.open(name);
	  for(int i  = 0; i < theta_dim_; i++)
	    {
	      for(int k = 0; k < num_sing_vals_; k++)
		{    
		  fout << (*V_)(i,k) << std::setw(20);
		}
	      fout << " " << std::endl;
	    }
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

    
  };
  
}

#endif
