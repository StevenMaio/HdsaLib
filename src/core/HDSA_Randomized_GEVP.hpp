#ifndef HDSA_RANDOMIZED_GEPV_HPP
#define HDSA_RANDOMIZED_GEPV_HPP

namespace HDSA
{

  template <class RealT>
  class Randomized_GEVP{
  private:
    bool reduced_space_sen_;
    std::string gevp_formulation_;
    int num_sing_vals_, p_, kpp_, q_, nonzero_z_dim_, theta_dim_, n_nonzero_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int sample_index_;
    HDSA::Ptr<HDSA::Vector<RealT> > W_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_;
    
  public:
    
    Randomized_GEVP(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
		    const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
		    const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, const int & sample_index): 
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), sample_index_(sample_index)
    {
      reduced_space_sen_ = parlist_sensitivity->sublist("Formulation").get("Reduced Space Sensitivities",true);
      gevp_formulation_ = parlist_sensitivity->sublist("GSVD Solver").get("GSVD Formulation","Jordan-Wielandt");
      num_sing_vals_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Singular Values",1); 
      p_ = parlist_sensitivity->sublist("GSVD Solver").get("Oversampling Factor",20);
      q_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Subspace Iterations",0);
      if(gevp_formulation_ == "Gram")
	{
	  kpp_ = num_sing_vals_ + p_;
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  kpp_ = 2*num_sing_vals_ + p_;
	}
      theta_dim_ = theta->dimension();
    }

    void Compute(void)
    {
      // This function executes a randomized generalized eigenvalue solver to estimate the GSVD
      
      std::clock_t timer_eigen = std::clock();
      
      HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist = HDSA::makePtr<HDSA::Processor_Distribution<RealT> >(comm_,kpp_);
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
      n_nonzero_ = nonzero_z_dim_ + theta_dim_;

      HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op_subcomm;
      if(reduced_space_sen_)
      	{
      	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
      	}
      else
      	{
	  Sen_Op_subcomm = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects_subcomm,Nom_subcomm,subcomm, proc_dist->Get_Comm_Split_Ranks());
      	}

      int dim = 0;
      if(gevp_formulation_ == "Gram")
	{
	  dim = theta_dim_;
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  dim = n_nonzero_;
	}     

      // Sketch the range of D
      if(comm_->getRank() == 0)
	{
	  std::cout << "Beginning the first loop of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}

      comm_->barrier();
      std::clock_t timer_Loop_1 = std::clock();
      // Instantiate matrix to interface with LA routines     
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,kpp_); 

      if(gevp_formulation_ == "Gram")
	{
	  comm_->barrier();    
	  for(int k = 0; k < kpp_; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
		{
		  HDSA::Ptr<HDSA::Vector<RealT> >  theta_vec_random = Sen_Op_subcomm->Generate_Random_theta_Vector();
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();

		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm_1,theta_vec_random);
		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_2,z_vec_subcomm_1);
		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm,z_vec_subcomm_2);
		  weight_matrices_subcomm->Apply_theta_Weight_Mat_Inverse(theta_vec_random,theta_vec_subcomm);  

		  Y->Write_Vector_to_Column(k,theta_vec_random);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the first set of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
		    }
		}
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  for(int k = 0; k < kpp_; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
	      {     
		  HDSA::Ptr<HDSA::Vector<RealT> >  theta_vec_random = Sen_Op_subcomm->Generate_Random_theta_Vector();
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_random = Sen_Op_subcomm->Generate_Random_z_Vector();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm = OP_Objects_subcomm->z->Clone();

		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm,theta_vec_random);

        	  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      Y->Replace_Element(i,k,(*z_vec_subcomm)(z_vec_subcomm->Get_map_reduced_to_full(i)));
		    }
	  
		  z_vec_subcomm->zero();
		  theta_vec_random->zero();

		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm,z_vec_random);
		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_random,z_vec_subcomm);
		  weight_matrices_subcomm->Apply_theta_Weight_Mat_Inverse(theta_vec_subcomm, theta_vec_random);
		  for(int i = nonzero_z_dim_; i < n_nonzero_; i++)
		    {
		      Y->Replace_Element(i,k,(*theta_vec_subcomm)(i-nonzero_z_dim_));
		    }

		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the first set of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
		    }
		}
	    }
	}
      
      comm_->barrier();
      proc_dist->Broadcast_Matrix(Y);
      RealT Time_Loop_1 = static_cast<RealT>(std::clock()-timer_Loop_1)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      // Orthogonalize
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,kpp_);
      std::string type;
      if(gevp_formulation_ == "Gram")
	{
	  type = "theta";
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  type = "joint";
	}
      weight_matrices_subcomm->CholQR(Q,Y,type,OP_Objects_subcomm);
      comm_->barrier();
      
      // Project D onto the sketched subspace
      if(comm_->getRank() == 0)
	{
	  std::cout << "Beginning the second loop of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}
      comm_->barrier();
      std::clock_t timer_Loop_2 = std::clock();
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_proj = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,kpp_);      

      if(gevp_formulation_ == "Gram")
	{
	  comm_->barrier();
	  for(int k = 0; k < kpp_; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))
	      {	  
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();

		  Q->Write_Column_to_Vector(k,theta_vec_subcomm);

		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm_1,theta_vec_subcomm);
		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_2,z_vec_subcomm_1);
		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm,z_vec_subcomm_2);
		  
		  Q_proj->Write_Vector_to_Column(k,theta_vec_subcomm);
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the second set of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
		    }
		}
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  for(int k = 0; k < kpp_; k++)
	    {
	      if(proc_dist->Does_Processor_Own_Vector(k))	
	      { 
		  HDSA::Ptr<HDSA::Vector<RealT> > theta_vec_subcomm = OP_Objects_subcomm->theta->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_1 = OP_Objects_subcomm->z->Clone();
		  HDSA::Ptr<HDSA::Vector<RealT> > z_vec_subcomm_2 = OP_Objects_subcomm->z->Clone();

		  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      z_vec_subcomm_1->Replace_Element(z_vec_subcomm_1->Get_map_reduced_to_full(i),(*Q)(i,k));
		    }

		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_2,z_vec_subcomm_1);
		  Sen_Op_subcomm->Apply_Sensitivity_Operator_Transpose(theta_vec_subcomm,z_vec_subcomm_2);
		  for(int i = nonzero_z_dim_; i < n_nonzero_; i++)
		    {
		      Q_proj->Replace_Element(i,k,(*theta_vec_subcomm)(i-nonzero_z_dim_));
		    }

		  theta_vec_subcomm->zero();
		  for(int i = 0; i < theta_dim_; i++)
		    {
		      theta_vec_subcomm->Replace_Element(i,(*Q)(i+nonzero_z_dim_,k));
		    }
		  z_vec_subcomm_1->zero();
		  z_vec_subcomm_2->zero();
		  Sen_Op_subcomm->Apply_Sensitivity_Operator(z_vec_subcomm_2,theta_vec_subcomm);
		  weight_matrices_subcomm->Apply_z_Weight_Mat(z_vec_subcomm_1,z_vec_subcomm_2);
		  
		  for(int i = 0; i < nonzero_z_dim_; i++)
		    {
		      Q_proj->Replace_Element(i,k,(*z_vec_subcomm_1)(z_vec_subcomm_1->Get_map_reduced_to_full(i)));
		    }
       
		  if(subcomm->getRank() == 0)
		    {
		      std::cout << "Completed iteration " << k+1 << " out of " << kpp_ << " in the second set of operator applies with local sensitivity sample number " << sample_index_ << std::endl;
		    }
		}
	    }
	}
      
      comm_->barrier();
      proc_dist->Broadcast_Matrix(Q_proj);
      RealT Time_Loop_2 = static_cast<RealT>(std::clock()-timer_Loop_2)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > T = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      RealT val = 0.0;
      
      for(int j = 0; j < kpp_; j++)
	{
	  for(int i = 0; i < kpp_; i++)
	    {
	      if(j <= i)
		{
		  
		  for(int k = 0; k < dim; k++)
		    {
		      val += ((*Q)(k,i))*(*Q_proj)(k,j);
		    }
		  T->Replace_Element(i,j,val);
		  val = 0.0;
		} else
		{
		  T->Replace_Element(i,j,(*T)(j,i));
		}
	    }
	}
      
      // Compute eigenvalue decomposition for the small matrix T
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp_,kpp_);
      W_ = HDSA::makePtr<Std_Vector<RealT> >(kpp_);
      HDSA::Linear_Algebra::Symmetric_Eig_Decomposition(T,V,W_);
      U_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,kpp_);
      Q->Multiply(U_,V);

      // Take square root of eigenvalues for Gram formulation so that W_ contains the singular values
      if(gevp_formulation_ == "Gram")
	{
	  for(int k = 0; k < kpp_; k++)
	    {
	      W_->Replace_Element(k,std::sqrt((*W_)(k)));
	    }
	}

      RealT Time_eigen = static_cast<RealT>(std::clock()-timer_eigen)/static_cast<RealT>(CLOCKS_PER_SEC);
      
      if(comm_->getRank() == 0)
	{
	  std::cout << " " << std::endl;
	  std::cout << "Total time for eigenvalue computation: " << Time_eigen << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by the first loop: " << Time_Loop_1 << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by the second loop: " << Time_Loop_2 << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << "Time contributed by everything else: " << Time_eigen-Time_Loop_1-Time_Loop_2 << " seconds with local sensitivity sample number " << sample_index_ << std::endl;
	  std::cout << " " << std::endl;
	}

      Write_Solution(Nom_subcomm, weight_matrices_subcomm);  

    }

    void Write_Solution(HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices)
    {
      std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim_,0.0);
      std::clock_t timer_SI = std::clock();
      // Estimate sensitivity indices
      std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_sing_vecs;
      theta_sing_vecs.resize(num_sing_vals_);

      if(gevp_formulation_ == "Gram")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      theta_sing_vecs[k] = theta_->Clone();
	      for(int i = 0; i < theta_dim_; i++)
		{
		  theta_sing_vecs[k]->Replace_Element(i,(*U_)(i,k));
		}
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      theta_sing_vecs[k] = theta_->Clone();
	      for(int i = 0; i < theta_dim_; i++)
		{
		  theta_sing_vecs[k]->Replace_Element(i,(*U_)(nonzero_z_dim_+i,k));
		}
	    }
	}
      // Normalization back to proper inner product space
      for(int k = 0; k < num_sing_vals_; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > theta_matvec = theta_->Clone();
	  weight_matrices->Apply_theta_Weight_Mat(theta_matvec,theta_sing_vecs[k]);
	  RealT normalize = std::sqrt(theta_matvec->dot(*theta_sing_vecs[k]));
	  theta_sing_vecs[k]->scale(1.0/normalize);
	}

      HDSA::Ptr<HDSA::Vector<RealT> > theta_i = theta_->Clone();	
      HDSA::Ptr<HDSA::Vector<RealT> > theta_matvec = theta_->Clone();
      for(int i = 0; i < theta_dim_; i++)
	{
	  theta_matvec->zero();
	  theta_i->basis(i);
	  weight_matrices->Apply_theta_Weight_Mat(theta_matvec,theta_i);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      sensitivity_indices[i] += std::pow((*W_)(k),2)*std::pow(theta_sing_vecs[k]->dot(*theta_matvec),2);
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

      name = "theta_Singular_Vector_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int i  = 0; i < theta_dim_; i++)
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {    
	      fout << (*theta_sing_vecs[k])(i) << std::setw(20);
	    }
	  fout << " " << std::endl;
	}
      fout.close();
	    
      name = "Singular_Values_" + std::to_string(sample_index_) + ".txt";
      fout.open(name);
      for(int k = 0; k < num_sing_vals_; k++)
	{
	  fout << (*W_)(k) << std::setw(20);
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

  };

}

#endif
