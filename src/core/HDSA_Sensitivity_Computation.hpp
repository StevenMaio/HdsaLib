#ifndef HDSA_SENSITIVITY_COMPUTATION_HPP
#define HDSA_SENSITIVITY_COMPUTATION_HPP

// This class performs the sensitivity computation at a fixed parameter sample

namespace HDSA
{

  template <class RealT>
  class Sensitivity_Computation{
  private:
   
    HDSA::Ptr<HDSA::Vector<RealT> > theta_;
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    const HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    std::vector<int> Comm_Split_Ranks_;

    std::string gsvd_formulation_;
    std::string gsvd_solver_;
    bool compute_gsvd_sensitivities_;
    bool compute_direct_sensitivities_;
    bool compute_model_error_sensitivities_;
    bool reduced_space_sen_;
    bool randomized_hessian_eigenvector_projection_;
    bool randomized_LIS_eigenvector_projection_;
    bool randomized_LIS_construct_B_;
    bool FD_check_;
    bool construct_weight_mat_, construct_K_, construct_B_, construct_B_transpose_, construct_misfit_hessian_, construct_regularization_hessian_, construct_model_error_objects_, check_solution_read_;
    int sample_index_;
    
  public:
    
    Sensitivity_Computation(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
			    const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
			    const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, std::vector<int> & Comm_Split_Ranks, const int & myRank_position):
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory), Comm_Split_Ranks_(Comm_Split_Ranks)
    {
      FD_check_ = parlist_sensitivity_->sublist("Formulation").get("Finite Difference Check", false);
      construct_weight_mat_  = parlist_sensitivity_->sublist("Formulation").get("Construct Weight Matrices",false);
      construct_K_  = parlist_sensitivity_->sublist("Formulation").get("Construct K",false);
      construct_B_  = parlist_sensitivity_->sublist("Formulation").get("Construct B",false);
      construct_B_transpose_  = parlist_sensitivity_->sublist("Formulation").get("Construct B Transpose",false);
      construct_misfit_hessian_  = parlist_sensitivity_->sublist("Formulation").get("Construct Misfit Hessian",false);
      construct_regularization_hessian_  = parlist_sensitivity_->sublist("Formulation").get("Construct Regularization Hessian",false);
      construct_model_error_objects_  = parlist_sensitivity_->sublist("Formulation").get("Construct Model Error Objects",false);
      check_solution_read_ = parlist_sensitivity_->sublist("Formulation").get("Check Solution",false);

      reduced_space_sen_ = parlist_sensitivity_->sublist("Formulation").get("Reduced Space Sensitivities",true);
      compute_gsvd_sensitivities_ = parlist_sensitivity_->sublist("Formulation").get("Compute GSVD Sensitivities", true);
      compute_direct_sensitivities_ = parlist_sensitivity_->sublist("Formulation").get("Compute Direct Sensitivities", false);
      compute_model_error_sensitivities_ = parlist_sensitivity_->sublist("Formulation").get("Compute Model Error Sensitivities", false);
      randomized_hessian_eigenvector_projection_ = parlist_sensitivity_->sublist("Formulation").get("Randomized Hessian Eigenvector Projection", false);
      randomized_LIS_eigenvector_projection_ = parlist_sensitivity_->sublist("Formulation").get("Randomized Likelihood Informed Subspace Eigenvector Projection", false);
      randomized_LIS_construct_B_ = parlist_sensitivity_->sublist("Randomized Likelihood Informed Subspace EVP").get("Construct B", false);

      gsvd_formulation_ = parlist_sensitivity_->sublist("GSVD Solver").get("GSVD Formulation", "GSVD");
      gsvd_solver_ = parlist_sensitivity_->sublist("GSVD Solver").get("GSVD Solver", "Randomized");

      sample_index_ = myRank_position + 1;

      if(FD_check_ || construct_weight_mat_ || construct_K_ || construct_B_ || construct_B_transpose_ || construct_misfit_hessian_ || construct_regularization_hessian_ ||
	 construct_model_error_objects_ || check_solution_read_)
	{
	  Debugging_Tools();
	}
      else
	{
	  Solver_Call();
	}
    }

    void Debugging_Tools(void)
    {
      if(construct_model_error_objects_)
	{
	  HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  Opt_Problem_Objects_Model_Error<RealT> &eOP_Objects = dynamic_cast<Opt_Problem_Objects_Model_Error<RealT>&>(*OP_Objects);
	  eOP_Objects.Construct_Model_Error_Objects_Test();
	}

      if(construct_weight_mat_)
	{
	  HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices_ = weight_matrices_factory_->Construct_Weight_Matrices(theta_,comm_);
	  weight_matrices_->Construct_Weight_Matrix_Test(OP_Objects->theta,OP_Objects->z);
	}

      if(construct_K_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
	  HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op;
	  if(reduced_space_sen_)
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  else
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  Sen_Op->Construct_K_Test();
	}

      if(construct_misfit_hessian_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
	  HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op;
	  if(reduced_space_sen_)
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  else
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  Sen_Op->Construct_Misfit_Hessian_Test();
	}

      if(construct_regularization_hessian_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
	  HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op;
	  if(reduced_space_sen_)
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  else
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  Sen_Op->Construct_Regularization_Hessian_Test();
	}
      
      if(construct_B_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
	  HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op;
	  if(reduced_space_sen_)
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  else
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  Sen_Op->Construct_B_Test();
	}

      if(construct_B_transpose_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom = HDSA::makePtr<HDSA::Nominal_Data<RealT> >(parlist_sensitivity_,OP_Objects);
	  HDSA::Ptr<HDSA::Sensitivity_Operators<RealT> > Sen_Op;
	  if(reduced_space_sen_)
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_RS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  else
	    {
	      Sen_Op = HDSA::makePtr<HDSA::Sensitivity_Operators_FS<RealT> >(OP_Objects, Nom, comm_, Comm_Split_Ranks_);
	    }
	  Sen_Op->Construct_B_Transpose_Test();
	}

      if(FD_check_)
	{
	  HDSA::Ptr<HDSA::Finite_Difference_Sensitivity_Indices<RealT> > solver = HDSA::makePtr<HDSA::Finite_Difference_Sensitivity_Indices<RealT> >(theta_, parlist_sensitivity_, 
																		     comm_, OP_Objects_Factory_, weight_matrices_factory_);
	  solver->Compute();
	}

      if(check_solution_read_)
	{
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
	  OP_Objects->Load_Optimal_Solution();
	  bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
	  if(enforce_z_zeros_)
	    {
	      OP_Objects->z->Enforce_Zeros();
	    }
	  if(reduced_space_sen_)
	    {
	      RealT obj_val = OP_Objects->rs_obj->value(*OP_Objects->z,*OP_Objects->theta,true);
	      HDSA::Ptr<HDSA::Vector<RealT> > grad = OP_Objects->z->Clone();
	      OP_Objects->rs_obj->gradient_z(*grad,*OP_Objects->z,*OP_Objects->theta,true);
	      RealT grad_norm = grad->norm();
	      if(comm_->getRank() == Comm_Split_Ranks_[0])
		{
		  std::cout << "Objective function = " << obj_val << " and gradient norm = " << grad_norm << std::endl;
		}
	      // Write to text files
	      std::string name;
	      std::ofstream fout;
	      name = "gradient.txt";
	      fout.open(name);
	      for(int k = 0; k < grad->dimension(); k++)
		{
		  fout << (*grad)(k) << std::setw(20);
		}
	      fout.close();
	    }
	}

    }

    void Solver_Call(void)
    {     
      // Randomized hessian eigenvalue problem
      if(randomized_hessian_eigenvector_projection_)
	{
	  HDSA::Ptr<HDSA::Randomized_Hessian_EVP<RealT> > solver = HDSA::makePtr<HDSA::Randomized_Hessian_EVP<RealT> >(theta_, parlist_sensitivity_, comm_, 
														       OP_Objects_Factory_, weight_matrices_factory_,sample_index_);
	  solver->Compute();
	}
      
      // Randomized likelihood informed subspace eigenvalue problem
      if(randomized_LIS_eigenvector_projection_)
	{
	  HDSA::Ptr<HDSA::Randomized_Likelihood_Informed_Subspace_EVP<RealT> > solver = HDSA::makePtr<HDSA::Randomized_Likelihood_Informed_Subspace_EVP<RealT> >(theta_, parlist_sensitivity_, comm_, 
																				 OP_Objects_Factory_, weight_matrices_factory_,
																				 sample_index_);
	  if(randomized_LIS_construct_B_)
	    {
	      solver->Construct_B_Parallel();
	    }
	  else
	    {
	      solver->Compute();
	    }
	}

      // Direct sensitivities
      if(compute_direct_sensitivities_)
	{
	  HDSA::Ptr<HDSA::Direct_Sensitivity_Computation<RealT> > solver = HDSA::makePtr<HDSA::Direct_Sensitivity_Computation<RealT> >(theta_, parlist_sensitivity_, comm_, 
																       OP_Objects_Factory_, weight_matrices_factory_, sample_index_);
	  solver->Compute();
	}
      
      // sensitivity index estimatation via GSVD
      if(compute_gsvd_sensitivities_)
	{
	  // GSVD directly
	  if(gsvd_formulation_ == "GSVD")
	    {
	      if(gsvd_solver_ != "Randomized")
		{
		  std::cout << "Error: GSVD Formulation = GSVD requires GSVD Solver = Randomized" << std::endl;
		}
	      else
		{	    
		  HDSA::Ptr<HDSA::Randomized_GSVD<RealT> > solver = HDSA::makePtr<HDSA::Randomized_GSVD<RealT> >(theta_, parlist_sensitivity_, comm_, OP_Objects_Factory_,
														 weight_matrices_factory_, sample_index_);
		  solver->Compute();
		}
	    }
	  else
	    {
	      if(gsvd_solver_ != "Randomized")
		{
		  HDSA::Ptr<HDSA::Iterative_GEVP<RealT> > solver = HDSA::makePtr<HDSA::Iterative_GEVP<RealT> >(theta_, parlist_sensitivity_, comm_, OP_Objects_Factory_, weight_matrices_factory_, 
													       Comm_Split_Ranks_,sample_index_);
		  solver->Compute();
		}
	      else
		{	    
		  // GSVD via randomized generalized eigenvalue problem
		  HDSA::Ptr<HDSA::Randomized_GEVP<RealT> > solver = HDSA::makePtr<HDSA::Randomized_GEVP<RealT> >(theta_, parlist_sensitivity_, comm_, OP_Objects_Factory_, 
														 weight_matrices_factory_, sample_index_);
		  solver->Compute();
		}  
	    }
	}

      // sensitivity index with respect to model error
      if(compute_model_error_sensitivities_)
	{	    
	  HDSA::Ptr<HDSA::Randomized_GSVD_Model_Error<RealT> > solver = HDSA::makePtr<HDSA::Randomized_GSVD_Model_Error<RealT> >(theta_, parlist_sensitivity_, comm_, OP_Objects_Factory_,
																 weight_matrices_factory_, sample_index_);
	  solver->Compute();
	}

    }
    
  };

}

#endif
