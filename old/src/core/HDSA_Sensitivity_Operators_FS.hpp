#ifndef HDSA_SENSITIVITY_OPERATORS_FS_HPP
#define HDSA_SENSITIVITY_OPERATORS_FS_HPP

// This class contains the various operators needed throughout the computation

namespace HDSA
{
  
  template <class RealT>
  class Sensitivity_Operators_FS: public HDSA::Sensitivity_Operators<RealT>{
  protected:

    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_;
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_KKT_;
    HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_;
    Ptr<const HDSA::Comm<int> > comm_;
    std::vector<int> Comm_Split_Ranks_;
    RealT h_;
    HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal_u_, grad_nominal_z_, grad_nominal_lambda_;
    int z_dim_, theta_dim_;
    bool verbosity_;
    bool conserve_memory_;
    
  public:
   
    Sensitivity_Operators_FS(const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, const HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, 
			     const HDSA::Ptr<const HDSA::Comm<int> > & comm, const std::vector<int> & Comm_Split_Ranks): Sensitivity_Operators<RealT>(OP_Objects, Nom, Comm_Split_Ranks) 						 
    {
      OP_Objects_ = OP_Objects;
      Nom_ = Nom;
      comm_ = comm;
      Comm_Split_Ranks_ = Comm_Split_Ranks;
      z_dim_ = OP_Objects_->z->dimension();
      theta_dim_ = OP_Objects_->theta->dimension();
      
      // Set verbosity = false for basic output, verbosity = true for output from the status of each operator apply
      verbosity_ = Nom_->Get_parlist_sensitivity()->sublist("Formulation").get("Verbosity",false);

     // If true, construct an Opt_Problem_Objects which is devoted to the KKT solve which avoids reassembly when the parameters are perturbed for evaluations of B and B^T
      conserve_memory_ = Nom_->Get_parlist_sensitivity()->sublist("Formulation").get("Conserve Memory",false);     
      if(!conserve_memory_)
	{
	  OP_Objects_KKT_ = OP_Objects->Construct_Opt_Problem_Objects(OP_Objects_->theta,comm_);
	  OP_Objects_KKT_->u->set(*OP_Objects_->u);
	  OP_Objects_KKT_->z->set(*OP_Objects_->z);
	  OP_Objects_KKT_->lambda->set(*OP_Objects_->lambda);
	}

      h_ = Nom_->Get_parlist_sensitivity()->sublist("Parameter Derivative").get("Finite Difference Step",1.e-4); // Read in finite difference step size from Sensitivity_input.xml
     
      grad_nominal_u_ = OP_Objects_->u->Clone();
      grad_nominal_z_ = OP_Objects_->z->Clone();
      grad_nominal_lambda_ = OP_Objects_->lambda->Clone();
      Grad_at_Nominal_Solution(); // Evaluate the gradient at the nominal solution          
    }
    
    ~Sensitivity_Operators_FS()
    { }
    
    // Compute gradient of Lagrangian at nominal solution
    void Grad_at_Nominal_Solution(void)
    {
     if(conserve_memory_)
	{
	  // Compute the state gradient
	  HDSA::Ptr<HDSA::Vector<RealT> > obj_grad_u = OP_Objects_->u->Clone();
	  OP_Objects_->fs_obj->gradient_u(*obj_grad_u,*OP_Objects_->u,*OP_Objects_->z, *OP_Objects_->theta);
	  OP_Objects_->con->jacobian_u_adjoint(*grad_nominal_u_, *OP_Objects_->lambda,*OP_Objects_->u,*OP_Objects_->z,*OP_Objects_->theta);
	  grad_nominal_u_->plus(*obj_grad_u);
	  
	  // Compute the control gradient
	  HDSA::Ptr<HDSA::Vector<RealT> > obj_grad_z = OP_Objects_->z->Clone();
	  OP_Objects_->fs_obj->gradient_z(*obj_grad_z,*OP_Objects_->u,*OP_Objects_->z, *OP_Objects_->theta);
	  OP_Objects_->con->jacobian_z_adjoint(*grad_nominal_z_, *OP_Objects_->lambda,*OP_Objects_->u,*OP_Objects_->z,*OP_Objects_->theta);
	  grad_nominal_z_->plus(*obj_grad_z);
	  
	  // Compute the adjoint gradient
	  OP_Objects_->con->value(*grad_nominal_lambda_,*OP_Objects_->u,*OP_Objects_->z,*OP_Objects_->theta);
	}
      else
	{
	  // Compute the state gradient
	  HDSA::Ptr<HDSA::Vector<RealT> > obj_grad_u = OP_Objects_KKT_->u->Clone();
	  OP_Objects_KKT_->fs_obj->gradient_u(*obj_grad_u,*OP_Objects_KKT_->u,*OP_Objects_KKT_->z, *OP_Objects_KKT_->theta);
	  OP_Objects_KKT_->con->jacobian_u_adjoint(*grad_nominal_u_, *OP_Objects_KKT_->lambda,*OP_Objects_KKT_->u,*OP_Objects_KKT_->z,*OP_Objects_KKT_->theta);
	  grad_nominal_u_->plus(*obj_grad_u);
	  
	  // Compute the control gradient
	  HDSA::Ptr<HDSA::Vector<RealT> > obj_grad_z = OP_Objects_KKT_->z->Clone();
	  OP_Objects_KKT_->fs_obj->gradient_z(*obj_grad_z,*OP_Objects_KKT_->u,*OP_Objects_KKT_->z, *OP_Objects_KKT_->theta);
	  OP_Objects_KKT_->con->jacobian_z_adjoint(*grad_nominal_z_, *OP_Objects_KKT_->lambda,*OP_Objects_KKT_->u,*OP_Objects_KKT_->z,*OP_Objects_KKT_->theta);
	  grad_nominal_z_->plus(*obj_grad_z);
	  
	  // Compute the adjoint gradient
	  OP_Objects_KKT_->con->value(*grad_nominal_lambda_,*OP_Objects_KKT_->u,*OP_Objects_KKT_->z,*OP_Objects_KKT_->theta);
	}
    }

    // Apply the operator B 
    void Apply_B(HDSA::Ptr<HDSA::Vector<RealT> > & delta_u, HDSA::Ptr<HDSA::Vector<RealT> > & delta_z, HDSA::Ptr<HDSA::Vector<RealT> > & delta_lambda, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_theta)
    {
      // u
      OP_Objects_->fs_obj->hessVec_u_theta(*delta_u, *delta_theta, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta);
      HDSA::Ptr<HDSA::Vector<RealT> > Ju = delta_u->Clone();
      OP_Objects_->con->hessian_u_theta_adjoint(*Ju, *OP_Objects_->lambda, *delta_theta, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta);
      delta_u->plus(*Ju);
      delta_u->scale(-1.0);

      // z 
      OP_Objects_->fs_obj->hessVec_z_theta(*delta_z, *delta_theta, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta);
      HDSA::Ptr<HDSA::Vector<RealT> > Jz = delta_z->Clone();
      OP_Objects_->con->hessian_z_theta_adjoint(*Jz, *OP_Objects_->lambda, *delta_theta, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta);
      delta_z->plus(*Jz);
      delta_z->scale(-1.0);

      // lambda
      OP_Objects_->con->jacobian_theta(*delta_lambda, *delta_theta, *OP_Objects_->u, *OP_Objects_->z, *OP_Objects_->theta, true, grad_nominal_lambda_);
      delta_lambda->scale(-1.0);
    }
    
    // Apply the operator B^T 
    void Apply_B_Transpose(HDSA::Ptr<HDSA::Vector<RealT> > & delta_theta, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_u, 
			   const HDSA::Ptr<HDSA::Vector<RealT> > & delta_z, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_lambda)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > theta_pert = delta_theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > u_pert = delta_u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_pert = delta_z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > lambda_pert = delta_lambda->Clone();
      
      for(int i = 0; i < theta_dim_; i++)
	{
	  theta_pert->basis(i);
	  Apply_B(u_pert,z_pert,lambda_pert,theta_pert);
	  RealT val = u_pert->dot(*delta_u) + z_pert->dot(*delta_z) + lambda_pert->dot(*delta_lambda);
	  delta_theta->Replace_Element(i,val);
	}
    }

    // KKT system solve via solving a quadratic linear program
    void Apply_K_Inverse(const HDSA::Ptr<HDSA::Vector<RealT> > & sol_u, const HDSA::Ptr<HDSA::Vector<RealT> > & sol_z, const HDSA::Ptr<HDSA::Vector<RealT> > & sol_lambda,
			 const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_u, const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_z, const HDSA::Ptr<HDSA::Vector<RealT> > & rhs_lambda)
    {
      int maxits = Nom_->Get_parlist_sensitivity()->sublist("KKT Solve").get("Maximum Iterations",OP_Objects_->z->dimension());  
      RealT tol = Nom_->Get_parlist_sensitivity()->sublist("KKT Solve").get("Tolerance",1.e-5); 

      HDSA::nullstream bhs;
      HDSA::Ptr<std::ostream> outStream;
      bool verbose = Nom_->Get_parlist_sensitivity()->sublist("KKT Solve").get("Verbosity",false);
      if( verbose && (Comm_Split_Ranks_[comm_->getRank()] == 0) )
     	{
     	  outStream = HDSA::makePtrFromRef(std::cout);
     	}
      else
     	{
     	  outStream = HDSA::makePtrFromRef(bhs);
     	}

      if(conserve_memory_)
	{
	  Grad_at_Nominal_Solution(); // This is necessary to set the state and adjoint solves so that we do not need updates inside the linear solver
	  HDSA::KKT_Solver::QP_Solve<RealT>(sol_u, sol_z, sol_lambda, rhs_u, rhs_z, rhs_lambda, OP_Objects_, maxits, tol, outStream);
	}
      else
	{
	  // By using OP_Objects_KKT_ we avoid the call to Grad_at_Nominal_Solution, thus saving a forward and adjoint solve for each call to Apply_K_Inverse
	  HDSA::KKT_Solver::QP_Solve<RealT>(sol_u, sol_z, sol_lambda, rhs_u, rhs_z, rhs_lambda, OP_Objects_KKT_, maxits, tol, outStream);
	}
 
    }

    void Apply_Sensitivity_Operator(HDSA::Ptr<HDSA::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Vector<RealT> > & theta)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_1 = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > lambda_vec_1 = OP_Objects_->lambda->Clone();    
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_2 = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > lambda_vec_2 = OP_Objects_->lambda->Clone();  

      if(verbosity_)
	{
	  std::clock_t timer = std::clock();
	  Apply_B(u_vec_1,z_vec_1,lambda_vec_1,theta);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_B in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	  timer = std::clock();
	  Apply_K_Inverse(u_vec_2,z_vec_2,lambda_vec_2,u_vec_1,z_vec_1,lambda_vec_1);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_K_inverse in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	  z->set(*z_vec_2);
	 }
       else
	 {
	   Apply_B(u_vec_1,z_vec_1,lambda_vec_1,theta);
	   Apply_K_Inverse(u_vec_2,z_vec_2,lambda_vec_2,u_vec_1,z_vec_1,lambda_vec_1);
	   z->set(*z_vec_2);
	 }
    }
    
    void Apply_Sensitivity_Operator_Transpose(HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::Vector<RealT> > & z)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_1 = OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_1 = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > lambda_vec_1 = OP_Objects_->lambda->Clone();    
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_2 = OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec_2 = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > lambda_vec_2 = OP_Objects_->lambda->Clone();  

      u_vec_1->zero();
      z_vec_1->set(*z);
      lambda_vec_1->zero();

      if(verbosity_)
	{
	  std::clock_t timer = std::clock();
	  Apply_K_Inverse(u_vec_2,z_vec_2,lambda_vec_2,u_vec_1,z_vec_1,lambda_vec_1);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_K_inverse in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	  timer = std::clock();
	  Apply_B_Transpose(theta,u_vec_2,z_vec_2,lambda_vec_2);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_B_Transpose in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	}
      else
	{
	  Apply_K_Inverse(u_vec_2,z_vec_2,lambda_vec_2,u_vec_1,z_vec_1,lambda_vec_1);
	  Apply_B_Transpose(theta,u_vec_2,z_vec_2,lambda_vec_2);
	}
    }
    
    void Construct_K_Test(void)
    {
      std::cout << "Construct_K_Test is not supported in full space" << std::endl;
    }

    void Construct_Misfit_Hessian_Test(void)
    {
      std::cout << "Construct_Misfit_Hessian_Test is not supported in full space" << std::endl;
    }

    void Construct_Regularization_Hessian_Test(void)
    {
      std::cout << "Construct_Regularization_Hessian_Test is not supported in full space" << std::endl;
    }

    void Construct_B_Test(void) 
    {
      std::cout << "Construct_B_Test is not supported in full space" << std::endl;
    }

    void Construct_B_Transpose_Test(void)
    {
      std::cout << "Construct_B_Transpose_Test is not supported in full space" << std::endl;
    }

  };

}

#endif
