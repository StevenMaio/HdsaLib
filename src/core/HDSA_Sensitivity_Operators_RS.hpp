#ifndef HDSA_SENSITIVITY_OPERATORS_RS_HPP
#define HDSA_SENSITIVITY_OPERATORS_RS_HPP

namespace HDSA
{

  template <class RealT>
  class Sensitivity_Operators_RS: public HDSA::Sensitivity_Operators<RealT>{
  private:

    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_;
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_KKT_;
    HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_;
    Ptr<const HDSA::Comm<int> > comm_;
    std::vector<int> Comm_Split_Ranks_;
    RealT h_;
    HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal_;
    int z_dim_, theta_dim_;
    bool verbosity_;
    bool conserve_memory_;
    bool use_preconstructed_operators_;
    std::vector<std::vector<RealT> > preconstructed_B_;
    std::vector<std::vector<RealT> > preconstructed_K_;

  public:
    
    Sensitivity_Operators_RS(const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, const HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, 
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

      // option to load preconstructed operators
      use_preconstructed_operators_ = Nom_->Get_parlist_sensitivity()->sublist("Formulation").get("Use Preconstructed Operators",false);
      if(use_preconstructed_operators_)
	{
	  Load_Preconstructed_Operators();
	}

      // If true, construct an Opt_Problem_Objects which is devoted to the KKT solve which avoids reassembly when the parameters are perturbed for evaluations of B and B^T
      conserve_memory_ = Nom_->Get_parlist_sensitivity()->sublist("Formulation").get("Conserve Memory",false);  
      if(!conserve_memory_)
	{
	  OP_Objects_KKT_ = OP_Objects->Construct_Opt_Problem_Objects(OP_Objects_->theta,comm_);
	  OP_Objects_KKT_->z->set(*OP_Objects_->z);
	}

      h_ = Nom_->Get_parlist_sensitivity()->sublist("Parameter Derivative").get("Finite Difference Step",1.e-4); // Read in finite difference step size from Sensitivity_input.xml
      grad_nominal_ = OP_Objects_->z->Clone();
      Grad_at_Nominal_Solution(); // Evaluate the gradient at the nominal solution        
    }
    
    ~Sensitivity_Operators_RS()
    { }
    
    // Load operators for testing on small problems
    void Load_Preconstructed_Operators(void)
    {
      preconstructed_B_.resize(z_dim_);
      preconstructed_K_.resize(z_dim_);
      for(int i = 0; i < z_dim_; i++)
	{
	  preconstructed_B_[i].resize(theta_dim_);
	  preconstructed_K_[i].resize(z_dim_);
	}

      // read in B
      std::ifstream inputFileB("B.txt");          
      RealT value;
      // read the elements in the file into a vector  
      // test file open   
      if (inputFileB) {   
	for(int i = 0; i < z_dim_; i++)
	  {
	    for(int j = 0; j < theta_dim_; j++)
	      {
		inputFileB >> value;
		preconstructed_B_[i][j] = value;
	      }
	  }
      }
      else
	{
	  std::cout << "Error loading B" << std::endl;
	}  
      // read in K
      std::ifstream inputFileK("K.txt");          
      // read the elements in the file into a vector  
      // test file open   
      if (inputFileK) {   
	for(int i = 0; i < z_dim_; i++)
	  {
	    for(int j = 0; j < z_dim_; j++)
	      {
		inputFileK >> value;
		preconstructed_K_[i][j] = value;
	      }
	  }
      }
      else
	{
	  std::cout << "Error loading K" << std::endl;
	}  

    }

    // Compute gradient at the nominal solution
    void Grad_at_Nominal_Solution(void)
    {
      if(conserve_memory_)
	{
	  OP_Objects_->rs_obj->gradient_z(*grad_nominal_,*OP_Objects_->z,*OP_Objects_->theta);
	}
      else
	{
	  OP_Objects_KKT_->rs_obj->gradient_z(*grad_nominal_,*OP_Objects_KKT_->z,*OP_Objects_KKT_->theta);
	}
    }    
   
    // Apply the operator B 
    void Apply_B(HDSA::Ptr<HDSA::Vector<RealT> > & delta_z, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_theta)
    {
      if(!use_preconstructed_operators_)
	{
	  OP_Objects_->rs_obj->hessVec_z_theta(*delta_z, *delta_theta, *OP_Objects_->z, *OP_Objects_->theta,true,grad_nominal_);
	  delta_z->scale(-1.0);
	}
      else
	{
	  for(int i = 0; i < z_dim_; i++)
	    {
	      RealT val = 0.0;
	      for(int j = 0; j < theta_dim_; j++)
		{
		  val += preconstructed_B_[i][j]*(*delta_theta)(j);
		}
	      delta_z->Replace_Element(i,val);
	    }
	}
    }

    // Apply the operator B^T 
    void Apply_B_Transpose(HDSA::Ptr<HDSA::Vector<RealT> > & delta_theta, const HDSA::Ptr<HDSA::Vector<RealT> > & delta_z)
    {
      if(!use_preconstructed_operators_)
	{
	  OP_Objects_->rs_obj->hessVec_theta_z(*delta_theta, *delta_z, *OP_Objects_->z, *OP_Objects_->theta,true,grad_nominal_);
	  delta_theta->scale(-1.0);
	}
      else
	{
	  for(int i = 0; i < theta_dim_; i++)
	    {
	      RealT val = 0.0;
	      for(int j = 0; j < z_dim_; j++)
		{
		  val += preconstructed_B_[j][i]*(*delta_z)(j);
		}
	      delta_theta->Replace_Element(i,val);
	    }
	}
    }
            
    // KKT system solve via CG solver in reduced space
    void Apply_K_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & x_star, const HDSA::Ptr<HDSA::Vector<RealT> > & b)
    {
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > hessian_op;
      if(!use_preconstructed_operators_)
	{
	  // If conserve_memory_=false, use OP_Objects_KKT_ to avoid reassembly after theta perturbations introduced in B and B^T
	  if(conserve_memory_)
	    {
	      Grad_at_Nominal_Solution(); // This is necessary to set the state and adjoint solves so that we do not need updates inside the linear solver
	      hessian_op = HDSA::makePtr<Hessian_Operator<RealT> >(OP_Objects_,Nom_,grad_nominal_);
	    }
	  else
	    {
	      // By using OP_Objects_KKT_ we avoid the call to Grad_at_Nominal_Solution, thus saving a forward and adjoint solve for each call to Apply_K_Inverse
	      hessian_op = HDSA::makePtr<Hessian_Operator<RealT> >(OP_Objects_KKT_,Nom_,grad_nominal_);
	    }
	}
      else
	{
	  hessian_op = HDSA::makePtr<Hessian_Operator_Preconstructed<RealT> >(preconstructed_K_);
	}
      RealT tol = Nom_->Get_parlist_sensitivity()->sublist("KKT Solve").get("Tolerance",1.e-5);
      std::string solver = Nom_->Get_parlist_sensitivity()->sublist("KKT Solve").get("Solver","CG");
      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(x_star,b,hessian_op,tol);
    }
    
    void Apply_Sensitivity_Operator(HDSA::Ptr<HDSA::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Vector<RealT> > & theta)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec = z->Clone();
      
      if(verbosity_)
	{
	  std::clock_t timer = std::clock();
	  Apply_B(z_vec, theta);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_B in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	  timer = std::clock();
	  Apply_K_Inverse(z,z_vec);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_K_inverse in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	}
      else
	{
	  Apply_B(z_vec, theta);
	  Apply_K_Inverse(z,z_vec);
	}
    }

    void Apply_Sensitivity_Operator_Transpose(HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::Vector<RealT> > & z)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_vec = z->Clone();
      
      if(verbosity_)
	{
	  std::clock_t timer = std::clock();
	  Apply_K_Inverse(z_vec,z);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_K_inverse in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	  timer = std::clock();
	  Apply_B_Transpose(theta,z_vec);
	  std::cout << "Processor = " << Comm_Split_Ranks_[comm_->getRank()] << " completed Apply_B_Transpose in " << static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC) << " seconds" << std::endl;
	}
      else
	{
	  Apply_K_Inverse(z_vec,z);
	  Apply_B_Transpose(theta,z_vec);
	}
    }

    // Construct K
    void Construct_K_Test(void)
    {
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > A = HDSA::makePtr<Hessian_Operator<RealT> >(OP_Objects_KKT_,Nom_,grad_nominal_);
      HDSA::Ptr<HDSA::Vector<RealT> > v_in = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > v_out = OP_Objects_->z->Clone();
      int d = v_in->dimension();
      std::vector<std::vector<RealT> > H;
      H.resize(d);
      for(int k = 0; k < d; k++)
	{
	  H[k].resize(d);
	}
      for(int k = 0; k < d; k++)
	{
	  std::cout << "Computing column " << k+1 << " out of " << d << "." << std::endl;
	  v_in->basis(k);
	  v_out->zero();
	  A->matvec(v_out,v_in);
	  for(int i = 0; i < d; i++)
	    {
	      H[i][k] = (*v_out)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "K.txt";
      fout.open(name);
      for(int i = 0; i < d; i++)
	{
	  for(int j = 0; j < d; j++)
	    {
	      fout << std::setprecision(16) << H[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 
    }
    
    // Construct Mistfit_Hessian
    void Construct_Misfit_Hessian_Test(void)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > v_in = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > v_out = OP_Objects_->z->Clone();
      int d = v_in->dimension();
      std::vector<std::vector<RealT> > H;
      H.resize(d);
      for(int k = 0; k < d; k++)
	{
	  H[k].resize(d);
	}
      for(int k = 0; k < d; k++)
	{
	  std::cout << "Computing column " << k+1 << " out of " << d << "." << std::endl;
	  v_in->basis(k);
	  v_out->zero();
	  OP_Objects_->rs_obj->Misfit_hessVec_z_z(*v_out,*v_in,*OP_Objects_->z,*OP_Objects_->theta, true, grad_nominal_);	  
	  for(int i = 0; i < d; i++)
	    {
	      H[i][k] = (*v_out)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "Misfit_Hessian.txt";
      fout.open(name);
      for(int i = 0; i < d; i++)
	{
	  for(int j = 0; j < d; j++)
	    {
	      fout << std::setprecision(16) << H[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 
    }
    
    // Construct Regularization_Hessian
    void Construct_Regularization_Hessian_Test(void)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > v_in = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > v_out = OP_Objects_->z->Clone();
      int d = v_in->dimension();
      std::vector<std::vector<RealT> > H;
      H.resize(d);
      for(int k = 0; k < d; k++)
	{
	  H[k].resize(d);
	}
      for(int k = 0; k < d; k++)
	{
	  std::cout << "Computing column " << k+1 << " out of " << d << "." << std::endl;
	  v_in->basis(k);
	  v_out->zero();
	  OP_Objects_->rs_obj->Regularization_hessVec_z_z(*v_out,*v_in,*OP_Objects_->z,*OP_Objects_->theta, true, grad_nominal_);	  
	  for(int i = 0; i < d; i++)
	    {
	      H[i][k] = (*v_out)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "Regularization_Hessian.txt";
      fout.open(name);
      for(int i = 0; i < d; i++)
	{
	  for(int j = 0; j < d; j++)
	    {
	      fout << std::setprecision(16) << H[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close(); 
    }

    // Construct B
    void Construct_B_Test(void)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > e = OP_Objects_->theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > Be = OP_Objects_->z->Clone();
      std::vector<std::vector<RealT> > B;
      B.resize(z_dim_);
      for(int i = 0; i < z_dim_; i++)
	{
	  B[i].resize(theta_dim_);
	}
      for(int k = 0; k < theta_dim_; k++)
	{
	  std::cout << "Computing column " << k+1 << " out of " << theta_dim_ << "." << std::endl;
	  e->basis(k);
	  Be->zero();
	  Apply_B(Be,e);
	  for(int i = 0; i < z_dim_; i++)
	    {
	      B[i][k] = (*Be)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "B.txt";
      fout.open(name);
      for(int i = 0; i < z_dim_; i++)
	{
	  for(int j = 0; j < theta_dim_; j++)
	    {
	      fout << std::setprecision(16) << B[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();
    }
    
    // Construct B Transpose
    void Construct_B_Transpose_Test(void)
    {
      HDSA::Ptr<HDSA::Vector<RealT> > e = OP_Objects_->z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > Bte = OP_Objects_->theta->Clone();
      std::vector<std::vector<RealT> > Bt;
      Bt.resize(theta_dim_);
      for(int i = 0; i < theta_dim_; i++)
	{
	  Bt[i].resize(z_dim_);
	}
      for(int k = 0; k < z_dim_; k++)
	{
	  std::cout << "Computing column " << k+1 << " out of " << z_dim_ << "." << std::endl;
	  e->basis(k);
	  Bte->zero();
	  Apply_B_Transpose(Bte,e);
	  for(int i = 0; i < theta_dim_; i++)
	    {
	      Bt[i][k] = (*Bte)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "B_Transpose.txt";
      fout.open(name);
      for(int i = 0; i < theta_dim_; i++)
	{
	  for(int j = 0; j < z_dim_; j++)
	    {
	      fout << std::setprecision(16) << Bt[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();
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

    // Overload HDSA::Linear_Operator to take matrix vector products
    template <class ScalarType>
    class Hessian_Operator_Preconstructed : public HDSA::Linear_Operator<ScalarType>
    {
      std::vector<std::vector<ScalarType> > preconstructed_K_;
      int z_dim_;

    public:
      
      Hessian_Operator_Preconstructed(std::vector<std::vector<ScalarType> > & preconstructed_K): preconstructed_K_(preconstructed_K)
      { 
	z_dim_ = preconstructed_K.size();
      }
      
      //! Dtor
      ~Hessian_Operator_Preconstructed()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	for(int i = 0; i < z_dim_; i++)
	  {
	    RealT val = 0.0;
	    for(int j = 0; j < z_dim_; j++)
	      {
		val += preconstructed_K_[i][j]*(*x)(j);
	      }
	    y->Replace_Element(i,val);
	  }
      }
      
    };

  };

}

#endif
