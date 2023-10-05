#ifndef HDSA_ITERATIVE_GEVP_HPP
#define HDSA_ITERATIVE_GEVP_HPP

namespace HDSA
{

  template <class RealT>
  class Iterative_GEVP{

  private:
    bool reduced_space_sen_;
    int num_sing_vals_, theta_dim_;
    HDSA::Ptr<HDSA::Vector<RealT> > theta_; 
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    HDSA::Ptr<const HDSA::Comm<int> > comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    std::vector<int> Comm_Split_Ranks_;
    int sample_index_;
    std::string gevp_formulation_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > evecs_;
    std::vector<RealT> evals_;

  public:
    Iterative_GEVP(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
		   const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
		   const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory, std::vector<int> & Comm_Split_Ranks, const int & sample_index)
      : theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory),
      Comm_Split_Ranks_(Comm_Split_Ranks), sample_index_(sample_index)
    {
      reduced_space_sen_ = parlist_sensitivity->sublist("Formulation").get("Reduced Space Sensitivities",true);
      num_sing_vals_ = parlist_sensitivity->sublist("GSVD Solver").get("Number of Singular Values",1); 
      gevp_formulation_ = parlist_sensitivity->sublist("GSVD Solver").get("GSVD Formulation","Jordan-Wielandt");
      theta_dim_ = theta->dimension();
    }
    
    ~Iterative_GEVP()
    {}

    void Compute(void)
    {
      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_);
      OP_Objects->Load_Optimal_Solution();
      bool enforce_z_zeros_ = parlist_sensitivity_->sublist("Formulation").get("Enforce z Zeros",false);
      if(enforce_z_zeros_)
	{
	  OP_Objects->z->Enforce_Zeros();
	}
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = weight_matrices_factory_->Construct_Weight_Matrices(theta_,comm_); 
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
      
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > A;
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > B;
      bool is_joint = false;

      if(gevp_formulation_ == "Gram")
	{
	  A = HDSA::makePtr<Gram_Operator<RealT> >(Sen_Op,weight_matrices);
	  B = HDSA::makePtr<Gram_Weight_Mat_Operator<RealT> >(weight_matrices);
	  evals_.resize(num_sing_vals_);
	  evecs_.resize(num_sing_vals_);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      evecs_[k] = theta_->Clone();
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  is_joint = true;
	  A = HDSA::makePtr<JW_Operator<RealT> >(Sen_Op,weight_matrices);
	  B = HDSA::makePtr<JW_Weight_Mat_Operator<RealT> >(weight_matrices);
	  evals_.resize(num_sing_vals_);
	  evecs_.resize(num_sing_vals_);
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      evecs_[k] = HDSA::makePtr<HDSA::Joint_Vector<RealT> >(OP_Objects->z,theta_);
	    }
	}

      RealT tol = parlist_sensitivity_->sublist("Anasazi GEVP").get("Eigen Solver Tolerance",1.0e-3);
      int numBlocks = parlist_sensitivity_->sublist("Anasazi GEVP").get("Number of Blocks",2);
      int blocksize = parlist_sensitivity_->sublist("Anasazi GEVP").get("Block Size",num_sing_vals_);
      std::string solver = parlist_sensitivity_->sublist("Anasazi GEVP").get("Eigen Solver","LOBPCG"); 
      HDSA::Linear_Algebra::Iterative_GEVP_Solver<RealT>(evecs_,evals_,A,B,tol,numBlocks,blocksize,solver,is_joint);
      Write_Solution(Nom, weight_matrices);
    }
    

    void Write_Solution(HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices)
    {
      std::vector<RealT> sensitivity_indices = std::vector<RealT>(theta_dim_,0.0);
      std::clock_t timer_SI = std::clock();
      // Estimate sensitivity indices
      std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > theta_sing_vecs;
      theta_sing_vecs.resize(num_sing_vals_);

      std::vector<RealT> sing_vals = std::vector<RealT>(num_sing_vals_,0.0);
      if(gevp_formulation_ == "Gram")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      sing_vals[k] = std::sqrt(evals_[k]);
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      sing_vals[k] = evals_[k];
	    }
	} 

      if(gevp_formulation_ == "Gram")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      theta_sing_vecs[k] = evecs_[k];
	    }
	}
      else if(gevp_formulation_ == "Jordan-Wielandt")
	{
	  for(int k = 0; k < num_sing_vals_; k++)
	    {
	      HDSA::Joint_Vector<RealT>* ex;
	      ex = dynamic_cast<HDSA::Joint_Vector<RealT>* >(&(*evecs_[k])); 
	      theta_sing_vecs[k]= ex->Get_Component_Vector_2();
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
	      sensitivity_indices[i] += std::pow(sing_vals[k],2)*std::pow(theta_sing_vecs[k]->dot(*theta_matvec),2);
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
	  fout << sing_vals[k] << std::setw(20);
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

    template <class ScalarType>
    class Gram_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      HDSA::Ptr<Sensitivity_Operators<ScalarType> > Sen_Op_;
      HDSA::Ptr<Weight_Matrices<ScalarType> > weight_matrices_;
      
    public:
      Gram_Operator(HDSA::Ptr<Sensitivity_Operators<ScalarType> > & Sen_Op, HDSA::Ptr<Weight_Matrices<ScalarType> > & weight_matrices): Sen_Op_(Sen_Op), weight_matrices_(weight_matrices)
      { }
      
      ~Gram_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const
      {
	HDSA::Ptr<HDSA::Vector<ScalarType> > z_vec = Sen_Op_->Generate_Random_z_Vector();
	HDSA::Ptr<HDSA::Vector<ScalarType> > z_matvec = Sen_Op_->Generate_Random_z_Vector();
	z_vec->zero();
	z_matvec->zero();
	Sen_Op_->Apply_Sensitivity_Operator(z_vec,x);
	weight_matrices_->Apply_z_Weight_Mat(z_matvec,z_vec);
	Sen_Op_->Apply_Sensitivity_Operator_Transpose(y,z_matvec);
      }
    };

    template <class ScalarType>
    class Gram_Weight_Mat_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      HDSA::Ptr<Weight_Matrices<ScalarType> > weight_matrices_;
      
    public:
      Gram_Weight_Mat_Operator(HDSA::Ptr<Weight_Matrices<ScalarType> > & weight_matrices): weight_matrices_(weight_matrices)
      { }
      
      ~Gram_Weight_Mat_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const
      {
	weight_matrices_->Apply_theta_Weight_Mat(y,x);
      }
    };

    template <class ScalarType>
    class JW_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      HDSA::Ptr<Sensitivity_Operators<ScalarType> > Sen_Op_;
      HDSA::Ptr<Weight_Matrices<ScalarType> > weight_matrices_;
      
    public:
      JW_Operator(HDSA::Ptr<Sensitivity_Operators<ScalarType> > & Sen_Op, HDSA::Ptr<Weight_Matrices<ScalarType> > & weight_matrices): Sen_Op_(Sen_Op), weight_matrices_(weight_matrices)
      { }
     
      ~JW_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const
      {
	const HDSA::Joint_Vector<RealT>* ex = dynamic_cast<const HDSA::Joint_Vector<RealT>* >(&(*x));
	HDSA::Joint_Vector<RealT>* ey = dynamic_cast<HDSA::Joint_Vector<RealT>* >(&(*y));
	HDSA::Ptr<HDSA::Vector<RealT> > ex_z = ex->Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > ex_theta = ex->Get_Component_Vector_2();
	HDSA::Ptr<HDSA::Vector<RealT> > ey_z = ey->Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > ey_theta = ey->Get_Component_Vector_2();
	HDSA::Ptr<HDSA::Vector<ScalarType> > z_vec = Sen_Op_->Generate_Random_z_Vector();
	z_vec->zero();
	Sen_Op_->Apply_Sensitivity_Operator(z_vec,ex_theta);
	weight_matrices_->Apply_z_Weight_Mat(ey_z,z_vec); 
	z_vec->zero();
	weight_matrices_->Apply_z_Weight_Mat(z_vec,ex_z);
	Sen_Op_->Apply_Sensitivity_Operator_Transpose(ey_theta,z_vec);
      }
    };
    
    template <class ScalarType>
    class JW_Weight_Mat_Operator : public HDSA::Linear_Operator<ScalarType>
    {
    private:
      HDSA::Ptr<Weight_Matrices<ScalarType> > weight_matrices_;
      
    public:
      JW_Weight_Mat_Operator(HDSA::Ptr<Weight_Matrices<ScalarType> > & weight_matrices): weight_matrices_(weight_matrices)
      { }
      
      ~JW_Weight_Mat_Operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<RealT> > & y, const HDSA::Ptr<HDSA::Vector<RealT> > & x) const
      {
	const HDSA::Joint_Vector<RealT>* ex = dynamic_cast<const HDSA::Joint_Vector<RealT>* >(&(*x));
	HDSA::Joint_Vector<RealT>* ey = dynamic_cast<HDSA::Joint_Vector<RealT>* >(&(*y));
	HDSA::Ptr<HDSA::Vector<RealT> > ex_z = ex->Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > ex_theta = ex->Get_Component_Vector_2();
	HDSA::Ptr<HDSA::Vector<RealT> > ey_z = ey->Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > ey_theta = ey->Get_Component_Vector_2();
	weight_matrices_->Apply_z_Weight_Mat(ey_z,ex_z);
	weight_matrices_->Apply_theta_Weight_Mat(ey_theta,ex_theta);
      }
    };

  };

}


#endif
