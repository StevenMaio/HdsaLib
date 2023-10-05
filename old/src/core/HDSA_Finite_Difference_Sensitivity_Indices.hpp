#ifndef HDSA_FINITE_DIFFERENCE_SENSITIVITY_HPP
#define HDSA_FINITE_DIFFERENCE_SENSITIVITY_HPP

namespace HDSA
{

  template <class RealT>
  class Finite_Difference_Sensitivity_Indices
  {
  private:
    HDSA::Ptr<HDSA::Vector<RealT> > theta_;
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;
    const HDSA::Ptr<const HDSA::Comm<int> >  comm_;
    HDSA::Ptr<Opt_Problem_Objects<RealT> > OP_Objects_Factory_;
    HDSA::Ptr<Weight_Matrices<RealT> > weight_matrices_factory_;
    int theta_dim_;

  public:

    Finite_Difference_Sensitivity_Indices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity,
					  const HDSA::Ptr<const HDSA::Comm<int> > & comm, const HDSA::Ptr<Opt_Problem_Objects<RealT> > & OP_Objects_Factory,
					  const HDSA::Ptr<Weight_Matrices<RealT> > & weight_matrices_factory):
      theta_(theta), parlist_sensitivity_(parlist_sensitivity), comm_(comm), OP_Objects_Factory_(OP_Objects_Factory), weight_matrices_factory_(weight_matrices_factory)
    { }

    void Compute(void)
    {
      theta_dim_ = theta_->dimension();

      HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_,comm_); 
      HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = weight_matrices_factory_->Construct_Weight_Matrices(theta_,comm_); 
      comm_->barrier();
      OP_Objects->Solve_Optimization_Problem();
      comm_->barrier();
      HDSA::Ptr<HDSA::Vector<RealT> > z_nominal = OP_Objects->z->Clone();
      z_nominal->set(*OP_Objects->z);
      
      std::vector<RealT> FD_sensitivity_indices = std::vector<RealT>(theta_dim_);
      RealT h = parlist_sensitivity_->sublist("Formulation").get("Finite Difference Check Step",1.e-3);
      
      for(int k = 0; k < theta_dim_; k++)
	{
	  HDSA::Ptr<HDSA::Vector<RealT> > theta_1 = theta_->Clone();
	  HDSA::Ptr<HDSA::Vector<RealT> > theta_2 = theta_->Clone();
	  HDSA::Ptr<HDSA::Vector<RealT> > z_1 = z_nominal->Clone();
	  
	  theta_1->set(*theta_);
	  theta_1->Replace_Element(k,(*theta_)(k)+h);
	  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_k = OP_Objects_Factory_->Construct_Opt_Problem_Objects(theta_1,comm_); 
	  
	  // Solve optimization problem
	  if(comm_->getRank()==0)
	    {
	      std::cout << "Solving the optimization problem for perturbation " << k+1 << std::endl;
	    }
	  comm_->barrier();
	  OP_Objects_k->Solve_Optimization_Problem();
	  comm_->barrier();
	  
	  theta_1->basis(k);
	  weight_matrices->Apply_theta_Weight_Mat(theta_2,theta_1);
	  RealT theta_norm = std::sqrt(theta_2->dot(*theta_1));
	  
	  OP_Objects_k->z->axpy(-1.0,*z_nominal);
	  OP_Objects_k->z->scale(1.0/h);
	  weight_matrices->Apply_z_Weight_Mat(z_1,OP_Objects_k->z);
	  RealT z_norm = std::sqrt(OP_Objects_k->z->dot(*z_1));
	  
	  FD_sensitivity_indices[k] = z_norm/theta_norm;
	}
      
      comm_->barrier();
      Write_Solution(FD_sensitivity_indices);
    }

    void Write_Solution(std::vector<RealT> & FD_sensitivity_indices)
    {
      // Write solutions to text files
      std::string name;
      std::ofstream fout;
      name = "Finite_Difference_Sensitivity_Indices.txt";
      fout.open(name);
      for(int k = 0; k < theta_dim_; k++)
	{
	  fout << FD_sensitivity_indices[k] << std::setw(20);
	}
      fout.close();
    }

  };

}

#endif
