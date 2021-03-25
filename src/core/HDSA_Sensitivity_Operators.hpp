#ifndef HDSA_SENSITIVITY_OPERATORS_HPP
#define HDSA_SENSITIVITY_OPERATORS_HPP

// This class contains the various operators needed throughout the computation

namespace HDSA
{
  
  template <class RealT>
  class Sensitivity_Operators{
  protected:

    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects_;   
    HDSA::Ptr<HDSA::Nominal_Data<RealT> > Nom_;
    int theta_dim_, z_dim_;
    unsigned seed_;
    std::default_random_engine generator_;
    std::normal_distribution<RealT> distribution_;
    
  public:
    
    Sensitivity_Operators(const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, const HDSA::Ptr<HDSA::Nominal_Data<RealT> > & Nom, const std::vector<int> & Comm_Split_Ranks):
      OP_Objects_(OP_Objects), Nom_(Nom)
    { 
      bool time_seed = Nom->Get_parlist_sensitivity()->sublist("Formulation").get("System Time Seed",false);
      if(time_seed)
	{
	  seed_ = time(NULL)*(1+Comm_Split_Ranks[0]);
	}
      else
	{
	  seed_ = 239*(1+Comm_Split_Ranks[0]);
	}
      generator_.seed(seed_);
      distribution_ = std::normal_distribution<RealT>(0.0,1.0); 
      theta_dim_ = OP_Objects->theta->dimension();
      z_dim_ = OP_Objects->z->dimension();
    }
    
    virtual ~Sensitivity_Operators()
    { }
    
    HDSA::Ptr<HDSA::Vector<RealT> > Generate_Random_theta_Vector(void)
    {
      // Populate vectors with standard normal samples
      HDSA::Ptr<HDSA::Vector<RealT> > theta = OP_Objects_->theta->Clone();
      for(int l = 0; l < theta_dim_; l++)
	{
	  theta->Replace_Element(l,distribution_(generator_));
	}
      return theta;
    }

    HDSA::Ptr<HDSA::Vector<RealT> > Generate_Random_u_Vector(void)
    {
      // Populate vectors with standard normal samples
      HDSA::Ptr<HDSA::Vector<RealT> > u = OP_Objects_->u->Clone();
      for(int l = 0; l < u->dimension(); l++)
	{
	  u->Replace_Element(l,distribution_(generator_));
	}
      return u;
    }
    
    HDSA::Ptr<HDSA::Vector<RealT> > Generate_Random_z_Vector(void)
    {
      // Populate vectors with standard normal samples
      HDSA::Ptr<HDSA::Vector<RealT> > z = OP_Objects_->z->Clone();
      for(int l = 0; l < Nom_->Get_nonzero_z_dim(); l++)
	{
	  z->Replace_Element(z->Get_map_reduced_to_full(l),distribution_(generator_));
	}
      return z;
    }
    
    virtual void Apply_Sensitivity_Operator(HDSA::Ptr<HDSA::Vector<RealT> > & z, const HDSA::Ptr<HDSA::Vector<RealT> > & theta) = 0;
    
    virtual void Apply_Sensitivity_Operator_Transpose(HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::Vector<RealT> > & z) = 0;
    
    virtual void Construct_K_Test(void) = 0;

    virtual void Construct_Misfit_Hessian_Test(void) = 0;

    virtual void Construct_Regularization_Hessian_Test(void) = 0;

    virtual void Construct_B_Test(void) = 0;

    virtual void Construct_B_Transpose_Test(void) = 0;

  };

}

#endif
