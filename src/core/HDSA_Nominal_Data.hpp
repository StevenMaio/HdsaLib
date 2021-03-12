#ifndef HDSA_NOM_DATA_HPP
#define HDSA_NOM_DATA_HPP

// Nominal_Data stores data to efficient pass through sensitivity computation
namespace HDSA
{

  template <class RealT>
  class Nominal_Data{

  private:
    bool first_order_regularization_update_;
    bool second_order_regularization_update_;
    int update_rank_;
    RealT first_order_update_alpha_;
    RealT eig_val_min_;
    std::vector<RealT> second_order_update_eig_vals_;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > second_order_update_eig_vecs_;
    HDSA::Ptr<HDSA::Vector<RealT> > grad_nominal_;
    RealT coeff_;

    int nonzero_z_dim_, theta_dim_, z_dim_;
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

  public:

    Nominal_Data(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects): parlist_sensitivity_(parlist_sensitivity)
    {
      z_dim_ = OP_Objects->z->dimension();
      theta_dim_ = OP_Objects->theta->dimension();
      nonzero_z_dim_ = OP_Objects->z->Get_nonzero_dim();

      first_order_regularization_update_ = parlist_sensitivity_->sublist("Regularization Update").get("1st Order Update",false);
      second_order_regularization_update_ = parlist_sensitivity_->sublist("Regularization Update").get("2nd Order Update",false);
      bool reduced_space_sen = parlist_sensitivity_->sublist("Formulation").get("Reduced Space Sensitivities",true);
      if(first_order_regularization_update_)
	{
	  if(!reduced_space_sen)
	    {
	      std::cout << "Error: 1st Order Update = true requires Reduced Space Sensitivities = true" << std::endl;
	    }
	  first_order_update_alpha_ = parlist_sensitivity_->sublist("Regularization Update").get("alpha",1.0);
	  eig_val_min_ = parlist_sensitivity_->sublist("Regularization Update").get("Target Minimum Eigenvalue",1.e-6);
	}
      
      if(second_order_regularization_update_)
	{
	  if(!reduced_space_sen)
	    {
	      std::cout << "Error: 2nd Order Update = true requires Reduced Space Sensitivities = true" << std::endl;
	    }
	  update_rank_ = parlist_sensitivity_->sublist("Regularization Update").get("2nd Order Update Rank",1);
	  second_order_update_eig_vals_.resize(update_rank_);
	  second_order_update_eig_vecs_.resize(update_rank_);
	  
	  // read in data
	  std::ifstream in_evals("Second_Order_Update_evals.txt");          
	  if (in_evals) 
	    {   
	      for(int i = 0; i < update_rank_; i++)
		{
		  in_evals >> second_order_update_eig_vals_[i];
		}
	    }
	  
	  // read in data
	  std::ifstream in_evecs("Second_Order_Update_evecs.txt"); 
	  RealT value = 0.0;
	  if (in_evecs) 
	    {   
	      for(int i = 0; i < update_rank_; i++)
		{
		  second_order_update_eig_vecs_[i] = OP_Objects->z->Clone();
		  for(int j = 0; j < z_dim_; j++)
		    {
		      in_evecs >> value;
		      second_order_update_eig_vecs_[i]->Replace_Element(j,value);
		    }
		  second_order_update_eig_vecs_[i]->scale(1.0/second_order_update_eig_vecs_[i]->norm());
		}
	    }
	  
	}
	  
    }

    // Regularization update functions
    void Set_grad_nominal(HDSA::Ptr<HDSA::Vector<RealT> > & grad_nominal)
    {
      grad_nominal_ = grad_nominal;
      coeff_ = 1.0/(first_order_update_alpha_*grad_nominal_->norm());
    }

    void Apply_Regularization_Update(const HDSA::Ptr<HDSA::Vector<RealT> > & matvec, const HDSA::Ptr<const HDSA::Vector<RealT> > & vec)
    {
      // compute regularization update applied to vec and add it to matvec
      if(first_order_regularization_update_)
	{
	  RealT update_coeff = coeff_*(vec->dot(*grad_nominal_));
	  matvec->axpy(update_coeff,*grad_nominal_);
	}
	    
      if(second_order_regularization_update_)
	{
	  for(int i = 0; i < update_rank_; i++)
	    {
	      RealT update_coeff = (eig_val_min_ + std::abs(second_order_update_eig_vals_[i]))*(vec->dot(*second_order_update_eig_vecs_[i]));
	      matvec->axpy(update_coeff,*second_order_update_eig_vecs_[i]);
	    }
	}
    }

    // Accessor functions
    int Get_nonzero_z_dim(void)
    {
      return nonzero_z_dim_;
    }

    int Get_theta_dim(void)
    {
      return theta_dim_;
    }

    int Get_z_dim(void)
    {
      return z_dim_;
    }

    HDSA::Ptr<HDSA::ParameterList> Get_parlist_sensitivity(void)
    {
      return parlist_sensitivity_;
    }

  };

}

#endif
