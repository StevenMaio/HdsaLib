#ifndef HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Elliptic_u_Prior_Interface : public HDSA::MD_u_Prior_Interface<RealT> {
  
  private:
    HDSA::Ptr<HDSA::Randomized_GSVD<RealT> > gsvd_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > sing_vecs_input_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > sing_vecs_output_;
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > sing_vals_;
    RealT alpha_u_;

  public:
    MD_Elliptic_u_Prior_Interface(RealT & alpha_u): alpha_u_(alpha_u)  
    { }

    virtual ~MD_Elliptic_u_Prior_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                            
    // Pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                      

    virtual void Apply_E_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

    virtual void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

    virtual void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                  
    // User interface to GSVD                                    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////     

    void Compute_E_u_Inverse_GSVD(int & num_sing_vals, int & oversampling, int & num_subspace_iters, HDSA::Vector<RealT> & u_vec)
    {
      gsvd_ = HDSA::makePtr<Elliptic_GSVD<RealT> >(this,u_vec);
      sing_vecs_output_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_sing_vals,u_vec);
      sing_vecs_input_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_sing_vals,u_vec);
      sing_vals_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_sing_vals,1);
      gsvd_->Compute_GSVD(*sing_vecs_input_,*sing_vecs_output_,*sing_vals_, num_sing_vals, oversampling, num_subspace_iters); 
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                             
    // Implementation of base class Virtual functions 
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    virtual void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT & beta) const 
    {
      u_out.zeros();
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > c = sing_vecs_output_->MatVec(u_in);
      for(int k = 0; k < c->numRows(); k++)
	{
	  RealT coeff = (*c)(k,0)*std::pow((*sing_vals_)(k,0),2.0)/(1.0+beta*std::pow((*sing_vals_)(k,0),2.0));
	  u_out.axpy(coeff,*(*sing_vecs_output_)[k]);
	}
      u_out.scale(alpha_u_);
    }

    virtual void Apply_W_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
    {
      u_out.zeros();
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > c = sing_vecs_output_->MatVec(u_in);
      for(int k = 0; k < c->numRows(); k++)
	{
	  RealT coeff = (*c)(k,0)*std::pow((*sing_vals_)(k,0),2.0);
	  u_out.axpy(coeff,*(*sing_vecs_output_)[k]);
	}
      u_out.scale(alpha_u_);
    }

    // Factorize W_u^{-1} = F*F^T and compute the matvec u_out = F*u_in
    virtual void Apply_W_u_Inverse_Factor(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
    {

    }

    // Factorize (W_u+scalar*M_u)^{-1} = F*F^T and compute the matvec u_out = F*u_in                                                                                                              
    virtual void Apply_W_u_Plus_scalar_M_u_Inverse_Factor(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT & beta) const
    {

    }

    template <class ScalarType>
    class Elliptic_GSVD : public HDSA::Randomized_GSVD<ScalarType>
    {
      
    private:
      HDSA::MD_Elliptic_u_Prior_Interface<ScalarType>* u_prior_interface_;

    public:
      Elliptic_GSVD(HDSA::MD_Elliptic_u_Prior_Interface<ScalarType>* u_prior_interface, const HDSA::Vector<ScalarType> & u): HDSA::Randomized_GSVD<ScalarType>(u,u)
      { 
	u_prior_interface_ = u_prior_interface;
      }
     
      virtual ~Elliptic_GSVD()
      {}

      void Apply_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const 
      {
	u_prior_interface_->Apply_E_u_Inverse(vec_out,vec_in);
      }
      
      void Apply_Operator_Transpose(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const 
      {
	u_prior_interface_->Apply_E_u_Inverse_Transpose(vec_out,vec_in);
      }
      
      void Apply_Input_Weighting_Operator_Inverse(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const
      {
	u_prior_interface_->Apply_M_u(vec_out,vec_in);
      }
      
      void Apply_Output_Weighting_Operator(HDSA::Vector<ScalarType> & vec_out, const HDSA::Vector<ScalarType> & vec_in) const 
      {
	u_prior_interface_->Apply_M_u(vec_out,vec_in);
      }

    };


  };

}

#endif

