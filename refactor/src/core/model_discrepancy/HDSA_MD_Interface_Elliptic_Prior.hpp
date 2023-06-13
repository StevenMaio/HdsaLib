#ifndef HDSA_MD_INTERFACE_ELLIPTIC_PRIOR_HPP
#define HDSA_MD_INTERFACE_ELLIPTIC_PRIOR_HPP

namespace HDSA
{

template <class RealT>
class Model_Discrepancy_Interface_Elliptic_Prior : public HDSA::Model_Discrepancy_Interface<RealT> {

public:
  HDSA::Ptr<HDSA::Randomized_GSVD<RealT> > gsvd_;
  HDSA::Ptr<HDSA::MultiVector<RealT> > U_;
  HDSA::Ptr<HDSA::MultiVector<RealT> > V_;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Sigma_;

  Model_Discrepancy_Interface_Elliptic_Prior()
  {  }

  virtual ~Model_Discrepancy_Interface_Elliptic_Prior()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Call "Compute_Elliptic_GSVD from the derived class constructor
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  void Compute_Elliptic_GSVD(int num_sing_vals, int oversampling, int num_subspace_iters, HDSA::Vector<RealT> & u_vec)
  {
    gsvd_ = HDSA::makePtr<Elliptic_GSVD<RealT> >(this,u_vec);
    U_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_sing_vals,u_vec);
    V_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_sing_vals,u_vec);
    Sigma_ = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(num_sing_vals,1);
    gsvd_->Compute_GSVD(*V_,*U_,*Sigma_, num_sing_vals, oversampling, num_subspace_iters); 
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis
  // Note that the elliptic operator may be symmetric and hence the _Transpose function is redundant, but this
  // ensures generality for problems where boundary conditions break the symmetry
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Apply_u_Elliptic_Operator_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Elliptic_Operator_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Mass_Mat(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_u_Mass_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis (from the base class)
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const =0;

  virtual void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Implementation of some of the base class pure virtual functions
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void Apply_L_Plus_beta_Identity_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT beta) const 
  {
    u_out.zeros();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > c = U_->MatVec(u_in);
    for(int k = 0; k < c->numRows(); k++)
      {
	RealT coeff = (*c)(k,0)*std::pow((*Sigma_)(k,0),2.0)/(1.0+beta*std::pow((*Sigma_)(k,0),2.0));
	u_out.axpy(coeff,*(*U_)[k]);
      }
  }

  void Apply_L_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const  
  {
    u_out.zeros();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > c = U_->MatVec(u_in);
    for(int k = 0; k < c->numRows(); k++)
      {
	u_out.axpy((*c)(k,0)*std::pow((*Sigma_)(k,0),2.0),*(*U_)[k]);
      }
  }

  template <class ScalarType>
    class Elliptic_GSVD : public HDSA::Randomized_GSVD<ScalarType>
    {
      
    private:
      Model_Discrepancy_Interface_Elliptic_Prior<RealT>* elliptic_prior_;

    public:
      Elliptic_GSVD(Model_Discrepancy_Interface_Elliptic_Prior<RealT>* elliptic_prior, const HDSA::Vector<RealT> & u): HDSA::Randomized_GSVD<ScalarType>(u,u)
      { 
	elliptic_prior_ = elliptic_prior;
      }
     
      ~Elliptic_GSVD()
      {}

      void Apply_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
      {
	elliptic_prior_->Apply_u_Elliptic_Operator_Inverse(vec_out,vec_in);
      }
      
      void Apply_Operator_Transpose(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
      {
	elliptic_prior_->Apply_u_Elliptic_Operator_Inverse_Transpose(vec_out,vec_in);
      }
      
      void Apply_Input_Weighting_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
      {
	elliptic_prior_->Apply_u_Mass_Mat_Inverse(vec_out,vec_in);
      }
      
      void Apply_Input_Weighting_Operator_Inverse(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const
      {
	elliptic_prior_->Apply_u_Mass_Mat(vec_out,vec_in);
      }
      
      void Apply_Output_Weighting_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const 
      {
	vec_out.set(vec_in);
      }

    };

};

}

#endif


