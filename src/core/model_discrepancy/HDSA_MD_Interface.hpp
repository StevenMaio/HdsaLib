#ifndef HDSA_MD_INTERFACE_HPP
#define HDSA_MD_INTERFACE_HPP

namespace HDSA
{

template <class RealT>
class Model_Discrepancy_Interface {  

private:
  HDSA::Ptr<HDSA::Hessian_Inversion<RealT> > hess_invert_;

public:

  Model_Discrepancy_Interface()
  {  
    hess_invert_ = HDSA::makePtr<Model_Discrepancy_Hessian_Inversion<RealT> >(this);
  }

  virtual ~Model_Discrepancy_Interface()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

  virtual void Apply_L_Plus_beta_Identity_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT beta) const = 0;

  virtual void Apply_L_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const = 0;

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Implemention of inverse hessian via calling an iterative solver
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  virtual void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
  {
    hess_invert_->Apply_RS_Hessian_Inverse(z_out,z_in,z);
  }

  template <class ScalarType>
  class Model_Discrepancy_Hessian_Inversion : public HDSA::Hessian_Inversion<ScalarType>
  {
  private:
    const Model_Discrepancy_Interface<ScalarType>* md_interface_;

  public:
    Model_Discrepancy_Hessian_Inversion(const Model_Discrepancy_Interface<ScalarType>* md_interface): md_interface_(md_interface)
    { }
    
    ~Model_Discrepancy_Hessian_Inversion()
    {}
    
    void Apply_RS_Hessian(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const 
    {
      md_interface_->Apply_RS_Hessian(z_out,z_in,z);
    }

  };

};

}

#endif


