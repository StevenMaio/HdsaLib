#ifndef HDSA_MD_INTERFACE_HPP
#define HDSA_MD_INTERFACE_HPP

namespace HDSA
{

template <class RealT>
class Model_Discrepancy_Interface {  

public:

  Model_Discrepancy_Interface()
  {  }

  virtual ~Model_Discrepancy_Interface()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define on a problem-to-problem basis
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  virtual void Apply_Gamma_Mat_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

  virtual void Apply_L_Plus_beta_Identity_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT beta) const = 0;

  virtual void Apply_L_Mat_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

  virtual void Apply_Solution_Operator_z_Jacobian_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_RS_Hessian_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Misfit_Gradient(HDSA::Vector<RealT> & u_grad, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual void Apply_Misfit_Hessian(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & z) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

  virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

  virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const = 0;

};

}

#endif


