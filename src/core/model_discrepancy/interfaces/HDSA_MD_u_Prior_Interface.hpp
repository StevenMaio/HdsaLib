#ifndef HDSA_MD_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_u_Prior_Interface {

  private:

  public:
    MD_u_Prior_Interface()
    { }

    virtual ~MD_u_Prior_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                            
    // Pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                               

    virtual void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

    virtual void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT & beta) const = 0;

    virtual void Apply_W_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                
    // Virtual functions which must be implemented to enable posterior sampling                                                                                                                                                                                                     
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////  
    
    // Factorize W_u^{-1} = F*F^T and compute the matvec u_out = F*u_in
    virtual void Apply_W_u_Inverse_Factor(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const 
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
			       "The method Apply_W_u_Inverse_Factor was called, but not implemented" << std::endl);
    }

    // Factorize (W_u+scalar*M_u)^{-1} = F*F^T and compute the matvec u_out = F*u_in                                                                                                                                                                                                           
    virtual void Apply_W_u_Plus_scalar_M_u_Inverse_Factor(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in, const RealT & beta) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Apply_W_u_Plus_scalar_M_u_Inverse_Factor must be implemented to use the sampling algorithms" << std::endl);
    }

  };

}

#endif

