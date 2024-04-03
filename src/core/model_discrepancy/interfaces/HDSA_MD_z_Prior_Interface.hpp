#ifndef HDSA_MD_Z_PRIOR_INTERFACE_HPP
#define HDSA_MD_Z_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_z_Prior_Interface {

  private:

  public:
    MD_z_Prior_Interface()
    { }

    virtual ~MD_z_Prior_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                            
    // Pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                               

    virtual void Apply_W_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                
    // Virtual functions which must be implemented to enable posterior sampling and the Hessian GEVP                                                                                                                                                                                                     
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////  
    
    // Factorize W_z^{-1} = F*F^T and compute the matvec z_out = F*z_in
    virtual void Apply_W_z_Inverse_Factor(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
			       "Apply_W_z_Inverse_Factor must be implemented to use sampling algorithms" << std::endl);
    }

    virtual void Apply_W_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Apply_W_z must be implemented to use the Hessian GEVD implemented" << std::endl);
    }

  };

}

#endif

