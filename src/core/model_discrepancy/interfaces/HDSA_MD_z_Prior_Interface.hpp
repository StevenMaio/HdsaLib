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

    // Compute samples from a mean zero Gaussian with covariance W_z^{-1}
    virtual void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> & samples) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms" << std::endl);
    }
        
    virtual void Apply_W_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Apply_W_z must be implemented to use the Hessian GEVD implementation" << std::endl);
    }

  };

}

#endif

