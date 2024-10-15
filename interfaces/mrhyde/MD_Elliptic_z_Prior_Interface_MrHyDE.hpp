#ifndef HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP
#define HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP

  template <class RealT>
  class MD_Elliptic_z_Prior_Interface_MrHyDE : public HDSA::MD_Elliptic_z_Prior_Interface<RealT> {

  private:

  public:
    MD_Elliptic_z_Prior_Interface_MrHyDE(RealT & alpha_z): HDSA::MD_Elliptic_z_Prior_Interface<RealT>(alpha_z)
    { }

    virtual ~MD_Elliptic_z_Prior_Interface_MrHyDE()
    { }

    void Apply_E_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {} //identical to u

    void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {} //identical u if opt var lives on mesh

    void Apply_M_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {} //identical to u


    // Compute samples from a mean zero Gaussian with covariance W_z^{-1}                                                                                                                                                              
    virtual void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> & samples) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "MD_Elliptic_z_Prior_Interface::Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms" << std::endl);
    }   
    
    virtual void Apply_E_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
			       "MD_Elliptic_z_Prior_Interface::Apply_E_z must be implemented to use the Hessian GEVP" << std::endl);
        //could call mass matrix
    }

    virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "MD_Elliptic_z_Prior_Interface::Apply_E_z_Transpose must be implemented to use the Hessian GEVP" << std::endl);
         //could call mass matrix
    }

    virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "MD_Elliptic_z_Prior_Interface::Apply_M_z_Inverse must be implemented to use the Hessian GEVP" << std::endl);
         // identity
    }

  };

#endif

