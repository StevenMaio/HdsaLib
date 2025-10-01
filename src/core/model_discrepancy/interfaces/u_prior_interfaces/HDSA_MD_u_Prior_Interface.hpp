#ifndef HDSA_MD_U_PRIOR_INTERFACE_HPP
#define HDSA_MD_U_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_u_Prior_Interface
  {

  private:
  public:
    MD_u_Prior_Interface()
    {
    }

    virtual ~MD_u_Prior_Interface()
    {
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual void Apply_M_u(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const = 0;

    virtual void Apply_W_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in) const = 0;

    virtual void Apply_W_u_Plus_scalar_M_u_Inverse(HDSA::Vector<RealT> &u_out, const HDSA::Vector<RealT> &u_in, const RealT &scalar) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Virtual functions which must be implemented to enable posterior sampling
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    // Compute samples from a mean zero Gaussian with covariance W_u^{-1}
    virtual void Sample_with_Covariance_W_u_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "MD_u_Prior_Interface::Sample_with_Covariance_W_u_Inverse must be implemented to use sampling algorithms" << std::endl);
    }

    // Compute samples from a mean zero Gaussian with covariance W_u^{-1}
    virtual void Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(HDSA::MultiVector<RealT> &samples, const RealT &scalar) const
    {
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse must be implemented to use sampling algorithms" << std::endl);
    }
  };

}

#endif
