#ifndef HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_HPP
#define HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Elliptic_z_Prior_Interface : public HDSA::MD_z_Prior_Interface<RealT> {

  private:
    RealT alpha_z_;

  public:
    MD_Elliptic_z_Prior_Interface(RealT & alpha_z): alpha_z_(alpha_z)
    { }

    virtual ~MD_Elliptic_z_Prior_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                                                                    
    // Pure virtual functions                                                                                                                                                                                                                                                         
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                                                                    

    virtual void Apply_E_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

    virtual void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

    virtual void Apply_M_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                
    // Virtual functions which must be implemented to enable the Hessian GEVP 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////  
    
    virtual void Apply_E_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
			       "Apply_E_z must be implemented to use the Hessian GEVD implementation" << std::endl);
    }

    virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Apply_E_z_Transpose must be implemented to use the Hessian GEVD implementation" << std::endl);
    }

    virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA_TEST_FOR_EXCEPTION( true, std::logic_error,
                               "Apply_M_z_Inverse must be implemented to use the Hessian GEVD implementation" << std::endl);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                                                                    
    // Implementation of base class Virtual functions                                                                                                                                                                                                                                 
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////  

    virtual void Apply_W_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = z_in.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = z_in.clone();
      Apply_E_z_Inverse_Transpose(*z_tmp1,z_in);
      Apply_M_z(*z_tmp2,*z_tmp1);
      Apply_E_z_Inverse(z_out,*z_tmp2);
      z_out.scale(alpha_z_);
    }

    virtual void Apply_W_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = z_in.clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = z_in.clone();
      Apply_E_z(*z_tmp1,z_in);
      Apply_M_z_Inverse(*z_tmp2,*z_tmp1);
      Apply_E_z_Transpose(z_out,*z_tmp2);
      z_out.scale(1.0/alpha_z_);
    }

    // Factorize W_z^{-1} = F*F^T and compute the matvec z_out = F*z_in                                                                                                                                                                                                               
    virtual void Apply_W_z_Inverse_Factor(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
    {

    }


  };

}

#endif

