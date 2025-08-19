#ifndef HDSA_MD_VECTOR_Z_PRIOR_INTERFACE_HPP
#define HDSA_MD_VECTOR_Z_PRIOR_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Vector_z_Prior_Interface : public HDSA::MD_Scaled_z_Prior_Interface<RealT>
  {

  private:
  public:
    MD_Vector_z_Prior_Interface(RealT alpha_z): HDSA::MD_Scaled_z_Prior_Interface<RealT>(alpha_z)
    {
    }

    virtual ~MD_Vector_z_Prior_Interface()
    {
    }

    void Apply_M_z(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      z_out.set(z_in);
    }

    void Apply_W_z_Acute_Inverse(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      z_out.set(z_in);
    }

    void Sample_with_Covariance_W_z_Acute_Inverse(HDSA::MultiVector<RealT> &samples) const
    {
      samples.randomize_standard_normal();
    }

    void Apply_W_z_Acute(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in) const
    {
      z_out.set(z_in);
    }
  };

}

#endif
