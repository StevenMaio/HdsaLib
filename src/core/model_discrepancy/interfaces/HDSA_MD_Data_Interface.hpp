#ifndef HDSA_MD_DATA_INTERFACE_HPP
#define HDSA_MD_DATA_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Data_Interface {

  public:
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt;
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z;
    HDSA::Ptr<HDSA::MultiVector<RealT> > D;

    MD_Data_Interface()
    { }

    virtual ~MD_Data_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                
    // Pure virtual functions                                                                                                                                                               

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////                                                                                                                  
    virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

    virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const = 0;

    void Load_Data(void)
    {
      u_opt = Load_Optimal_u();
      z_opt = Load_Optimal_z();
      Z = Load_Z_Data();
      D = Load_D_Data();
    }

  };

}

#endif
