#ifndef HDSA_MD_DATA_INTERFACE_HPP
#define HDSA_MD_DATA_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Data_Interface {

  public:
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt_;
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > D_;

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
      u_opt_ = Load_Optimal_u();
      z_opt_ = Load_Optimal_z();
      Z_ = Load_Z_Data();
      D_ = Load_D_Data();
    }

  };

}

#endif
