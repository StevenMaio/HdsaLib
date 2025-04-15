#ifndef HDSA_MD_DATA_INTERFACE_HPP
#define HDSA_MD_DATA_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Data_Interface {

  private:
    HDSA::Ptr<const HDSA::Vector<RealT> > u_opt_;
    HDSA::Ptr<const HDSA::Vector<RealT> > z_opt_;
    HDSA::Ptr<const HDSA::MultiVector<RealT> > Z_;
    HDSA::Ptr<const HDSA::MultiVector<RealT> > D_;
    HDSA::Ptr<HDSA::Vector<RealT> > data_shift_;
    bool is_u_opt_loaded_;
    bool is_z_opt_loaded_;
    bool is_Z_loaded_;
    bool is_D_loaded_;

  public:
    MD_Data_Interface()
    { 
      is_u_opt_loaded_ = false;
      is_z_opt_loaded_ = false;
      is_Z_loaded_ = false;
      is_D_loaded_ = false;
    }

    virtual ~MD_Data_Interface()
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // pure virtual functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
    
    virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const = 0;

    virtual HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // virtual functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    // Need to think about how to implement this in C++
    virtual void Separate_State_Components(HDSA::Vector<RealT> & u)
    { }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // accessor functions
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////// 

    HDSA::Ptr<const HDSA::Vector<RealT> > get_u_opt(void)
    {
      if(!is_u_opt_loaded_)
      {
        u_opt_ = Load_Optimal_u();
        data_shift_ = u_opt_->clone();
        data_shift_->zeros();
        is_u_opt_loaded_ = true;
      }
      return u_opt_;
    }

    HDSA::Ptr<const HDSA::Vector<RealT> > get_z_opt(void)
    {
      if(!is_z_opt_loaded_)
        {
          z_opt_ = Load_Optimal_z();
          is_z_opt_loaded_ = true;
        }
      return z_opt_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT> > get_Z(void)
    {
      if(!is_Z_loaded_)
	      {
          Z_ = Load_Z_Data();
          is_Z_loaded_ = true;
        }
      return Z_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT> > get_D(void)
    {
      if(!is_D_loaded_)
        {
          D_ = Load_D_Data();
          is_D_loaded_ = true;
        }
      return D_;
    }

  };

}

#endif
