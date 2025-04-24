#ifndef HDSA_MD_DATA_INTERFACE_HPP
#define HDSA_MD_DATA_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_Data_Interface {

  private:
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt_;
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z_;
    HDSA::Ptr<HDSA::MultiVector<RealT> > D_;
    HDSA::Ptr<HDSA::Vector<RealT> > data_shift_;
    bool is_data_loaded_;

  public:
    MD_Data_Interface()
    {  
      is_data_loaded_ = false;
    }

    virtual ~MD_Data_Interface()
    { }

    void Load_Data(void)
    {
      u_opt_ = u_opt_ = Load_Optimal_u();
      z_opt_ = Load_Optimal_z();
      Z_ = Load_Z_Data();
      D_ = Load_D_Data();
      data_shift_ = u_opt_->clone();
      data_shift_->zeros();
      is_data_loaded_ = true;

    }

    void Center_Data()
    {
        data_shift_->setScalar(1.0);
        RealT val = data_shift_->dot(*(*D_)[0])/static_cast<RealT>(data_shift_->dimension());
        data_shift_->setScalar(val);
        for (int k = 0; k < D_->Number_of_Vectors(); k++)
        {
          (*D_)[k]->axpy(-1.0,*data_shift_);
        }
    }

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
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return u_opt_;
    }

    HDSA::Ptr<const HDSA::Vector<RealT> > get_z_opt(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return z_opt_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT> > get_Z(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return Z_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT> > get_D(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return D_;
    }

    HDSA::Ptr<const HDSA::Vector<RealT> > get_data_shift(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return data_shift_;
    }

  };

}

#endif
