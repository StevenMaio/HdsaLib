#ifndef HDSA_MD_DATA_INTERFACE_HPP
#define HDSA_MD_DATA_INTERFACE_HPP

#include "HDSA_Vector.hpp"
#include "HDSA_MultiVector.hpp"
#include "HDSA_Transient_Vector.hpp"
#include "HDSA_Transient_Vector_Const.hpp"

namespace HDSA
{

  template <class RealT>
  class MD_Data_Interface
  {

  private:
    HDSA::Ptr<HDSA::Vector<RealT>> u_opt_;
    HDSA::Ptr<HDSA::Vector<RealT>> z_opt_;
    HDSA::Ptr<HDSA::MultiVector<RealT>> Z_;
    HDSA::Ptr<HDSA::MultiVector<RealT>> D_;
    HDSA::Ptr<HDSA::Vector<RealT>> data_shift_;
    bool is_data_loaded_;

  public:
    MD_Data_Interface()
    {
      is_data_loaded_ = false;
    }

    virtual ~MD_Data_Interface()
    {
    }

    void Load_Data(void)
    {
      u_opt_ = Load_Optimal_u();
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
      RealT val = data_shift_->dot(*(*D_)[0]) / static_cast<RealT>(data_shift_->dimension());
      data_shift_->setScalar(val);
      for (int k = 0; k < D_->Number_of_Vectors(); k++)
      {
        (*D_)[k]->axpy(-1.0, *data_shift_);
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // pure virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_u(void) const = 0;

    virtual HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_z(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Z_Data(void) const = 0;

    virtual HDSA::Ptr<HDSA::MultiVector<RealT>> Load_D_Data(void) const = 0;

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // virtual functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    virtual HDSA::Ptr<HDSA::MultiVector<RealT>> Read_Spatial_Node_Data() const
    {
      HDSA::Ptr<HDSA::MultiVector<RealT>> spatial_nodes;
      std::cout << "Error: HDSA::MD_Data_Interface::Read_Spatial_Node_Data is a virtual function that was called but never implemented" << std::endl;
      return spatial_nodes;
    }

    virtual HDSA::Ptr<const HDSA::Vector<RealT>> Extract_State_Component(const HDSA::Vector<RealT> &u, int component_id) const
    {
      HDSA::Ptr<const HDSA::Vector<RealT>> u_component = HDSA::makePtrFromRef(u);
      return u_component;
    }

    virtual void Set_State_Component(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &u_component, int component_id) const
    {
      u.set(u_component);
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // functions to manage abstraction for stationary and transient problems
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    HDSA::Ptr<const HDSA::Vector<RealT>> Extract_State_Component(const HDSA::Vector<RealT> &u, int component_id, bool check_transient) const
    {
      HDSA::Ptr<const HDSA::Vector<RealT>> u_component;
      if (check_transient)
      {
        if (const Transient_Vector<RealT> *u_trans = dynamic_cast<const Transient_Vector<RealT> *>(&(u)))
        {
          int n_t = u_trans->Get_n_t();
          std::vector<HDSA::Ptr<const HDSA::Vector<RealT>>> u_component_trans;
          u_component_trans.resize(n_t);
          for (int k = 0; k < n_t; k++)
          {
            u_component_trans[k] = Extract_State_Component(*(*u_trans)[k], component_id);
          }
          u_component = HDSA::makePtr<Transient_Vector_Const<RealT>>(u_component_trans);
        }
        else
        {
          u_component = Extract_State_Component(u, component_id);
        }
      }
      else
      {
        u_component = Extract_State_Component(u, component_id);
      }
      return u_component;
    }

    void Set_State_Component(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &u_component, int component_id, bool check_transient) const
    {
      if (check_transient)
      {
        if (Transient_Vector<RealT> *u_trans = dynamic_cast<Transient_Vector<RealT> *>(&(u)))
        {
          const Transient_Vector<RealT> *u_component_trans = dynamic_cast<const Transient_Vector<RealT> *>(&(u_component));
          int n_t = u_trans->Get_n_t();
          for (int k = 0; k < n_t; k++)
          {
            Set_State_Component(*(*u_trans)[k], *(*u_component_trans)[k], component_id);
          }
        }
        else
        {
          Set_State_Component(u, u_component, component_id);
        }
      }
      else
      {
        Set_State_Component(u, u_component, component_id);
      }
    }

    //////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // accessor functions
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////

    HDSA::Ptr<const HDSA::Vector<RealT>> get_u_opt(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return u_opt_;
    }

    HDSA::Ptr<const HDSA::Vector<RealT>> get_z_opt(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return z_opt_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT>> get_Z(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return Z_;
    }

    HDSA::Ptr<const HDSA::MultiVector<RealT>> get_D(void)
    {
      if (!is_data_loaded_)
      {
        Load_Data();
      }
      return D_;
    }

    HDSA::Ptr<const HDSA::Vector<RealT>> get_data_shift(void)
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
