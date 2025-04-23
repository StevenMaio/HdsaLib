#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_HPP

namespace HDSA
{

  template <class RealT>
  class MD_z_Hyperparameter_Interface
  {

  private:
    std::string z_type_;
    int num_state_solves_;
    RealT discrepancy_percent_z_variation_;

    RealT alpha_z_;
    RealT beta_z_;
    RealT beta_t_;

    RealT z1_norm_sq_;
    std::vector<RealT> z_pert_norm_sq_;

  public:
    virtual HDSA::Ptr<HDSA::Vector<RealT>> Load_Spatial_Node_Data(void) const
    {
      std::cout << "Load_Spatial_Node_Data is required for hyperparameter algorithm-based initialization" << std::endl;
      HDSA::Ptr<HDSA::Vector<RealT> > vec;
      return vec;
    }

    virtual HDSA::Ptr<HDSA::Vector<RealT>> Load_Time_Node_Data(void) const
    {
      std::cout << "Load_Time_Node_Data is required for hyperparameter algorithm-based initialization" << std::endl;
      HDSA::Ptr<HDSA::Vector<RealT> > vec;
      return vec;
    }

    virtual void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
    {
      std::cout << "State_Solve is required to estimate alpha_z using low-fidelity solves" << std::endl;
    }

    MD_z_Hyperparameter_Interface(const std::string &z_type, const int &num_state_solves = 0) : z_type_(z_type), num_state_solves_(num_state_solves)
    {
      if (!(z_type == "spatial field" || z_type == "transient vector" || z_type == "vector"))
      {
        std::cout << "Error in MD_z_Hyperparameter_Interface: The input z_type should be either 'spatial field' 'transient vector' or 'vector' " << std::endl;
      }

      discrepancy_percent_z_variation_ = 1.0;

      alpha_z_ = 0.0;
      beta_z_ = 0.0;
      beta_t_ = 0.0;
    }

    virtual ~MD_z_Hyperparameter_Interface()
    {}

    void Set_alpha_z(RealT & alpha_z_new)
    {
      alpha_z_ = alpha_z_new;
    }

    void Set_beta_z(RealT & beta_z_new)
    {
      beta_z_ = beta_z_new;
    }

    void Set_beta_t(RealT & beta_t_new)
    {
      beta_t_ = beta_t_new;
    }
    
  };

}

#endif
