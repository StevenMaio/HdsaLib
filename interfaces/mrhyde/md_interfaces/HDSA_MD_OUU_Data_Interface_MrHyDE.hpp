#ifndef HDSA_MD_OUU_DATA_INTERFACE_MRHYDE_HPP
#define HDSA_MD_OUU_DATA_INTERFACE_MRHYDE_HPP

#include "HDSA_MD_OUU_Data_Interface.hpp"
#include "HDSA_MD_Data_Interface_MrHyDE.hpp"

template <class RealT>
class MD_OUU_Data_Interface_MrHyDE : public HDSA::MD_OUU_Data_Interface<RealT>
{

private:
  std::vector<HDSA::Ptr<MD_Data_Interface_MrHyDE<RealT>>> data_interface_mrhyde_;
  HDSA::Ptr<ROL::SampleGenerator<RealT>> sampler_;
  int ens_size_;

public:
  MD_OUU_Data_Interface_MrHyDE(std::vector<HDSA::Ptr<MD_Data_Interface_MrHyDE<ScalarT>>> &data_interface_mrhyde, HDSA::Ptr<ROL::SampleGenerator<RealT>> &sampler, int ens_size) : HDSA::MD_OUU_Data_Interface<RealT>(ens_size), data_interface_mrhyde_(data_interface_mrhyde), sampler_(sampler), ens_size_(ens_size)
  {
  }

  virtual ~MD_OUU_Data_Interface_MrHyDE()
  {
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_us(int s) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> us_opt = data_interface_mrhyde_[s]->Load_Optimal_u();
    return us_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_z(void) const
  {
    HDSA::Ptr<HDSA::Vector<RealT>> z_opt = data_interface_mrhyde_[0]->Load_Optimal_z();
    return z_opt;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Z_Data(void) const
  {
    HDSA::Ptr<HDSA::MultiVector<RealT>> Z = data_interface_mrhyde_[0]->Load_Z_Data();
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Ds_Data(int s) const
  {
    std::vector<RealT> param = sampler_->getMyPoint(s);
    data_interface_mrhyde_[s]->Get_Parameter_Manager()->updateParams(param, "stochastic");
    HDSA::Ptr<HDSA::MultiVector<RealT>> Ds = data_interface_mrhyde_[s]->Load_D_Data();
    return Ds;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Read_Spatial_Node_Data() const
  {
    HDSA::Ptr<HDSA::MultiVector<RealT>> spatial_nodes = data_interface_mrhyde_[0]->Read_Spatial_Node_Data();
    return spatial_nodes;
  }

};
#endif
