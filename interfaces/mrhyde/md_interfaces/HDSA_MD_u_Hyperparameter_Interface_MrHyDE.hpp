/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_MD_U_HYPERPARAMETER_INTERFACE_MRHYDE_HPP
#define HDSA_MD_U_HYPERPARAMETER_INTERFACE_MRHYDE_HPP

#include "HDSA_MD_Data_Interface_MrHyDE.hpp"
#include "HDSA_MD_OUU_Data_Interface_MrHyDE.hpp"

template <class RealT>
class MD_u_Hyperparameter_Interface_MrHyDE : public HDSA::MD_u_Hyperparameter_Interface<RealT>
{

private:
  HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
  Teuchos::RCP<Teuchos::MpiComm<int>> comm_;

public:
  std::vector<std::vector<RealT>> Spatial_Domain_Bounds(void) const override
  {

    HDSA::Ptr<HDSA::MultiVector<ScalarT>> spatial_coords = data_interface_->Read_Spatial_Node_Data();
    int dim = spatial_coords->Number_of_Vectors();

    std::vector<std::vector<RealT>> vec;
    vec.resize(dim);
    for (int i = 0; i < dim; i++)
    {
      vec[i].resize(2);
      HDSA::Ptr<HDSA::Vector<RealT>> veci_hdsa = (*spatial_coords)[i];
      HDSA::Tpetra_Vector<RealT> eveci = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(*veci_hdsa);
      HDSA::Ptr<Tpetra::MultiVector<RealT>> tpetra_vec = eveci.getVector();

      // Compute the minimum and maximum
      RealT minVal = std::numeric_limits<RealT>::max();
      RealT maxVal = std::numeric_limits<RealT>::lowest();

      // Reduce to find the min and max
      int Number_of_Rows = tpetra_vec->getLocalLength();
      auto localView = tpetra_vec->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
      for (size_t i = 0; i < Number_of_Rows; ++i)
      {
        minVal = std::min(minVal, localView(i, 0));
        maxVal = std::max(maxVal, localView(i, 0));
      }

      // Use a global reduction to find the overall min and max
      RealT globalMin, globalMax;
      Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MIN, 1, &minVal, &globalMin);
      Teuchos::reduceAll(*comm_, Teuchos::REDUCE_MAX, 1, &maxVal, &globalMax);

      vec[i][0] = globalMin;
      vec[i][1] = globalMax;
    }

    return vec;
  }

  MD_u_Hyperparameter_Interface_MrHyDE(const Teuchos::RCP<Teuchos::MpiComm<int>> &comm, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface, const bool &is_transient, const bool &center_data = false, const bool &adapt_time_variance = false, const int &component_id = 0)
      : HDSA::MD_u_Hyperparameter_Interface<RealT>(is_transient, center_data, adapt_time_variance, component_id)
  {
    comm_ = comm;
    data_interface_ = data_interface;
  }
};

#endif
