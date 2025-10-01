#ifndef HDSA_MD_Z_HYPERPARAMETER_INTERFACE_MRHYDE_HPP
#define HDSA_MD_Z_HYPERPARAMETER_INTERFACE_MRHYDE_HPP

#include "HDSA_MD_Data_Interface_MrHyDE.hpp"
#include "HDSA_MD_OUU_Data_Interface_MrHyDE.hpp"

template <class RealT>
class MD_z_Hyperparameter_Interface_MrHyDE : public HDSA::MD_z_Hyperparameter_Interface<RealT>
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
      HDSA_Tpetra_Vector<RealT> eveci = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(*veci_hdsa);
      HDSA::Ptr<Tpetra::MultiVector<RealT>> tpetra_vec = eveci.getVector();

      // Compute the minimum and maximum
      RealT minVal = std::numeric_limits<RealT>::max();
      RealT maxVal = std::numeric_limits<RealT>::lowest();

      // Reduce to find the min and max
      int numRows = tpetra_vec->getLocalLength();
      auto localView = tpetra_vec->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
      for (size_t i = 0; i < numRows; ++i)
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

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const override
  {
    if (const MD_Data_Interface_MrHyDE<RealT> *data_interface_myhyde = dynamic_cast<const MD_Data_Interface_MrHyDE<RealT> *>(&(*data_interface_)))
    {
      data_interface_myhyde->State_Solve(u, z);
    }
    else if (const MD_OUU_Data_Interface_MrHyDE<RealT> *ouu_data_interface_myhyde = dynamic_cast<const MD_OUU_Data_Interface_MrHyDE<RealT> *>(&(*data_interface_)))
    {
      std::vector<HDSA::Ptr<MD_Data_Interface_MrHyDE<RealT>>> data_interface_myhyde_std = ouu_data_interface_myhyde->Get_Data_Interface_MrHyDE();
      data_interface_myhyde_std[0]->State_Solve(u, z);
    }
    else
    {
      std::cout << "Error: MD_z_Hyperparameter_Interface_MrHyDE was unable to execute the requested State_Solve" << std::endl;
    }
  }

  MD_z_Hyperparameter_Interface_MrHyDE(const Teuchos::RCP<Teuchos::MpiComm<int>> &comm, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(z_type, num_state_solves)
  {
    comm_ = comm;
    data_interface_ = data_interface;
  }

  MD_z_Hyperparameter_Interface_MrHyDE(const Teuchos::RCP<Teuchos::MpiComm<int>> &comm, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface, int seed, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(seed, z_type, num_state_solves)
  {
    comm_ = comm;
    data_interface_ = data_interface;
  }

  MD_z_Hyperparameter_Interface_MrHyDE(const Teuchos::RCP<Teuchos::MpiComm<int>> &comm, const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator, const std::string &z_type, const int &num_state_solves = 0) : HDSA::MD_z_Hyperparameter_Interface<RealT>(random_number_generator, z_type, num_state_solves)
  {
    comm_ = comm;
    data_interface_ = data_interface;
  }
};

#endif
