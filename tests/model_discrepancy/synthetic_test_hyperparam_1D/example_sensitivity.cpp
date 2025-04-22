#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "../../../src/source_file.hpp"
#include "MD_Data_Interface_synthetic_test.hpp"
#include "MD_Opt_Prob_Interface_synthetic_test.hpp"
#include "MD_u_Hyperparameter_Interface_synthetic_test.hpp"
#include "MD_z_Hyperparameter_Interface_synthetic_test.hpp"

typedef double RealT;

int main(int argc, char *argv[])
{

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int>> comm = HDSA::makePtr<HDSA::Comm<int>>();

  int num_random_numbers = 1.e5;
  std::string random_number_file = "random_numbers.txt";
  HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator = HDSA::makePtr<HDSA::Random_Number_Generator<RealT>>(num_random_numbers, random_number_file);

  HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface = HDSA::makePtr<MD_Data_Interface_synthetic_test<RealT>>(random_number_generator);
  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> opt_prob_interface = HDSA::makePtr<MD_Opt_Prob_Interface_synthetic_test<RealT>>();
  HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface = HDSA::makePtr<MD_u_Hyperparameter_Interface_synthetic_test<RealT>>();
  HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> z_hyperparam_interface = HDSA::makePtr<MD_z_Hyperparameter_Interface_synthetic_test<RealT>>();

  // Define the mass matrix
  const int m = 51;
  auto map = Tpetra::createUniformContigMap<Tpetra::Map<>::local_ordinal_type, Tpetra::Map<>::global_ordinal_type>(m, comm->Get_Teuchos_Communicator());
  HDSA::Ptr<Tpetra::CrsMatrix<RealT,Tpetra::Map<>::local_ordinal_type,Tpetra::Map<>::global_ordinal_type> > M = HDSA::makePtr<Tpetra::CrsMatrix<RealT,Tpetra::Map<>::local_ordinal_type,Tpetra::Map<>::global_ordinal_type> >(map, 3); // 3 is the maximum number of non-zero entries per row
  RealT h = 1.0 / static_cast<RealT>(m - 1);
  Teuchos::Array<Tpetra::Map<>::global_ordinal_type> cols0 = {0, 1};
  Teuchos::Array<RealT> vals0 = {h / 3.0, h / 6.0};
  M->insertGlobalValues(0, cols0(), vals0());
  for (int i = 1; i < m - 1; ++i)
  {
    Teuchos::Array<Tpetra::Map<>::global_ordinal_type> cols = {i - 1, i, i + 1};
    Teuchos::Array<RealT> vals = {h / 6.0, 2.0 * h / 3.0, h / 6.0};
    M->insertGlobalValues(i, cols(), vals());
  }
  Teuchos::Array<Tpetra::Map<>::global_ordinal_type> colsm = {m - 2, m - 1};
  Teuchos::Array<RealT> valsm = {h / 6.0, h / 3.0};
  M->insertGlobalValues(m - 1, colsm(), valsm());
  M->fillComplete();

  return 0;
}
