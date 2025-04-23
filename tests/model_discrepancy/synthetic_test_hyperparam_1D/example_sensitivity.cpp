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

  HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface = HDSA::makePtr<MD_Data_Interface_synthetic_test<RealT>>(random_number_generator,comm);
  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT>> opt_prob_interface = HDSA::makePtr<MD_Opt_Prob_Interface_synthetic_test<RealT>>(comm);
  HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> u_hyperparam_interface = HDSA::makePtr<MD_u_Hyperparameter_Interface_synthetic_test<RealT>>();
  HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> z_hyperparam_interface = HDSA::makePtr<MD_z_Hyperparameter_Interface_synthetic_test<RealT>>();

  HDSA::Ptr<MD_Opt_Prob_Interface_synthetic_test<RealT> > opt_prob_interface_st = HDSA::dynamicPtrCast<MD_Opt_Prob_Interface_synthetic_test<RealT> >(opt_prob_interface);
  HDSA::Ptr<HDSA::Sparse_Matrix<RealT> > M = opt_prob_interface_st->Get_Mass_Matrix();
  HDSA::Ptr<HDSA::Sparse_Matrix<RealT> > S = opt_prob_interface_st->Get_Stiffness_Matrix();

  HDSA::Ptr<HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT> > u_prior_interface = HDSA::makePtr<HDSA::MD_Numeric_Laplacian_u_Prior_Interface<RealT> >(S,M,data_interface,u_hyperparam_interface); 
  HDSA::Ptr<HDSA::MD_Numeric_Laplacian_z_Prior_Interface<RealT> > z_prior_interface = HDSA::makePtr<HDSA::MD_Numeric_Laplacian_z_Prior_Interface<RealT> >(S,M,data_interface,z_hyperparam_interface,u_prior_interface); 
  
  return 0;
}
