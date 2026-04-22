/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "HDSA_Stream.hpp"
#include "HDSA_Comm.hpp"
#include "HDSA_Ptr.hpp"
#include "HDSA_Std_Vector.hpp"
#include "Matrix_Sqrt_test.hpp"

typedef double RealT;

int main(int argc, char *argv[])
{

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int>> comm = HDSA::makePtr<HDSA::Comm<int>>();

  int m = 100;
  HDSA::Ptr<HDSA::Matrix_Sqrt<RealT>> mat_sqrt = HDSA::makePtr<Matrix_Sqrt_test<RealT>>(m);
  HDSA::Ptr<HDSA::Vector<RealT>> vec_in = HDSA::makePtr<HDSA::Std_Vector<RealT>>(m, comm);
  vec_in->Randomize_Standard_Normal();
  HDSA::Ptr<HDSA::Vector<RealT>> vec_out_1 = HDSA::makePtr<HDSA::Std_Vector<RealT>>(m, comm);
  HDSA::Ptr<HDSA::Vector<RealT>> vec_out_2 = HDSA::makePtr<HDSA::Std_Vector<RealT>>(m, comm);
  HDSA::Ptr<HDSA::Vector<RealT>> vec_out_3 = HDSA::makePtr<HDSA::Std_Vector<RealT>>(m, comm);

  mat_sqrt->Matrix_Sqrt_Apply(*vec_out_1, *vec_in);
  mat_sqrt->Matrix_Sqrt_Apply(*vec_out_2, *vec_out_1);
  mat_sqrt->Apply(*vec_out_3, *vec_in);
  vec_out_2->Scaled_Plus(-1.0, *vec_out_3);
  std::cout << "Error in Matrix squart root = " << vec_out_2->Norm() << std::endl;

  return 0;
}
