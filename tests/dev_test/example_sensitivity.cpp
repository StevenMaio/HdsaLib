#include "Teuchos_GlobalMPISession.hpp"
#include <fstream>

#include "../../src/source_file.hpp"
#include "MD_Data_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_Opt_Prob_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_u_Prior_Interface_model_discrepancy_synthetic_test.hpp"
#include "MD_z_Prior_Interface_model_discrepancy_synthetic_test.hpp"

typedef double RealT;

int main(int argc, char *argv[]) {

  HDSA::nullstream bhs;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &bhs);
  HDSA::Ptr<const HDSA::Comm<int> > comm = HDSA::makePtr<HDSA::Comm<int> >();
 
  HDSA::Ptr<HDSA::MD_Data_Interface<RealT> > data_interface = HDSA::makePtr<MD_Data_Interface_model_discrepancy_synthetic_test<RealT> >();
  data_interface->Load_Data();

  HDSA::Ptr<HDSA::MD_Opt_Prob_Interface<RealT> > opt_prob_interface = HDSA::makePtr<MD_Opt_Prob_Interface_model_discrepancy_synthetic_test<RealT> >();

  HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT> > u_prior_interface = HDSA::makePtr<MD_u_Prior_Interface_model_discrepancy_synthetic_test<RealT> >();
  HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT> > z_prior_interface = HDSA::makePtr<MD_z_Prior_Interface_model_discrepancy_synthetic_test<RealT> >();

  HDSA::Ptr<HDSA::MD_Prior_Sampling<RealT> > prior_sampling = HDSA::makePtr<HDSA::MD_Prior_Sampling<RealT> >(data_interface,u_prior_interface,z_prior_interface);

  HDSA::Ptr<HDSA::MultiVector<RealT> > z = HDSA::makePtr<HDSA::MultiVector<RealT> >(3,*data_interface->z_opt_);
  HDSA::Ptr<HDSA::Vector<RealT> > z0 = (*z)[0];
  HDSA::Ptr<HDSA::Vector<RealT> > z1 = (*z)[1];
  HDSA::Ptr<HDSA::Vector<RealT> > z2 = (*z)[2];
  Std_Vector<RealT> z0_std = dynamic_cast<Std_Vector<RealT>&>(*z0);
  Std_Vector<RealT> z1_std = dynamic_cast<Std_Vector<RealT>&>(*z1);
  Std_Vector<RealT> z2_std = dynamic_cast<Std_Vector<RealT>&>(*z2);
  int m = z0->dimension();
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > x = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,1);
  for(int k = 0; k < m; k++)
    {
      x->Replace_Element(k,0,static_cast<RealT>(k)/static_cast<RealT>(m-1));
    }
  for(int k = 0; k < m; k++)
    {
      z0_std.Replace_Element(k,(*x)(k,0)+1.5);
      z1_std.Replace_Element(k,(*x)(k,0) + std::pow((*x)(k,0),2.0));
      z2_std.Replace_Element(k,(*x)(k,0) + std::pow((*x)(k,0),2.0) - 2.0);
    }

  int num_samples = 10;
  std::vector<HDSA::Ptr<HDSA::MultiVector<RealT> > > prior_samples = prior_sampling->Prior_Discrepancy_Samples(*z,num_samples);

  for(int i = 0; i < num_samples; i++)
    {
      std::string name = "delta_sample_" + std::to_string(i+1) + "_evaluated_at_z";
      prior_samples[i]->Write_to_File(name);
    }

  HDSA::Ptr<HDSA::MultiVector<RealT> > prior_samples_at_z_opt = prior_sampling->Prior_Discrepancy_Samples_at_z_opt(num_samples);
  std::string name = "delta_sample_evaluated_at_z_opt";
  prior_samples_at_z_opt->Write_to_File(name);
    

  return 0;
}
