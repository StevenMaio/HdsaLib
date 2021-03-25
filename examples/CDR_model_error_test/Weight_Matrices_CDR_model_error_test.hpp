#ifndef WEIGHT_MATRICES_CDR_MODEL_ERROR_TEST_HPP
#define WEIGHT_MATRICES_CDR_MODEL_ERROR_TEST_HPP

template <class RealT>
class Weight_Matrices_CDR_model_error_test : public HDSA::Weight_Matrices<RealT> {

private:
  int z_dim;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Mz;
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

public:

  Weight_Matrices_CDR_model_error_test(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity), parlist_sensitivity_(parlist_sensitivity)
  {  
    z_dim = 48;

    // Mz
    Mz = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(z_dim,z_dim);
    // read in data
    std::ifstream in_Mz("M_z.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Mz) {   
      for(int i = 0; i < z_dim; i++)
	{
	  for(int j = 0; j < z_dim; j++)
	    {
	      in_Mz >> val;
	      Mz->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from M_z.txt" << std::endl;
      } 
  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_CDR_model_error_test<RealT> >(parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    theta_out->set(*theta_in);
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > matvec = Mz->Multiply(*z_in,false);
    matvec->Write_Column_to_Vector(0,*z_out);
  }
    
};


#endif
