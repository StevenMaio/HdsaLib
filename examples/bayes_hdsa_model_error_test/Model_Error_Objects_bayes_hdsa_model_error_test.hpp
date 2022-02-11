#ifndef MODEL_ERROR_OBJECTS_BAYES_HDSA_MODEL_ERROR_TEST_HPP
#define MODEL_ERROR_OBJECTS_BAYES_HDSA_MODEL_ERROR_TEST_HPP

// Instantiation of Model_Error_Objects

template <class RealT>
class Model_Error_Objects_bayes_hdsa_model_error_test : public HDSA::Model_Error_Objects<RealT> {

private:
  int dim;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > L;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Linv;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gammainv;

public:

  Model_Error_Objects_bayes_hdsa_model_error_test(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
					   const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory):
    HDSA::Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory)
  {

  }

  virtual ~Model_Error_Objects_bayes_hdsa_model_error_test()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    dim = 51;

    // L
    L = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,dim);
    // read in data
    std::ifstream in_L("L.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_L) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < dim; j++)
	    {
	      in_L >> val;
	      L->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from L.txt" << std::endl;
      } 

    // Linv
    Linv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,dim);
    // read in data
    std::ifstream in_Linv("Linv.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in_Linv) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < dim; j++)
	    {
	      in_Linv >> val;
	      Linv->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Linv.txt" << std::endl;
      } 

    // Gamma
    Gamma = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,dim);
    // read in data
    std::ifstream in_Gamma("Gamma.txt");           
    // read the elements in the file into a vector  
    // test file open
    if (in_Gamma) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < dim; j++)
	    {
	      in_Gamma >> val;
	      Gamma->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Gamma.txt" << std::endl;
      }  

    // Gammainv
    Gammainv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(dim,dim);
    // read in data
    std::ifstream in_Gammainv("Gamma_inv.txt");           
    // read the elements in the file into a vector  
    // test file open
    if (in_Gammainv) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < dim; j++)
	    {
	      in_Gammainv >> val;
	      Gammainv->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Gamma_inv.txt" << std::endl;
      } 

  }
  
  void Apply_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > matvec = Gamma->Multiply(*z_in,false);
    matvec->Write_Column_to_Vector(0,*z_out);
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > matvec = Gammainv->Multiply(*z_in,false);
    matvec->Write_Column_to_Vector(0,*z_out);
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > matvec = L->Multiply(*u_in,false);
    matvec->Write_Column_to_Vector(0,*u_out);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > matvec = Linv->Multiply(*u_in,false);
    matvec->Write_Column_to_Vector(0,*u_out);
  }


};


#endif
