#ifndef BAYES_MODEL_ERROR_OBJECTS_STOKES_HPP
#define BAYES_MODEL_ERROR_OBJECTS_STOKES_HPP

// Instantiation of Bayes_Model_Error_Objects

template <class RealT>
class Bayes_Model_Error_Objects_stokes : public HDSA::Bayes_Model_Error_Objects<RealT> {

private:


public:

  Bayes_Model_Error_Objects_stokes(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
							const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory):
    HDSA::Bayes_Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory)
  {

  }

  virtual ~Bayes_Model_Error_Objects_stokes()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
   

  }
  
  void Apply_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    z_out->set(*z_in);
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    z_out->set(*z_in);
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    u_out->set(*u_in);
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    u_out->set(*u_in);
  }

  void Apply_Sqrt_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    z_out->set(*z_in);
  }

  void Apply_Sqrt_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    z_out->set(*z_in);
  }
  
  void Apply_Sqrt_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    u_out->set(*u_in);
  }
  
  void Apply_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    u_out->set(*u_in);
  }

  void Apply_Sqrt_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    u_out->set(*u_in);
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->z);
    // read in data
    std::ifstream in_Z("Z.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Z) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Z >> val;
	      (*Z)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Z.txt" << std::endl;
      } 
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    // read in data
    std::ifstream in_Y("Y.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Y) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Y >> val;
	      (*Y)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Y.txt" << std::endl;
      } 
    return Y;
  }

};


#endif
