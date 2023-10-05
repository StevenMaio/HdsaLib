#ifndef BAYES_MODEL_ERROR_OBJECTS_BAYES_HDSA_MODEL_ERROR_TEST_HPP
#define BAYES_MODEL_ERROR_OBJECTS_BAYES_HDSA_MODEL_ERROR_TEST_HPP

// Instantiation of Bayes_Model_Error_Objects

template <class RealT>
class Bayes_Model_Error_Objects_bayes_hdsa_model_error_test : public HDSA::Bayes_Model_Error_Objects<RealT> {

private:
  int dim;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > L;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Linv;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gamma;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Gammainv;

public:

  Bayes_Model_Error_Objects_bayes_hdsa_model_error_test(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
							const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory):
    HDSA::Bayes_Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory)
  {

  }

  virtual ~Bayes_Model_Error_Objects_bayes_hdsa_model_error_test()
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

  void Apply_Sqrt_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    int n = z_in->dimension();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(n);
    HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(Gamma, V, S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Vt_z = V->Multiply(z_in, true);
    for(int k = 0; k < n; k++)
      {
	RealT val = (*Vt_z)(k,0)*std::sqrt((*S)(k));
	Vt_z->Replace_Element(k,0,val);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_out_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    V->Multiply(z_out_mat, Vt_z); 
    z_out_mat->Write_Column_to_Vector(0,z_out);
  }

  void Apply_Sqrt_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    int n = z_in->dimension();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(n);
    HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(Gammainv, V, S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Vt_z = V->Multiply(z_in, true);
    for(int k = 0; k < n; k++)
      {
	RealT val = (*Vt_z)(k,0)*std::sqrt((*S)(k));
	Vt_z->Replace_Element(k,0,val);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z_out_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,1);
    V->Multiply(z_out_mat, Vt_z); 
    z_out_mat->Write_Column_to_Vector(0,z_out);
  }
  
  void Apply_Sqrt_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    int m = u_in->dimension();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,m);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(m);
    HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(Linv, V, S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Vt_u = V->Multiply(u_in, true);
    for(int k = 0; k < m; k++)
      {
	RealT val = (*Vt_u)(k,0)*std::sqrt((*S)(k));
	Vt_u->Replace_Element(k,0,val);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > u_out_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,1);
    V->Multiply(u_out_mat, Vt_u); 
    u_out_mat->Write_Column_to_Vector(0,u_out);
  }
  
  void Apply_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    int m = u_in->dimension();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,m);
    for(int i = 0; i < m; i++)
      {
	for(int j = 0; j < m; j++)
	  {
	    A->Replace_Element(i,j,(*L)(i,j));
	  }
	A->Replace_Element(i,i,(*L)(i,i) + beta);
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,m);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(m);
    HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(A, V, S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Vt_u = V->Multiply(u_in, true);
    for(int k = 0; k < m; k++)
      {
	RealT val = (*Vt_u)(k,0)*(1.0/(*S)(k));
	Vt_u->Replace_Element(k,0,val);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > u_out_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,1);
    V->Multiply(u_out_mat, Vt_u); 
    u_out_mat->Write_Column_to_Vector(0,u_out);
  }

  void Apply_Sqrt_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    int m = u_in->dimension();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,m);
    for(int i = 0; i < m; i++)
      {
	for(int j = 0; j < m; j++)
	  {
	    A->Replace_Element(i,j,(*L)(i,j));
	  }
	A->Replace_Element(i,i,(*L)(i,i) + beta);
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,m);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(m);
    HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(A, V, S);

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Vt_u = V->Multiply(u_in, true);
    for(int k = 0; k < m; k++)
      {
	RealT val = (*Vt_u)(k,0)*(1.0/std::sqrt((*S)(k)));
	Vt_u->Replace_Element(k,0,val);
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > u_out_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m,1);
    V->Multiply(u_out_mat, Vt_u); 
    u_out_mat->Write_Column_to_Vector(0,u_out);
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    int num_vecs = 1;
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
