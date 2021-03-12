#ifndef OPT_PROBLEM_OBJECTS_CDR_MODEL_ERROR_TEST_HPP
#define OPT_PROBLEM_OBJECTS_CDR_MODEL_ERROR_TEST_HPP

template <class RealT>
class Opt_Problem_Objects_CDR_model_error_test : public HDSA::Opt_Problem_Objects<RealT> {

private:


public:

  Opt_Problem_Objects_CDR_model_error_test( )
  {

  }

  Opt_Problem_Objects_CDR_model_error_test(const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    // Load in data
    int u_dim = 50;
    int z_dim = 48;

    /*************************************************************************/
    /***************** READ DATA *********************************************/
    /*************************************************************************/ 
    // A
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > A = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(u_dim,u_dim);
    // read in data
    std::ifstream in_A("A.txt");          
    RealT val;
    // read the elements in the file into a vector  
    // test file open   
    if (in_A) {   
      for(int i = 0; i < u_dim; i++)
	{
	  for(int j = 0; j < u_dim; j++)
	    {
	      in_A >> val;
	      A->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from A.txt" << std::endl;
      }  

    // B_source
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B_source = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(u_dim,z_dim);
    // read in data
     std::ifstream in_B_source("B_source.txt");           
    // read the elements in the file into a vector  
    // test file open   
     if (in_B_source) {   
      for(int i = 0; i < u_dim; i++)
	{
	  for(int j = 0; j < z_dim; j++)
	    {
	      in_B_source >> val;
	      B_source->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from B_source.txt" << std::endl;
      }  

    // D
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > D = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(u_dim,u_dim);
    // read in data
     std::ifstream in_D("D.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in_D) {   
      for(int i = 0; i < u_dim; i++)
	{
	  for(int j = 0; j < u_dim; j++)
	    {
	      in_D >> val;
	      D->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from D.txt" << std::endl;
      }  

    // R
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(z_dim,z_dim);
    // read in data
     std::ifstream in_R("R.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in_R) {   
      for(int i = 0; i < z_dim; i++)
	{
	  for(int j = 0; j < z_dim; j++)
	    {
	      in_R >> val;
	      R->Replace_Element(i,j,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from R.txt" << std::endl;
      } 

    // data
    HDSA::Ptr<HDSA::Vector<RealT> > data = HDSA::makePtr<Std_Vector<RealT> >(u_dim);
    // read in data
    std::ifstream in_data("data.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in_data) {   
      for(int i = 0; i < u_dim; i++)
	{
	  in_data >> val;
	  data->Replace_Element(i,val);
	}   
    }
    else
      {
	std::cout << "Error loading the data from data.txt" << std::endl;
      } 
    
    // beta
    RealT beta;
    // read in data
     std::ifstream in_beta("beta.txt");           
    // read the elements in the file into a vector  
    // test file open   
    if (in_beta) { 
      in_beta >> beta;
    }
    else
      {
	std::cout << "Error loading the data from beta.txt" << std::endl;
      } 
    
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in->Clone();
    HDSA::Opt_Problem_Objects<RealT>::theta->set(*theta_in);

    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<Std_Vector<RealT> >(z_dim);
    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<CDR_model_error_test_RS_Objective<RealT> >(A,B_source,D,R,data, beta);
    HDSA::Opt_Problem_Objects<RealT>::fs_obj = HDSA::makePtr<CDR_model_error_test_FS_Objective<RealT> >(D,R,data,beta);
    HDSA::Opt_Problem_Objects<RealT>::con = HDSA::makePtr<CDR_model_error_test_Constraint<RealT> >(A,B_source);
    HDSA::Opt_Problem_Objects<RealT>::u = HDSA::makePtr<Std_Vector<RealT> >(u_dim);
  }

  virtual ~Opt_Problem_Objects_CDR_model_error_test()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_CDR_model_error_test<RealT> >(theta,comm);
    return OP_Objects;
  }

  void Solve_Optimization_Problem() 
  {
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("opt_solution.txt");          
    RealT value;
    int count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile) {   
      while ( inputFile >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal solution" << std::endl;
      } 

    // read in solution and write to Opt_Problem_Objects<RealT>::u
    std::ifstream inputFile_state("opt_state.txt");          
    count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile_state) {   
      while ( inputFile_state >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::u->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal state" << std::endl;
      } 
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("opt_solution.txt");          
    RealT value;
    int count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile) {   
      while ( inputFile >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal solution" << std::endl;
      } 

    // read in solution and write to Opt_Problem_Objects<RealT>::u
    std::ifstream inputFile_state("opt_state.txt");          
    count = 0;
    // read the elements in the file into a vector  
    // test file open   
    if (inputFile_state) {   
      while ( inputFile_state >> value ) {
	HDSA::Opt_Problem_Objects<RealT>::u->Replace_Element(count,value);
	count += 1;
      }
    }
    else
      {
	std::cout << "Error loading the optimal state" << std::endl;
      } 
  }
  
  void Write_Optimal_Solution()
  {

  }

};


#endif
