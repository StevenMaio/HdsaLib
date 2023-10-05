template <class RealT>
class Opt_Problem_Objects_ode_control : public HDSA_Matlab_Opt_Problem_Objects<RealT> {

private:
  int z_dim_;

public:

  Opt_Problem_Objects_ode_control(const int & z_dim): HDSA_Matlab_Opt_Problem_Objects<RealT>(z_dim), z_dim_(z_dim)
  { }

  Opt_Problem_Objects_ode_control(const int & z_dim, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in,
				  const HDSA::Ptr<const Teuchos::Comm<int> > & comm): HDSA_Matlab_Opt_Problem_Objects<RealT>(z_dim, theta_in, comm), z_dim_(z_dim)
  { }


  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const Teuchos::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_ode_control<RealT> >(z_dim_,theta,comm);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    // read in solution and write to Opt_Problem_Objects<RealT>::z
    std::ifstream inputFile("control.txt");          
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
	std::cout << "Error loading the optimal control solution" << std::endl;
      }  
  }
  


};
