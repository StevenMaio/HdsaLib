template <class RealT>
class HDSA_Matlab_Opt_Problem_Objects : public HDSA::Opt_Problem_Objects<RealT> {

private:
  int z_dim_;

public:

  HDSA_Matlab_Opt_Problem_Objects(const int & z_dim): z_dim_(z_dim)
  { }

  HDSA_Matlab_Opt_Problem_Objects(const int & z_dim, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in,const HDSA::Ptr<const Teuchos::Comm<int> > & comm): z_dim_(z_dim)
  { 
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in->Clone();
    HDSA::Opt_Problem_Objects<RealT>::theta->set(*theta_in);

    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<Std_Vector<RealT> >(z_dim);

    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<Matlab_RS_Objective<RealT> >();
  }

  virtual HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const Teuchos::Comm<int> > & comm) =0;

  void Solve_Optimization_Problem()
  {
    std::cout << "The code current does not support solving the optimization problem" << std::endl;
  }

  virtual void Load_Optimal_Solution() =0;
  
  void Write_Optimal_Solution()
  { 
    std::cout << "The code current does not support writing the solution of the optimization problem" << std::endl;
  }

};
