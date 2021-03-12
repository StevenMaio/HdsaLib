#ifndef OPT_PROBLEM_OBJECTS_SYNTHETIC_TEST_HPP
#define OPT_PROBLEM_OBJECTS_SYNTHETIC_TEST_HPP

template <class RealT>
class Opt_Problem_Objects_synthetic_test : public HDSA::Opt_Problem_Objects<RealT> {

private:
  std::vector<RealT> a_;
  std::vector<RealT> b_;
  int z_dim_;

public:

  Opt_Problem_Objects_synthetic_test(std::vector<RealT> & a, std::vector<RealT> & b, int & z_dim): a_(a), b_(b), z_dim_(z_dim)
  { }

  Opt_Problem_Objects_synthetic_test(std::vector<RealT> & a, std::vector<RealT> & b,int & z_dim, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in, const HDSA::Ptr<const HDSA::Comm<int> > & comm):
    a_(a), b_(b), z_dim_(z_dim)
  { 
    HDSA::Opt_Problem_Objects<RealT>::theta = theta_in->Clone();
    HDSA::Opt_Problem_Objects<RealT>::theta->set(*theta_in);

    HDSA::Opt_Problem_Objects<RealT>::z = HDSA::makePtr<Std_Vector<RealT> >(z_dim);

    HDSA::Opt_Problem_Objects<RealT>::rs_obj = HDSA::makePtr<synthetic_test_RS_Objective<RealT> >(a,b);
  }

  virtual ~Opt_Problem_Objects_synthetic_test()
  { }

  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_synthetic_test<RealT> >(a_,b_,z_dim_,theta,comm);
    return OP_Objects;
  }

  void Solve_Optimization_Problem()
  {
    int a_dim = a_.size();
    HDSA::Opt_Problem_Objects<RealT>::z->zero();
    for(int k = 0; k < a_dim; k++)
      {
	HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(k,a_[k]/(1.0+b_[k]));
      }  
  }

  void Load_Optimal_Solution() 
  { 
    int a_dim = a_.size();
    HDSA::Opt_Problem_Objects<RealT>::z->zero();
    for(int k = 0; k < a_dim; k++)
      {
	HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(k,a_[k]/(1.0+b_[k]));
      }  
  }
  
  void Write_Optimal_Solution()
  { }

};


#endif
