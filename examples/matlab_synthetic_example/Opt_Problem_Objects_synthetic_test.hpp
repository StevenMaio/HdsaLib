template <class RealT>
class Opt_Problem_Objects_synthetic_test : public HDSA_Matlab_Opt_Problem_Objects<RealT> {

private:
  int z_dim_;
  std::vector<RealT> a_;

public:

  Opt_Problem_Objects_synthetic_test(const int & z_dim, std::vector<RealT> & a): HDSA_Matlab_Opt_Problem_Objects<RealT>(z_dim), z_dim_(z_dim), a_(a)
  { }

  Opt_Problem_Objects_synthetic_test(const int & z_dim, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in,
                                     const HDSA::Ptr<const Teuchos::Comm<int> > & comm, std::vector<RealT> & a): HDSA_Matlab_Opt_Problem_Objects<RealT>(z_dim, theta_in, comm), z_dim_(z_dim), a_(a)
  { }


  HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > Construct_Opt_Problem_Objects(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const Teuchos::Comm<int> > & comm) 
  { 
    HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > OP_Objects =  HDSA::makePtr<Opt_Problem_Objects_synthetic_test<RealT> >(z_dim_,theta,comm,a_);
    return OP_Objects;
  }

  void Load_Optimal_Solution() 
  { 
    int a_dim = a_.size();
    HDSA::Opt_Problem_Objects<RealT>::z->zero();
    for(int k = 0; k < a_dim; k++)
      {
        HDSA::Opt_Problem_Objects<RealT>::z->Replace_Element(k,a_[k]);
      }  
  }
  


};
