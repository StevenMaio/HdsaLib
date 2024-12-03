#ifndef HDSA_PRIOR_FE_OP_LA_BASE_HPP
#define HDSA_PRIOR_FE_OP_LA_BASE_HPP

template <class RealT,
	  class LO=Tpetra::Map<>::local_ordinal_type,
	  class GO=Tpetra::Map<>::global_ordinal_type,
	  class Node=Tpetra::Map<>::node_type >
class MD_Prior_FE_Op_LA_Base{
  
private:

public:

  HDSA::Ptr<Tpetra::CrsMatrix<RealT,LO,GO,Node>> A_;
  
  MD_Prior_FE_Op_LA_Base(HDSA::Ptr<Tpetra::CrsMatrix<RealT,LO,GO,Node>> &A)
  {
    A_=A;
  }

  virtual ~MD_Prior_FE_Op_LA_Base()
  { }                                                                              

  void Apply_A_Inverse(HDSA::Vector<RealT> & x, const HDSA::Vector<RealT> & b) {
    RealT tol = 1.0E-10;
    std::string solver = "CG";
    bool verbose = false;
    HDSA::Ptr<A_Operator<RealT> > A_op = HDSA::makePtr<A_Operator<RealT> >(this);
    HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(x,b,*A_op,tol,solver,verbose); 
  }
  
  template <class ScalarType>
  class A_Operator : public HDSA::Linear_Operator<ScalarType>
  {
  private:
    const MD_Prior_FE_Op_LA_Base<ScalarType>* A_invert_;

  public:
    A_Operator(const MD_Prior_FE_Op_LA_Base<ScalarType>* A_invert): A_invert_(A_invert)
    {}
    
    ~A_Operator()
    {}
      
    void matvec(HDSA::Vector<ScalarType> & y, const HDSA::Vector<ScalarType> & x) const
    {
      const Tpetra_Vector<ScalarType> &ex = dynamic_cast<const Tpetra_Vector<ScalarType>&>(x);
      Tpetra_Vector<ScalarType> &ey = dynamic_cast<Tpetra_Vector<ScalarType>&>(y);
      A_invert_->A_->apply(*ex.getVector(),*ey.getVector()); 
    }

  };
};

#endif
