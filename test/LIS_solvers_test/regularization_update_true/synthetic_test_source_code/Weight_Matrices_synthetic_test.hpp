#ifndef WEIGHT_MATRICES_SYNTHETIC_TEST_HPP
#define WEIGHT_MATRICES_SYNTHETIC_TEST_HPP

template <class RealT>
class Weight_Matrices_synthetic_test : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

public:

  Weight_Matrices_synthetic_test(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity)
  { 
    parlist_sensitivity_ = parlist_sensitivity;
  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_synthetic_test<RealT> >(parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    int dim = theta_in->dimension();
    for(int k = 0; k < dim; k++)
      {
	RealT val = (*theta_in)(k)*static_cast<RealT>(k+1);
	theta_out->Replace_Element(k,val);
      }
  }
    
  // Precondition solve for theta_Weight_Mat inverse via solve P*A*P*y = P*b and x=P*y where we assume P=P^T, typically P is diagonal
  void Apply_theta_Weight_Mat_Preconditioner(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const
  { 
    theta_out->set(*theta_in);
    for(int k = 0; k < 5; k++)
      {
	RealT val = (*theta_in)(k)*std::sqrt(1.0/static_cast<RealT>(k+1));
	theta_out->Replace_Element(k,val);
      }
  }

  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    int dim = z_in->dimension();
    RealT val = 2*(*z_in)(1) + (*z_in)(2); 
    z_out->Replace_Element(0,val);
    for(int k = 1; k < dim-1; k++)
      {
	val = (*z_in)(k-1) + 2*(*z_in)(k) + (*z_in)(k+1); 
	z_out->Replace_Element(k,val);
      }
    val = (*z_in)(dim-1) + 2*(*z_in)(dim); 
    z_out->Replace_Element(dim-1,val);
  }
    
};


#endif
