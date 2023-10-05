template <class RealT>
class Weight_Matrices_synthetic_test : public HDSA::Weight_Matrices<RealT> {

private:

public:

  Weight_Matrices_synthetic_test(const HDSA::Ptr<Teuchos::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity)
  { }

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
    for(int k = 0; k < dim; k++)
      {
	RealT val = (*z_in)(k)*static_cast<RealT>(k+2);
	z_out->Replace_Element(k,val);
      }
  }
    
};
