template <class RealT>
class Weight_Matrices_ode_control : public HDSA::Weight_Matrices<RealT> {

private:

public:

  Weight_Matrices_ode_control(const HDSA::Ptr<Teuchos::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity)
  { }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    theta_out->set(*theta_in);
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    z_out->set(*z_in);
  }
    
};
