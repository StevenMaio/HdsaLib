#ifndef WEIGHT_MATRICES_DARCY_FLOW_DETERMINISTIC_HPP
#define WEIGHT_MATRICES_DARCY_FLOW_DETERMINISTIC_HPP

template <class RealT>
class Weight_Matrices_darcy_flow_deterministic : public HDSA::Weight_Matrices<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

public:

  Weight_Matrices_darcy_flow_deterministic(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity):
    HDSA::Weight_Matrices<RealT>(parlist_sensitivity)
  {
    parlist_sensitivity_ = parlist_sensitivity;
  }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  {
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_darcy_flow_deterministic<RealT> >(parlist_sensitivity_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const 
  {
    theta_out->set(*theta_in);
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    z_out->set(*z_in);
  }
    
};


#endif
