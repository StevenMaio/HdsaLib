#ifndef HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_MRHYDE_HPP
#define HDSA_MD_ELLIPTIC_U_PRIOR_INTERFACE_MRHYDE_HPP

template <class RealT,
	  class LO=Tpetra::Map<>::local_ordinal_type,
	  class GO=Tpetra::Map<>::global_ordinal_type,
	  class Node=Tpetra::Map<>::node_type >
class MD_Elliptic_u_Prior_Interface_MrHyDE : public HDSA::MD_Elliptic_u_Prior_Interface<RealT> {

  typedef Tpetra::CrsMatrix<ScalarT,LO,GO,Node>   LA_CrsMatrix;
  typedef Teuchos::RCP<LA_CrsMatrix>              matrix_RCP;
    
private:
    HDSA::Ptr<HDSA_Prior_FE_Op_MrHyDE_Interface<RealT> > prior_fe_op_;
    std::vector<matrix_RCP> E_u_;
    std::vector<HDSA::Ptr<MD_Prior_FE_Op_LA_Base<RealT> > > E_u_solver_;
  
public:
  MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, RealT & beta_u, HDSA::Ptr<HDSA_Prior_FE_Op_MrHyDE_Interface<RealT>> &prior_fe_op, HDSA::Ptr<HDSA::Vector<RealT> > &uvec, int &prior_num_sing_vals, int &prior_oversampling, int &prior_num_subspace_iter):
    HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u), prior_fe_op_(prior_fe_op)
  {
      Construct_Elliptic_Operator(beta_u,uvec,prior_num_sing_vals,prior_oversampling,prior_num_subspace_iter);
  }

  MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, RealT & beta_u, HDSA::Ptr<HDSA_Prior_FE_Op_MrHyDE_Interface<RealT>> &prior_fe_op, HDSA::Ptr<HDSA::Vector<RealT> > &uvec, int &prior_num_sing_vals, int &prior_oversampling, int &prior_num_subspace_iter, int & seed):
    HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u,seed),prior_fe_op_(prior_fe_op)
  {
    Construct_Elliptic_Operator(beta_u,uvec,prior_num_sing_vals,prior_oversampling,prior_num_subspace_iter);
  }

  MD_Elliptic_u_Prior_Interface_MrHyDE(RealT & alpha_u, RealT & beta_u, HDSA::Ptr<HDSA_Prior_FE_Op_MrHyDE_Interface<RealT>> &prior_fe_op, HDSA::Ptr<HDSA::Vector<RealT> > &uvec, int &prior_num_sing_vals, int &prior_oversampling,int &prior_num_subspace_iter, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > & random_number_generator):
    HDSA::MD_Elliptic_u_Prior_Interface<RealT>(alpha_u,random_number_generator),prior_fe_op_(prior_fe_op)
  {
    Construct_Elliptic_Operator(beta_u,uvec,prior_num_sing_vals,prior_oversampling,prior_num_subspace_iter);
  }
    
  virtual ~MD_Elliptic_u_Prior_Interface_MrHyDE()
  { }

  void Construct_Elliptic_Operator(RealT &beta_u, HDSA::Ptr<HDSA::Vector<RealT> > &uvec, int &prior_num_sing_vals, int &prior_oversampling,int &prior_num_subspace_iter) {
      int k = prior_fe_op_->M.size();
      E_u_.resize(k);
      E_u_solver_.resize(k);
      
      for(int i=0;i<k;i++) {
	Tpetra::RowMatrix<RealT, LO, GO, Node> &S_tmp = dynamic_cast<Tpetra::RowMatrix<RealT, LO, GO, Node> &> (*prior_fe_op_->S[i]);
	Teuchos::RCP<Tpetra::RowMatrix<RealT, LO, GO, Node> > E_u_tmp = prior_fe_op_->M[i]->add(beta_u,S_tmp,1.0,prior_fe_op_->M[i]->getDomainMap(),prior_fe_op_->M[i]->getRangeMap(),Teuchos::null);
	E_u_[i] = HDSA::dynamicPtrCast<Tpetra::CrsMatrix<RealT, LO, GO, Node> > (E_u_tmp);
	E_u_solver_[i] = HDSA::makePtr<MD_Prior_FE_Op_LA_Base<RealT>>(E_u_[i]);
      }

     HDSA::MD_Elliptic_u_Prior_Interface<RealT>::Compute_E_u_Inverse_GSVD(prior_num_sing_vals,prior_oversampling,prior_num_subspace_iter,*uvec);

  }

  void Apply_E_u_Inverse(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {
    const Vector_MrHyDE_State<RealT> &eu_in = dynamic_cast<const Vector_MrHyDE_State<RealT>&>(u_in);  
    Vector_MrHyDE_State<RealT> &eu_out = dynamic_cast<Vector_MrHyDE_State<RealT>&>(u_out);  
    for (int i=0; i<eu_in.mrhyde_state_vec.size(); ++i){
      HDSA::Ptr<HDSA::Vector<RealT> > veci_in = HDSA::makePtr<Tpetra_Vector<RealT> > (eu_in.mrhyde_state_vec[i][0], eu_in.Get_random_number_generator());
      HDSA::Ptr<HDSA::Vector<RealT> > veci_out = HDSA::makePtr<Tpetra_Vector<RealT> > (eu_out.mrhyde_state_vec[i][0],eu_out.Get_random_number_generator()); 
      E_u_solver_[i]->Apply_A_Inverse(*veci_out,*veci_in);
    }
  }

  void Apply_E_u_Inverse_Transpose(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {
    Apply_E_u_Inverse(u_out, u_in);
  }

  void Apply_M_u(HDSA::Vector<RealT> & u_out, const HDSA::Vector<RealT> & u_in) const {
    const Vector_MrHyDE_State<RealT> &eu_in = dynamic_cast<const Vector_MrHyDE_State<RealT>&>(u_in);
    Vector_MrHyDE_State<RealT> &eu_out = dynamic_cast<Vector_MrHyDE_State<RealT>&>(u_out);
    for (int i=0; i<eu_in.mrhyde_state_vec.size(); ++i){
      prior_fe_op_->M[i]->apply(*eu_in.mrhyde_state_vec[i][0],*eu_out.mrhyde_state_vec[i][0]);
    }
  }
};

#endif

