#ifndef HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP
#define HDSA_MD_ELLIPTIC_Z_PRIOR_INTERFACE_MRHYDE_HPP

template <class RealT,
	  class LO=Tpetra::Map<>::local_ordinal_type,
	  class GO=Tpetra::Map<>::global_ordinal_type,
	  class Node=Tpetra::Map<>::node_type >
class MD_Elliptic_z_Prior_Interface_MrHyDE : public HDSA::MD_Elliptic_z_Prior_Interface<RealT> {

private:

  HDSA::Ptr<Prior_FE_Op_MrHyDE<RealT> > prior_fe_op_;
  Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LO,GO,Node> > E_z_;
  HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT> >  E_z_solver_;
  HDSA::Ptr<HDSA::Sparse_Matrix_Solver<RealT> >  M_z_solver_;
    
public:
  MD_Elliptic_z_Prior_Interface_MrHyDE(RealT & alpha_z, RealT & beta_z, HDSA::Ptr<Prior_FE_Op_MrHyDE<RealT>> &prior_fe_op):
    HDSA::MD_Elliptic_z_Prior_Interface<RealT>(alpha_z), prior_fe_op_(prior_fe_op)
  {
    Construct_Elliptic_Operator(beta_z);
  }
  
  virtual ~MD_Elliptic_z_Prior_Interface_MrHyDE()
  { }

  void Construct_Elliptic_Operator(RealT &beta_z) {
    Tpetra::RowMatrix<RealT, LO, GO, Node> &S_tmp = dynamic_cast<Tpetra::RowMatrix<RealT, LO, GO, Node> &> (*prior_fe_op_->S[0]);
    Teuchos::RCP<Tpetra::RowMatrix<RealT, LO, GO, Node> > E_z_tmp = prior_fe_op_->M[0]->add(beta_z,S_tmp,1.0,prior_fe_op_->M[0]->getDomainMap(),prior_fe_op_->M[0]->getRangeMap(),Teuchos::null);
    E_z_ = HDSA::dynamicPtrCast<Tpetra::CrsMatrix<RealT, LO, GO, Node> > (E_z_tmp);
    E_z_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(E_z_);
    M_z_solver_ = HDSA::makePtr<HDSA::Sparse_Matrix_Solver<RealT>>(prior_fe_op_->M[0]);
  }
  
  void Apply_E_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
    const Opt_Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const Opt_Vector_MrHyDE<RealT>&>(z_in);  
    Opt_Vector_MrHyDE<RealT> &ez_out = dynamic_cast<Opt_Vector_MrHyDE<RealT>&>(z_out);
    const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
    MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);
    HDSA::Ptr<HDSA::Vector<RealT> > vec_in = HDSA::makePtr<Tpetra_Vector_MrHyDE<RealT> > (eez_in.getField()[0]->getVector(),ez_in.Get_random_number_generator());
    HDSA::Ptr<HDSA::Vector<RealT> > vec_out = HDSA::makePtr<Tpetra_Vector_MrHyDE<RealT> > (eez_out.getField()[0]->getVector(),ez_out.Get_random_number_generator()); 
    E_z_solver_->Apply_A_Inverse(*vec_out,*vec_in); 
  }
  
  void Apply_E_z_Inverse_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
    Apply_E_z_Inverse(z_out,z_in);
  } 

  void Apply_M_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const {
    const Opt_Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const Opt_Vector_MrHyDE<RealT>&>(z_in);  
    Opt_Vector_MrHyDE<RealT> &ez_out = dynamic_cast<Opt_Vector_MrHyDE<RealT>&>(z_out);
    const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
    MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);
    prior_fe_op_->M[0]->apply(*eez_in.getField()[0]->getVector(),*eez_out.getField()[0]->getVector());
  }

  // Compute samples from a mean zero Gaussian with covariance W_z^{-1}
  virtual void Sample_with_Covariance_W_z_Inverse(HDSA::MultiVector<RealT> & samples) const
  {
    samples.zeros();
  }   
    
  virtual void Apply_E_z(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const 
  {
    const Opt_Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const Opt_Vector_MrHyDE<RealT>&>(z_in);  
    Opt_Vector_MrHyDE<RealT> &ez_out = dynamic_cast<Opt_Vector_MrHyDE<RealT>&>(z_out);
    const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
    MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);
    E_z_->apply(*eez_in.getField()[0]->getVector(),*eez_out.getField()[0]->getVector());
  }
  
  virtual void Apply_E_z_Transpose(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    Apply_E_z(z_out,z_in);
  }

  virtual void Apply_M_z_Inverse(HDSA::Vector<RealT> & z_out, const HDSA::Vector<RealT> & z_in) const
  {
    const Opt_Vector_MrHyDE<RealT> &ez_in = dynamic_cast<const Opt_Vector_MrHyDE<RealT>&>(z_in);  
    Opt_Vector_MrHyDE<RealT> &ez_out = dynamic_cast<Opt_Vector_MrHyDE<RealT>&>(z_out);
    const MrHyDE_OptVector &eez_in = dynamic_cast<const MrHyDE_OptVector&>(*ez_in.mrhyde_vec);  
    MrHyDE_OptVector &eez_out = dynamic_cast<MrHyDE_OptVector&>(*ez_out.mrhyde_vec);
    HDSA::Ptr<HDSA::Vector<RealT> > vec_in = HDSA::makePtr<Tpetra_Vector_MrHyDE<RealT> > (eez_in.getField()[0]->getVector(),ez_in.Get_random_number_generator());
    HDSA::Ptr<HDSA::Vector<RealT> > vec_out = HDSA::makePtr<Tpetra_Vector_MrHyDE<RealT> > (eez_out.getField()[0]->getVector(),ez_out.Get_random_number_generator()); 
    M_z_solver_->Apply_A_Inverse(*vec_out,*vec_in);
  }

};

#endif

