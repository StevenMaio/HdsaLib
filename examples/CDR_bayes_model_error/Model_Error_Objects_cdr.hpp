#ifndef BAYES_MODEL_ERROR_OBJECTS_CDR_HPP
#define BAYES_MODEL_ERROR_OBJECTS_CDR_HPP

// Instantiation of Bayes_Model_Error_Objects

template <class RealT>
class Bayes_Model_Error_Objects_CDR : public HDSA::Bayes_Model_Error_Objects<RealT> {

private:
  HDSA::Ptr<HDSA::ParameterList> parlist_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_AG_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_AL_;
  HDSA::Ptr<ROL::Constraint_SimOpt<RealT> > con_M_;
  HDSA::Ptr<HDSA::Vector<RealT> > M_row_sum_;
  HDSA::Ptr<HDSA::Vector<RealT> > Ainv_sing_vals_;
  HDSA::Ptr<HDSA::MultiVector<RealT> > Ainv_U_;
  HDSA::Ptr<HDSA::MultiVector<RealT> > Ainv_V_;

public:

  Bayes_Model_Error_Objects_CDR(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity, const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects_Factory, 
				const HDSA::Ptr<HDSA::Weight_Matrices<RealT> > & weight_matrices_factory, const HDSA::Ptr<HDSA::ParameterList> & parlist):
    HDSA::Bayes_Model_Error_Objects<RealT>(parlist_sensitivity,OP_Objects_Factory,weight_matrices_factory)
  {
    parlist_ = parlist;
  }

  virtual ~Bayes_Model_Error_Objects_CDR()
  { }

  void Construct_Objects(const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    HDSA::Ptr<std::ostream> outStream;
    HDSA::nullstream bhs;
    int myRank = comm->getRank();
    if(myRank == 0)
      {
	outStream = HDSA::makePtrFromRef(std::cout);
      }
    else
      {	
	outStream =  HDSA::makePtrFromRef(bhs);
      }
    HDSA::Ptr<MeshManager<RealT> > meshMgr = HDSA::makePtr<MeshManager_CDR<RealT> >(*parlist_);

    RealT Gamma_diff = parlist_->sublist("Problem").get("Prior Control Diff Coefficient",1.e-5);
    RealT Gamma_eye = parlist_->sublist("Problem").get("Prior Control Identity Coefficient",1.0);
    HDSA::Ptr<PDE<RealT> > pde_AG = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist_,Gamma_diff,Gamma_eye);
    con_AG_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_AG,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_AG_->setSolveParameters(*parlist_);

    RealT L_diff = parlist_->sublist("Problem").get("Prior Discrepancy Diff Coefficient",1.0);
    RealT L_eye = parlist_->sublist("Problem").get("Prior Discrepancy Identity Coefficient",1.0);
    HDSA::Ptr<PDE<RealT> > pde_AL = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist_,L_diff,L_eye);
    con_AL_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_AL,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_AL_->setSolveParameters(*parlist_);

    HDSA::Ptr<PDE<RealT> > pde_M = HDSA::makePtr<Elliptic_Op<RealT> >(*parlist_,0.0,1.0);
    con_M_ = HDSA::makePtr<Linear_PDE_Constraint<RealT> >(pde_M,meshMgr,comm->Get_Teuchos_Communicator(),*parlist_,*outStream);
    con_M_->setSolveParameters(*parlist_);

    bool checkDeriv = false;
    if ( checkDeriv ) {
      HDSA::Ptr<HDSA::Vector<RealT> > up_hdsa = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > zp_hdsa = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->Clone();
      HDSA::Ptr<ROL::Vector<RealT> > up = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(up_hdsa)->get_rol_vec();
      HDSA::Ptr<ROL::Vector<RealT> > zp = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(zp_hdsa)->get_rol_vec();
      HDSA::Ptr<ROL::Vector<RealT> > dup = up->clone(); dup->randomize(-1.0,1.0);
      HDSA::Ptr<ROL::Vector<RealT> > dzp = zp->clone(); dzp->randomize(-1.0,1.0);
      con_AG_->checkApplyJacobian_1(*up,*zp,*dup,*up,true,*outStream);
      con_AG_->checkInverseJacobian_1(*up,*up,*up,*zp,true,*outStream);
    }

    M_row_sum_ = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = M_row_sum_->Clone();
    RealT tol = 1.e-8;
    HDSA::Ptr<ROL::Vector<RealT> > v1 = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(M_row_sum_)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > v2 = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(tmp)->get_rol_vec();
    v2->setScalar(1.0);
    con_M_->applyJacobian_1(*v1,*v2,*v2,*v2,tol);

    int k = parlist_->sublist("Problem").get("Prior GSVD Rank",400);
    int q = parlist_->sublist("Problem").get("Prior GSVD Subspace Iterations",1);
    A_L_Decomposition(k,q,comm);

    // Define Gamma = A_G^{-1}*M*A_G^{-1}
    // Define L = A_L*M^{-1}*A_L
  }
  
  void A_L_Decomposition(int k, int q, const HDSA::Ptr<const HDSA::Comm<int> > & comm)
  {
    std::clock_t timer = std::clock();

    // Initial sketch
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(k,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    int m = (*Y)[0]->dimension();
    for(int i = 0; i < k; i++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > ui = (*Y)[i];
	HDSA::Ptr<HDSA::Vector<RealT> > u_vec = ui->Generate_Gaussian_Random_Vector();
	Apply_A_L_Mat_Inv(ui,u_vec);	
      }
    HDSA::Ptr<HDSA::MultiVector<RealT> > MQ = HDSA::makePtr<HDSA::MultiVector<RealT> >(k,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    HDSA::Ptr<HDSA::MultiVector<RealT> > Q = HDSA::makePtr<HDSA::MultiVector<RealT> >(k,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    CholQR_M(MQ,Q,Y);
    
    // Subspace iteration for improving sketch
    for(int j = 0; j < q; j++)
      {
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
	for(int i = 0; i < k; i++)
	  {
	    HDSA::Ptr<HDSA::Vector<RealT> > ui = (*Y)[i];
	    ui->zero();
	    HDSA::Ptr<HDSA::Vector<RealT> > MQi = (*MQ)[i];
	    Apply_A_L_Mat_Inv(ui,MQi);	
	    Y_mat->Write_Vector_to_Column(i,ui);
	  }
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);	
	HDSA::Linear_Algebra::QR_Factorization<RealT>(Y_mat,Q_Y);
	for(int i = 0; i < k; i++)
	  {
	    HDSA::Ptr<HDSA::Vector<RealT> > ui = (*Y)[i];
	    ui->zero();
	    HDSA::Ptr<HDSA::Vector<RealT> > MQi = (*MQ)[i];
	    Q_Y->Write_Column_to_Vector(i,MQi);
	    Apply_A_L_Mat_Inv(ui,MQi);	
	  }
	CholQR_M(MQ,Q,Y);
      }

    // Projection on sketched subspace
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    for(int i = 0; i < k; i++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > ui = (*Y)[i];
	ui->zero();
	HDSA::Ptr<HDSA::Vector<RealT> > MQi = (*MQ)[i];
	Apply_A_L_Mat_Inv(ui,MQi);
	B->Write_Vector_to_Column(i,ui);	
      }
    // QR of projected matrix
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    HDSA::Linear_Algebra::QR_Factorization<RealT>(B,Q_B,R_B);

    // SVD of R_B^T
    Ainv_sing_vals_ = HDSA::makePtr<Std_Vector<RealT> >(k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U_Bt = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V_B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    HDSA::Linear_Algebra::SVD(R_B, V_B,  U_Bt, Ainv_sing_vals_);

    // Mapping to high dimensional vectors
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    for(int i = 0; i < k; i++)
      {
	Q_mat->Write_Vector_to_Column(i,(*Q)[i]);
      }
    Q_mat->Multiply(U,U_Bt,false,true);
    Q_B->Multiply(V,V_B);
    Ainv_U_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(k,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    Ainv_V_ = HDSA::makePtr<HDSA::MultiVector<RealT> >(k,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    for(int i = 0; i < k; i++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > vi = (*Ainv_V_)[i];
	V->Write_Column_to_Vector(i,vi);
	HDSA::Ptr<HDSA::Vector<RealT> > ui = (*Ainv_U_)[i];
	U->Write_Column_to_Vector(i,ui);
      }

    RealT Time = static_cast<RealT>(std::clock()-timer)/static_cast<RealT>(CLOCKS_PER_SEC);
    std::cout << "GSVD took " << Time << " seconds to execute" << std::endl;

    std::cout << "First singular value = " << (*Ainv_sing_vals_)(0) << " and last singular value = " << (*Ainv_sing_vals_)(k-1) << std::endl;
  }

  void CholQR_M(HDSA::Ptr<HDSA::MultiVector<RealT> > & MQ, HDSA::Ptr<HDSA::MultiVector<RealT> > & Q, HDSA::Ptr<HDSA::MultiVector<RealT> > & Z)
  {
    int m = (*Z)[0]->dimension();
    int k = Z->Number_of_Vectors();
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Z_mat = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k); 
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m, k);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    for(int i = 0; i < k; i++)
      {
	Z_mat->Write_Vector_to_Column(i,(*Z)[i]);
      }
    HDSA::Linear_Algebra::QR_Factorization<RealT>(Z_mat,Q_Z,R_Z);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    for(int j = 0; j < k; j++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Std_Vector<RealT> >(m);
	for(int l = 0; l < m; l++)
	  {
	    vec->Replace_Element(l,(*Q_Z)(l,j)*(*M_row_sum_)(l));
	  }
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > val = Q_Z->Multiply(vec, true);
	for(int i = 0; i < k; i++)
	  {
	    C->Replace_Element(i,j,(*val)(i,0));
	  }
      }
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_W = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(k, k);
    HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(C,R_W);
    // MQ = M*Q_Z*R_W^{-1}
    for(int i = 0; i < k; i++)
      {
	HDSA::Ptr<HDSA::Vector<RealT> > x = HDSA::makePtr<Std_Vector<RealT> >(k);
	HDSA::Ptr<HDSA::Vector<RealT> > b = HDSA::makePtr<Std_Vector<RealT> >(k);
	b->basis(i);
	HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(x,b,R_W);
	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > val = Q_Z->Multiply(x);
	HDSA::Ptr<HDSA::Vector<RealT> > Q_vec = (*Q)[i];
	HDSA::Ptr<HDSA::Vector<RealT> > MQ_vec = (*MQ)[i];
	for(int l = 0; l < m; l++)
	  {
	    Q_vec->Replace_Element(l,(*val)(l,0));
	    MQ_vec->Replace_Element(l,(*val)(l,0)*(*M_row_sum_)(l));
	  }
      }
  }


  void Apply_A_G_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AG_->applyJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
    z_out->Set_Zeros();
  }

  void Apply_A_G_Mat_Inv(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > z_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > z_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(z_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AG_->applyInverseJacobian_1(*z_out_rol,*z_in_rol,*z_in_rol,*z_in_rol,tol);
    z_out->Set_Zeros();
  }

  void Apply_A_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AL_->applyJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_A_L_Mat_Inv(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const
  {
    HDSA::Ptr<ROL::Vector<RealT> > u_in_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_in)->get_rol_vec();
    HDSA::Ptr<ROL::Vector<RealT> > u_out_rol = HDSA::dynamicPtrCast<ROL_Vector<RealT> >(u_out)->get_rol_vec();
    RealT tol = 1.e-8;
    con_AL_->applyInverseJacobian_1(*u_out_rol,*u_in_rol,*u_in_rol,*u_in_rol,tol);
  }

  void Apply_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_in->Clone();
    Apply_A_G_Mat_Inv(tmp,z_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)*(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_G_Mat_Inv(z_out,tmp); 
  }

  void Apply_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = z_in->Clone();
    Apply_A_G_Mat(tmp,z_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)/(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_G_Mat(z_out,tmp); 
  }

  void Apply_L_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > tmp = u_in->Clone();
    Apply_A_L_Mat(tmp,u_in);
    for(int k = 0; k < tmp->dimension(); k++)
      {
	RealT val = (*tmp)(k)/(*M_row_sum_)(k);
	tmp->Replace_Element(k,val);
      }
    Apply_A_L_Mat(u_out,tmp); 
  }

  void Apply_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > w = Ainv_V_->MatVec(*u_in);
    int k = Ainv_sing_vals_->dimension();
    RealT val = 0.0;
    u_out->zero();
    for(int i = 0; i < k; i++)
      {
        val = (*w)(i)*std::pow((*Ainv_sing_vals_)(i),2);
	u_out->axpy(val,*(*Ainv_V_)[i]);
      }

    // HDSA::Ptr<HDSA::Vector<RealT> > tmp = u_in->Clone();
    // Apply_A_L_Mat_Inv(tmp,u_in);
    // for(int k = 0; k < tmp->dimension(); k++)
    //   {
    // 	RealT val = (*tmp)(k)*(*M_row_sum_)(k);
    // 	tmp->Replace_Element(k,val);
    //   }
    // Apply_A_L_Mat_Inv(u_out,tmp); 
  }

  void Apply_Sqrt_Gamma_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    for(int k = 0; k < z_in->dimension(); k++)
      {
	RealT val = (*z_in)(k)*std::sqrt((*M_row_sum_)(k));
	z_in->Replace_Element(k,val);
      }
    Apply_A_G_Mat_Inv(z_out,z_in);
  }

  void Apply_Sqrt_Gamma_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    for(int k = 0; k < z_in->dimension(); k++)
      {
	RealT val = (*z_in)(k)/std::sqrt((*M_row_sum_)(k));
	z_in->Replace_Element(k,val);
      }
    Apply_A_G_Mat(z_out,z_in);
  }
  
  void Apply_Sqrt_L_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > w = Ainv_V_->MatVec(*u_in);
    int k = Ainv_sing_vals_->dimension();
    RealT val = 0.0;
    u_out->zero();
    for(int i = 0; i < k; i++)
      {
        val = (*w)(i)*(*Ainv_sing_vals_)(i);
	u_out->axpy(val,*(*Ainv_V_)[i]);
      }

    // for(int k = 0; k < u_in->dimension(); k++)
    //   {
    // 	RealT val = (*u_in)(k)*std::sqrt((*M_row_sum_)(k));
    // 	u_in->Replace_Element(k,val);
    //   }
    // Apply_A_L_Mat_Inv(u_out,u_in);
  }
  
  void Apply_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > w = Ainv_V_->MatVec(*u_in);
    int k = Ainv_sing_vals_->dimension();
    RealT val = 0.0;
    u_out->zero();
    for(int i = 0; i < k; i++)
      {
        val = (*w)(i)*(std::pow((*Ainv_sing_vals_)(i),2)/(1.0+beta*std::pow((*Ainv_sing_vals_)(i),2)));
	u_out->axpy(val,*(*Ainv_V_)[i]);
      }
  }

  void Apply_Sqrt_L_Plus_Shift_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & u_out, const HDSA::Ptr<HDSA::Vector<RealT> > & u_in, const RealT beta) const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > w = Ainv_V_->MatVec(*u_in);
    int k = Ainv_sing_vals_->dimension();
    RealT val = 0.0;
    u_out->zero();
    for(int i = 0; i < k; i++)
      {
        val = (*w)(i)*((*Ainv_sing_vals_)(i)/std::sqrt(1.0+beta*std::pow((*Ainv_sing_vals_)(i),2)));
	u_out->axpy(val,*(*Ainv_V_)[i]);
      }
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->z->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->z);
    // read in data
    std::ifstream in_Z("Z.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Z) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Z >> val;
	      (*Z)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Z.txt" << std::endl;
      } 
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Y_Data(void) const 
  {
    int num_vecs = 1;
    int dim = HDSA::Model_Error_Objects<RealT>::OP_Objects_->u->dimension();
    HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(num_vecs,HDSA::Model_Error_Objects<RealT>::OP_Objects_->u);
    // read in data
    std::ifstream in_Y("Y.txt");           
    // read the elements in the file into a vector  
    // test file open   
    RealT val;
    if (in_Y) {   
      for(int i = 0; i < dim; i++)
	{
	  for(int j = 0; j < num_vecs; j++)
	    {
	      in_Y >> val;
	      (*Y)[j]->Replace_Element(i,val);
	    }
	}   
    }
    else
      {
	std::cout << "Error loading the data from Y.txt" << std::endl;
      } 
    return Y;
  }

};


#endif
