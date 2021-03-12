#ifndef HDSA_WEIGHT_MATRICES_MODEL_ERROR_HPP
#define HDSA_WEIGHT_MATRICES_MODEL_ERROR_HPP

namespace HDSA
{

template <class RealT>
class Weight_Matrices_Model_Error : public Weight_Matrices<RealT>{

  HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > model_error_objects_;
  int seed_;

public:

  Weight_Matrices_Model_Error(HDSA::Ptr<HDSA::Model_Error_Objects<RealT> > & model_error_objects):
    HDSA::Weight_Matrices<RealT>(model_error_objects->parlist_sensitivity_), model_error_objects_(model_error_objects)
  { 
    bool time_seed = model_error_objects_->parlist_sensitivity_->sublist("Formulation").get("System Time Seed",false);
    if(time_seed)
      {
	seed_ = time(NULL);
      }
    else
      {
	seed_ = 2435;
      }

  }

  virtual ~Weight_Matrices_Model_Error()
  { }

  HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) 
  { 
    model_error_objects_->Instantiate_Objects(theta,comm);
    HDSA::Ptr<HDSA::Weight_Matrices<RealT> > weight_matrices = HDSA::makePtr<Weight_Matrices_Model_Error<RealT> >(model_error_objects_);
    return weight_matrices;
  }

  void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const
  {
    theta_out->zero();
    
    std::default_random_engine generator;
    generator.seed(seed_);
    std::normal_distribution<RealT> distribution = std::normal_distribution<RealT>(0.0,1.0);
    
    const Vector_Model_Error<RealT> &etheta_in = dynamic_cast<const Vector_Model_Error<RealT>&>(*theta_in);
    Vector_Model_Error<RealT> &etheta_out = dynamic_cast<Vector_Model_Error<RealT>&>(*theta_out);
    
    // Compute Y = A*Omega for range finding
    int kpp = 2;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Omega;
    Omega.resize(kpp);
    for(int k = 0; k < kpp; k++)
      {
	Omega[k] = model_error_objects_->OP_Objects_->z->Clone();
	for(int j = 0; j < model_error_objects_->n_; j++)
	  {
	    Omega[k]->Replace_Element(j,distribution(generator));
	  }
      }
    
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Theta_Omega;
    Theta_Omega.resize(kpp);
    for(int k = 0; k < kpp; k++)
      {
	Theta_Omega[k] = model_error_objects_->OP_Objects_->u->Clone();
	for(int i = 0; i < model_error_objects_->m_; i++)
	  {
	    Theta_Omega[k]->Replace_Element(i,etheta_in.Get_Row(i)->dot(*Omega[k]));
	  }
      }
    
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    for(int k = 0; k < kpp; k++)
      {
	Y->Write_Vector_to_Column(k,Theta_Omega[k]);
      }
    
    // Orthogonalize columns of Q
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    HDSA::Linear_Algebra::QR_Factorization<RealT>(Y, Q);
    
    // B = Theta^T*Q
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->n_,kpp);
    std::vector<std::vector<RealT> > Theta_std;
    Theta_std.resize(model_error_objects_->m_);
    std::vector<int> Theta_indices = etheta_in.Get_Row(0)->Get_Indices_on_Processor();
    for(int i = 0; i < model_error_objects_->m_; i++)
      {
    	Theta_std[i] = etheta_in.Get_Row(i)->Get_Data_on_Processor();
      }
    for(unsigned int j = 0; j < Theta_indices.size(); j++)
      {
    	for(int k = 0; k < kpp; k++)
    	  {
    	    RealT val = 0.0;
    	    for(int i = 0; i < model_error_objects_->m_; i++)
    	      {
    		val += Theta_std[i][Theta_indices[j]]*(*Q)(i,k);
    	      }
    	    B->Replace_Element(Theta_indices[j],k,val);
    	  }
      }
    
    // SVD of B
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->n_,kpp);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(kpp);
    HDSA::Linear_Algebra::SVD<RealT>(B, U, VT, S);
    
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    Q->Multiply(Q_V,VT,false,true);
    
    HDSA::Ptr<HDSA::Vector<RealT> > L_Q_V_0 = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > L_Q_V_1 = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > E_U_0 = model_error_objects_->OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > E_U_1 = model_error_objects_->OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_vec = model_error_objects_->OP_Objects_->z->Clone();
    
    Q_V->Write_Column_to_Vector(0,u_vec);
    model_error_objects_->Apply_L_Mat(L_Q_V_0,u_vec);
    Q_V->Write_Column_to_Vector(1,u_vec);
    model_error_objects_->Apply_L_Mat(L_Q_V_1,u_vec);
    
    U->Write_Column_to_Vector(0,z_vec);
    Apply_E(E_U_0,z_vec);
    U->Write_Column_to_Vector(1,z_vec);
    Apply_E(E_U_1,z_vec);
    
    for(int i = 0; i < model_error_objects_->m_; i++)
      {
	z_vec->set(*E_U_0);
	z_vec->scale((*S)(0)*(*L_Q_V_0)(i));
	z_vec->axpy((*S)(1)*(*L_Q_V_1)(i),*E_U_1);
	etheta_out.Set_Row(*z_vec,i);
      }
    RealT coeff = std::pow(model_error_objects_->alpha_,2.0)*model_error_objects_->nom_state_inner_prod_;
    theta_out->scale(1.0/coeff);
  }
    
  void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    model_error_objects_->weight_matrices_->Apply_z_Weight_Mat(z_out,z_in);
  }

  virtual void Apply_theta_Weight_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const
  {
    theta_out->zero();
	
    std::default_random_engine generator;
    generator.seed(seed_);
    std::normal_distribution<RealT> distribution = std::normal_distribution<RealT>(0.0,1.0);
    
    const Vector_Model_Error<RealT> &etheta_in = dynamic_cast<const Vector_Model_Error<RealT>&>(*theta_in);
    Vector_Model_Error<RealT> &etheta_out = dynamic_cast<Vector_Model_Error<RealT>&>(*theta_out);
    
    // Compute Y = A*Omega for range finding
    int kpp = 2;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Omega;
    Omega.resize(kpp);
    for(int k = 0; k < kpp; k++)
      {
	Omega[k] = model_error_objects_->OP_Objects_->z->Clone();
	for(int j = 0; j < model_error_objects_->n_; j++)
	  {
	    Omega[k]->Replace_Element(j,distribution(generator));
	  }
      }

    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > Theta_Omega;
    Theta_Omega.resize(kpp);
    for(int k = 0; k < kpp; k++)
      {
	Theta_Omega[k] = model_error_objects_->OP_Objects_->u->Clone();
	for(int i = 0; i < model_error_objects_->m_; i++)
	  {
	    Theta_Omega[k]->Replace_Element(i,etheta_in.Get_Row(i)->dot(*Omega[k]));
	  }
      }

    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Y = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    for(int k = 0; k < kpp; k++)
      {
	Y->Write_Vector_to_Column(k,Theta_Omega[k]);
      }

    // Orthogonalize columns of Q
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    HDSA::Linear_Algebra::QR_Factorization<RealT>(Y, Q);   
    
    // B = Theta^T*Q
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > B = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->n_,kpp);
    std::vector<std::vector<RealT> > Theta_std;
    Theta_std.resize(model_error_objects_->m_);
    std::vector<int> Theta_indices = etheta_in.Get_Row(0)->Get_Indices_on_Processor();
    for(int i = 0; i < model_error_objects_->m_; i++)
      {
    	Theta_std[i] = etheta_in.Get_Row(i)->Get_Data_on_Processor();
      }
    for(unsigned int j = 0; j < Theta_indices.size(); j++)
      {
    	for(int k = 0; k < kpp; k++)
    	  {
    	    RealT val = 0.0;
    	    for(int i = 0; i < model_error_objects_->m_; i++)
    	      {
    		val += Theta_std[i][Theta_indices[j]]*(*Q)(i,k);
    	      }
    	    B->Replace_Element(Theta_indices[j],k,val);
    	  }
      }

    // SVD of B
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->n_,kpp);
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
    HDSA::Ptr<HDSA::Vector<RealT> > S = HDSA::makePtr<Std_Vector<RealT> >(kpp);
    HDSA::Linear_Algebra::SVD<RealT>(B, U, VT, S);
    
    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Q_V = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(model_error_objects_->m_,kpp);
    Q->Multiply(Q_V,VT,false,true);
    
    HDSA::Ptr<HDSA::Vector<RealT> > Linv_Q_V_0 = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > Linv_Q_V_1 = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > u_vec = model_error_objects_->OP_Objects_->u->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > Einv_U_0 = model_error_objects_->OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > Einv_U_1 = model_error_objects_->OP_Objects_->z->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_vec = model_error_objects_->OP_Objects_->z->Clone();
    
    HDSA::Ptr<HDSA::Vector<RealT> > gamma_inv_z_star = model_error_objects_->OP_Objects_->z->Clone();
    for(int k = 0; k < model_error_objects_->n_; k++)
      {
	RealT val = model_error_objects_->OP_Objects_->z->Get_Element(k)/model_error_objects_->z_cov_[k];
	gamma_inv_z_star->Replace_Element(k,val);
      }
    RealT z_star_gamma_inv_z_star = gamma_inv_z_star->dot(*model_error_objects_->OP_Objects_->z);

    Q_V->Write_Column_to_Vector(0,u_vec);
    model_error_objects_->Apply_L_Mat_Inverse(Linv_Q_V_0,u_vec);
    Q_V->Write_Column_to_Vector(1,u_vec);
    model_error_objects_->Apply_L_Mat_Inverse(Linv_Q_V_1,u_vec);
    
    U->Write_Column_to_Vector(0,z_vec);
    Apply_E_Inverse(Einv_U_0,z_vec,gamma_inv_z_star,z_star_gamma_inv_z_star);
    U->Write_Column_to_Vector(1,z_vec);
    Apply_E_Inverse(Einv_U_1,z_vec,gamma_inv_z_star,z_star_gamma_inv_z_star);

    for(int i = 0; i < model_error_objects_->m_; i++)
      {
	z_vec->set(*Einv_U_0);
	z_vec->scale((*S)(0)*(*Linv_Q_V_0)(i));
	z_vec->axpy((*S)(1)*(*Linv_Q_V_1)(i),*Einv_U_1);
	etheta_out.Set_Row(*z_vec,i);
      }
    RealT coeff = std::pow(model_error_objects_->alpha_,2.0)*model_error_objects_->nom_state_inner_prod_;
    theta_out->scale(coeff);
  }

  virtual void Apply_z_Weight_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
  {
    model_error_objects_->weight_matrices_->Apply_z_Weight_Mat_Inverse(z_out,z_in);
  }

  void CholQR(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Q,  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Y, const std::string & type, 
	      const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & WQ = HDSA::nullPtr, 
	      const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R_Trans = HDSA::nullPtr) const
  {
    if(type == "theta inverse")
      {
	HDSA::Ptr<HDSA::Linear_Operator<RealT> > W = HDSA::makePtr<theta_Weight_Mat_Inverse_Operator<RealT> >(this);
	HDSA::Ptr<HDSA::Vector<RealT> > vec = OP_Objects->theta->Clone();
	bool W_inv = false;

	HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R;
	if(R_Trans == HDSA::nullPtr)
	  {
	    HDSA::Linear_Algebra::CholQR<RealT>(Y,Q,W,vec,WQ,R_Trans,W_inv);
	  }
	else
	  {
	    HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = R_Trans->Clone();
	    HDSA::Linear_Algebra::CholQR<RealT>(Y,Q,W,vec,WQ,R,W_inv);
	    for(int i = 0; i < R->numRows(); i++)
	      {
		for(int j = 0; j < R->numCols(); j++)
		  {
		    R_Trans->Replace_Element(i,j,(*R)(j,i));	    
		  }
	      }
	  }  
      }
    else
      {
	HDSA::Weight_Matrices<RealT>::CholQR(Q,Y,type,OP_Objects,WQ,R_Trans);
      }
  }

private:

  void Apply_E(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp = z_in->Clone();
    Apply_z_Weight_Mat(z_tmp,z_in);
    
    RealT val_coeff = z_tmp->dot(*model_error_objects_->OP_Objects_->z);
    for(int k = 0; k < model_error_objects_->n_; k++)
      {
	RealT val = (model_error_objects_->z_cov_[k])*(*z_tmp)(k);
	z_tmp->Replace_Element(k,val);
      }
    z_tmp->axpy(val_coeff,*model_error_objects_->OP_Objects_->z);

    Apply_z_Weight_Mat(z_out,z_tmp);
  }

  void Apply_E_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in, HDSA::Ptr<HDSA::Vector<RealT> > & gamma_inv_z_star, RealT z_star_gamma_inv_z_star) const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp1 = z_in->Clone();
    HDSA::Ptr<HDSA::Vector<RealT> > z_tmp2 = z_in->Clone();
    Apply_z_Weight_Mat_Inverse(z_tmp1,z_in);

    for(int k = 0; k < model_error_objects_->n_; k++)
      {
	RealT val = z_tmp1->Get_Element(k)/model_error_objects_->z_cov_[k];
	z_tmp2->Replace_Element(k,val);
      }
    RealT coeff = -z_tmp2->dot(*model_error_objects_->OP_Objects_->z)/(1.0+z_star_gamma_inv_z_star);
    z_tmp2->axpy(coeff,*gamma_inv_z_star);

    Apply_z_Weight_Mat_Inverse(z_out,z_tmp2);
  }

  // Overload HDSA::Linear_Operator to take matrix vector products
  template <class ScalarType>
  class theta_Weight_Mat_Inverse_Operator : public HDSA::Linear_Operator<ScalarType>
  {
    const HDSA::Weight_Matrices_Model_Error<ScalarType>* weight_matrices_model_error_;
    
  public:
    
    theta_Weight_Mat_Inverse_Operator(const HDSA::Weight_Matrices_Model_Error<ScalarType>* weight_matrices_model_error): weight_matrices_model_error_(weight_matrices_model_error)
    { }
    
    //! Dtor
    ~theta_Weight_Mat_Inverse_Operator()
    {}
    
    void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
    {
      weight_matrices_model_error_->Apply_theta_Weight_Mat_Inverse(y,x);
    }
    
  };
 

};

}

#endif


