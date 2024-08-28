#ifndef HDSA_RANDOMIZED_GEVP_HPP
#define HDSA_RANDOMIZED_GEVP_HPP

// This class executes the randomized GEVP solver

namespace HDSA
{

  template <class RealT>
  class Randomized_GEVP
  {
  private:
    HDSA::Ptr<HDSA::Vector<RealT> > vec_;

  public:

    Randomized_GEVP(const HDSA::Vector<RealT> & vec)
    {
      vec_ = vec.clone();
    }

    virtual ~Randomized_GEVP()
    { }

    virtual void Apply_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const = 0;

    virtual void Apply_Weighting_Operator(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const = 0;

    virtual void Apply_Weighting_Operator_Inverse(HDSA::Vector<RealT> & vec_out, const HDSA::Vector<RealT> & vec_in) const = 0;

    virtual void Generate_Random_Samples(HDSA::MultiVector<RealT> & samples) const = 0;
  
    void Compute_GEVP(HDSA::MultiVector<RealT> & evecs, HDSA::Dense_Matrix<RealT> & evals, const int & num_evals, const int & oversampling)
    {
      int kpp = num_evals + oversampling;

      HDSA::Ptr<HDSA::MultiVector<RealT> > Y = HDSA::makePtr<HDSA::MultiVector<RealT> >(kpp,*vec_);   
      HDSA::Ptr<HDSA::MultiVector<RealT> > tmp = HDSA::makePtr<HDSA::MultiVector<RealT> >(kpp,*vec_);
      Generate_Random_Samples(*tmp);
      for(int k = 0; k < kpp; k++)
      	{
	  HDSA::Ptr<HDSA::Vector<RealT> >  vec_tmp1 = vec_->clone();
	  Apply_Operator(*vec_tmp1,*(*tmp)[k]);
	  Apply_Weighting_Operator_Inverse(*(*Y)[k],*vec_tmp1);
      	}

      HDSA::Ptr<HDSA::MultiVector<RealT> > Q = HDSA::makePtr<HDSA::MultiVector<RealT> >(kpp,*vec_);
      HDSA::Ptr<HDSA::MultiVector<RealT> > WQ = HDSA::makePtr<HDSA::MultiVector<RealT> >(kpp,*vec_);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
      std::string type = "weighting";
      CholQR(*Q, *WQ, *R, *Y, type);

      Y->zeros();
      for(int k = 0; k < kpp; k++)
	{
	  Apply_Operator(*(*Y)[k],*(*Q)[k]);
	}
      
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > T = Y->MatMat(*Q);

      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_T = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
      int info = HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(*T,*R_T);

      if(info == 0)
	{
	  HDSA::Ptr<HDSA::MultiVector<RealT> > M = HDSA::makePtr<HDSA::MultiVector<RealT> >(kpp,*vec_);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > I = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
	  for(int k = 0; k < kpp; k++)
	    {
	      I->Replace_Element(k,k,1.0);
	    }
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_T_inv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
	  HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*R_T_inv,*I,*R_T);
	  for(int k = 0; k < kpp; k++)
	    {
	      for(int i = 0; i < kpp; i++)
		{
		  (*M)[k]->axpy((*R_T_inv)(i,k),*(*Y)[i]);
		}
	    }
	  
	  type = "weighting_inverse";
	  Q->zeros();
	  WQ->zeros();
	  R->zeros();
	  CholQR(*Q, *WQ, *R, *M, type);
	  
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > U = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > VT = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,1);
	  HDSA::Linear_Algebra::SVD<RealT>(*R, *U, *VT, *S);
	  
	  evecs.zeros();
	  for(int k = 0; k < num_evals; k++)
	    {
	      evals.Replace_Element(k,0,std::pow((*S)(k,0),2.0));
	      for(int i = 0; i < kpp; i++)
		{
		  evecs[k]->axpy((*U)(i,k),*(*WQ)[i]);
		}
	      if( (*U)(0,k) < 0.0 )
		{
		  evecs[k]->scale(-1.0);
		}
	    }
	}
      else
	{
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > S = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,kpp);
	  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > Lambda = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(kpp,1);
	  HDSA::Linear_Algebra::Symmetric_Eig_Decomposition<RealT>(*T, *S, *Lambda); 

	  evecs.zeros();
	  for(int k = 0; k < num_evals; k++)
	    {
	      evals.Replace_Element(k,0,(*Lambda)(k,0));
	      for(int i = 0; i < kpp; i++)
		{
		  evecs[k]->axpy((*S)(i,k),*(*Q)[i]);
		}
	      if( (*S)(0,k) < 0.0 )
		{
		  evecs[k]->scale(-1.0);
		}
	    }
	}
      
    }

    void CholQR(HDSA::MultiVector<RealT> & Q, HDSA::MultiVector<RealT> & WQ, HDSA::Dense_Matrix<RealT> & R, const HDSA::MultiVector<RealT> & Z, std::string & type)
    {   
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > ZtZ = Z.MatMat(Z);
      int n = ZtZ->numCols();
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_Z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(*ZtZ,*R_Z);
      
      HDSA::Ptr<HDSA::MultiVector<RealT> > Q_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(n,*Z[0]);
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > I = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
      for(int k = 0; k < n; k++)
	{
	  I->Replace_Element(k,k,1.0);
	}
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_Z_inv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
      HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*R_Z_inv,*I,*R_Z);
      for(int k = 0; k < n; k++)
	{
	  for(int i = 0; i < n; i++)
	    {
	      (*Q_Z)[k]->axpy((*R_Z_inv)(i,k),*Z[i]);
	    }
	}
      // Results in the factorization Z = Q_Z*R_Z

      // Applying weighting matrix to Q_Z
      HDSA::Ptr<HDSA::MultiVector<RealT> > W_Q_Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(n,*Z[0]);
      for(int k = 0; k < n; k++)
	{
	  if(type=="weighting")
	    {
	      Apply_Weighting_Operator(*(*W_Q_Z)[k],*(*Q_Z)[k]);
	    }
	  else if(type=="weighting_inverse")
	    {
              Apply_Weighting_Operator_Inverse(*(*W_Q_Z)[k],*(*Q_Z)[k]);
	    }
	  else
	    {
	      std::cout << "Error in type specification for CholQR" << std::endl;
	    }
	}

      // Compute C = Z^T*W*Z = Z^T*X
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > C = Q_Z->MatMat(*W_Q_Z);

      // Compute R_C=chol(C)
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_C = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
      HDSA::Linear_Algebra::Cholesky_Factorization<RealT>(*C,*R_C);

      // Compute Q=Q_Z*R_C^{-1}
      HDSA::Ptr<HDSA::Dense_Matrix<RealT> > R_C_inv = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n,n);
      HDSA::Linear_Algebra::Upper_Tri_Solve<RealT>(*R_C_inv,*I,*R_C);
    
      for(int k = 0; k < n; k++)
	{
	  for(int i = 0; i < n; i++)
	    {
	      Q[k]->axpy((*R_C_inv)(i,k),*(*Q_Z)[i]);
	      WQ[k]->axpy((*R_C_inv)(i,k),*(*W_Q_Z)[i]);
	    }
	}

      R_C->Multiply(R,*R_Z);

    }

  };
    
}
  


#endif
