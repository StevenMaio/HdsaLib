#ifndef HDSA_WEIGHT_MATRICES_HPP
#define HDSA_WEIGHT_MATRICES_HPP

#include <fstream>

// This class contains the weight matrices which the user should overload

namespace HDSA
{

  template <class RealT>
  class Weight_Matrices{

  private:
    HDSA::Ptr<HDSA::ParameterList> parlist_sensitivity_;

  public:

    Weight_Matrices(const HDSA::Ptr<HDSA::ParameterList> & parlist_sensitivity): parlist_sensitivity_(parlist_sensitivity)
    { }
    
    virtual ~Weight_Matrices()
    { }
    
    virtual HDSA::Ptr<HDSA::Weight_Matrices<RealT> > Construct_Weight_Matrices(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<const HDSA::Comm<int> > & comm) = 0;

    virtual void Apply_theta_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const = 0;
    
    virtual void Apply_z_Weight_Mat(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const = 0;
    
    // Precondition solve for theta_Weight_Mat inverse via solve P*A*P*y = P*b and x=P*y where we assume P=P^T, typically P is diagonal
    virtual void Apply_theta_Weight_Mat_Preconditioner(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const
    { 
      // This defaults to the identity, the user may redefine the function
      theta_out->set(*theta_in);
    }
    
    virtual void Apply_theta_Weight_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & theta_out, const HDSA::Ptr<HDSA::Vector<RealT> > & theta_in) const
    {    
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > weight_mat_op = HDSA::makePtr<Preconditioned_theta_Operator<RealT> >(this);
      RealT tol = parlist_sensitivity_->sublist("Parameter Weight Matrix Solve").get("Tolerance",1.e-5);
      HDSA::Ptr<HDSA::Vector<RealT> > theta_in_tmp = theta_in->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_out_tmp = theta_out->Clone();
      Apply_theta_Weight_Mat_Preconditioner(theta_in_tmp, theta_in);
      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(theta_out_tmp,theta_in_tmp,weight_mat_op,tol);
      Apply_theta_Weight_Mat_Preconditioner(theta_out, theta_out_tmp);
    }

    virtual void Apply_z_Weight_Mat_Inverse(HDSA::Ptr<HDSA::Vector<RealT> > & z_out, const HDSA::Ptr<HDSA::Vector<RealT> > & z_in) const 
    {
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > weight_mat_op = HDSA::makePtr<z_weight_mat_operator<RealT> >(this);
      RealT tol = 1.e-10;
      HDSA::Linear_Algebra::Iterative_Linear_Solve<RealT>(z_out,z_in,weight_mat_op,tol);
    }

    virtual void CholQR(HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Q,  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & Y, const std::string & type, 
			const HDSA::Ptr<HDSA::Opt_Problem_Objects<RealT> > & OP_Objects, const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & WQ = HDSA::nullPtr, 
			const HDSA::Ptr<HDSA::Dense_Matrix<RealT> > & R_Trans = HDSA::nullPtr) const
    {
      // type may equal "z" or "theta" or "theta inverse" or "joint"
      HDSA::Ptr<HDSA::Linear_Operator<RealT> > W;
      HDSA::Ptr<HDSA::Vector<RealT> > vec;
      bool W_inv = false;

      if(type == "z")
	{
	  W = HDSA::makePtr<z_weight_mat_operator<RealT> >(this);
	  vec = OP_Objects->z->Clone();
	}
      else if(type == "theta")
	{
	  W = HDSA::makePtr<theta_weight_mat_operator<RealT> >(this);
	  vec = OP_Objects->theta->Clone();
	}
      else if(type == "theta inverse")
	{
	  W_inv = true;
	  W = HDSA::makePtr<theta_weight_mat_operator<RealT> >(this);
	  vec = OP_Objects->theta->Clone();
	}
      else if(type == "joint")
	{
	  W = HDSA::makePtr<joint_weight_mat_operator<RealT> >(this);
	  vec = HDSA::makePtr<HDSA::Joint_Vector<RealT> >(OP_Objects->z,OP_Objects->theta);
	}

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
    
    void Construct_Weight_Matrix_Test(const HDSA::Ptr<HDSA::Vector<RealT> > & theta, const HDSA::Ptr<HDSA::Vector<RealT> > & z) const
    {
      int theta_dim = theta->dimension();
      int z_dim = z->dimension();
      HDSA::Ptr<HDSA::Vector<RealT> > z_in = z->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > z_out = z->Clone();
      std::vector<std::vector<RealT> > Mz(z_dim);
      for(int i = 0; i < z_dim; i++)
	{
	  Mz[i].resize(z_dim);
	}
      
      for(int j = 0; j < z_dim; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << z_dim << " for z weight matrix." << std::endl;
	  z_in->basis(j);
	  z_out->zero();
	  Apply_z_Weight_Mat(z_out,z_in);
	  for(int i = 0; i < z_dim; i++)
	    {
	      Mz[i][j] = (*z_out)(i);
	    }
	}
      
      // Write solutions to text files
      std::string name;
      std::ofstream fout;      
      name = "Mz.txt";
      fout.open(name);
      for(int i = 0; i < z_dim; i++)
	{
	  for(int j = 0; j < z_dim; j++)
	    {
	      fout << std::setprecision(16) << Mz[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();

      HDSA::Ptr<HDSA::Vector<RealT> > theta_in = theta->Clone();
      HDSA::Ptr<HDSA::Vector<RealT> > theta_out = theta->Clone();
      std::vector<std::vector<RealT> > Mp(theta_dim);
      for(int i = 0; i < theta_dim; i++)
	{
	  Mp[i].resize(theta_dim);
	}
      
      for(int j = 0; j < theta_dim; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << theta_dim << " for theta weight matrix." << std::endl;
	  theta_in->basis(j);
	  theta_out->zero();
	  Apply_theta_Weight_Mat(theta_out,theta_in);
	  for(int i = 0; i < theta_dim; i++)
	    {
	      Mp[i][j] = (*theta_out)(i);
	    }
	}
    
      // Write solutions to text files
      name = "Mp.txt";
      fout.open(name);
      for(int i = 0; i < theta_dim; i++)
	{
	  for(int j = 0; j < theta_dim; j++)
	    {
	      fout << std::setprecision(16) << Mp[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();

      std::vector<std::vector<RealT> > Mpinv(theta_dim);
      for(int i = 0; i < theta_dim; i++)
	{
	  Mpinv[i].resize(theta_dim);
	}
      
      for(int j = 0; j < theta_dim; j++)
	{
	  std::cout << "Computing column " << j+1 << " out of " << theta_dim << " for theta weight matrix inverse." << std::endl;
	  theta_in->basis(j);
	  theta_out->zero();
	  Apply_theta_Weight_Mat_Inverse(theta_out,theta_in);
	  for(int i = 0; i < theta_dim; i++)
	    {
	      Mpinv[i][j] = (*theta_out)(i);
	    }
	}
    
      // Write solutions to text files
      name = "Mp_inv.txt";
      fout.open(name);
      for(int i = 0; i < theta_dim; i++)
	{
	  for(int j = 0; j < theta_dim; j++)
	    {
	      fout << std::setprecision(16) << Mpinv[i][j] << "  ";
	    }
	  fout << "  " << std::endl;
	}
      fout.close();
          
    }
    
    // Overload Linear Operator to take matrix vector products for the parameter mass matrix inversion.
    template <class ScalarType>
    class Preconditioned_theta_Operator : public HDSA::Linear_Operator<ScalarType>
    {
      const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op_;
      
    public:
 
      Preconditioned_theta_Operator(const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op): Weight_Mat_Op_(Weight_Mat_Op)
      { }
     
      //! Dtor
      ~Preconditioned_theta_Operator()
      {}
      
   void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	HDSA::Ptr<HDSA::Vector<RealT> > temp1 = x->Clone();
	HDSA::Ptr<HDSA::Vector<RealT> > temp2 = x->Clone();
	Weight_Mat_Op_->Apply_theta_Weight_Mat_Preconditioner(temp1,x);
	Weight_Mat_Op_->Apply_theta_Weight_Mat(temp2,temp1);
	Weight_Mat_Op_->Apply_theta_Weight_Mat_Preconditioner(y,temp2);
      }
      
    };

    // Overload Linear Operator
    template <class ScalarType>
    class theta_weight_mat_operator : public HDSA::Linear_Operator<ScalarType>
    {
      const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op_;
      
    public:
 
      theta_weight_mat_operator(const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op): Weight_Mat_Op_(Weight_Mat_Op)
      { }
      
      //! Dtor
      ~theta_weight_mat_operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	Weight_Mat_Op_->Apply_theta_Weight_Mat(y,x);
      }
      
    };

    // Overload Linear Operator
    template <class ScalarType>
    class z_weight_mat_operator : public HDSA::Linear_Operator<ScalarType>
    {
      const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op_;
      
    public:
 
      z_weight_mat_operator(const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op): Weight_Mat_Op_(Weight_Mat_Op)
      { }
      
      //! Dtor
      ~z_weight_mat_operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	Weight_Mat_Op_->Apply_z_Weight_Mat(y,x);
      }
      
    };

    // Overload Linear Operator
    template <class ScalarType>
    class joint_weight_mat_operator : public HDSA::Linear_Operator<ScalarType>
    {
      const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op_;
      
    public:
 
      joint_weight_mat_operator(const HDSA::Weight_Matrices<ScalarType>* Weight_Mat_Op): Weight_Mat_Op_(Weight_Mat_Op)
      { }
      
      //! Dtor
      ~joint_weight_mat_operator()
      {}
      
      void matvec(HDSA::Ptr<HDSA::Vector<ScalarType> > & y, const HDSA::Ptr<HDSA::Vector<ScalarType> > & x) const 
      {
	Joint_Vector<RealT> &ex = dynamic_cast<Joint_Vector<RealT>&>(*x);
	Joint_Vector<RealT> &ey = dynamic_cast<Joint_Vector<RealT>&>(*y);
	HDSA::Ptr<HDSA::Vector<RealT> > x_z = ex.Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > x_theta = ex.Get_Component_Vector_2();
	HDSA::Ptr<HDSA::Vector<RealT> > y_z = ey.Get_Component_Vector_1();
	HDSA::Ptr<HDSA::Vector<RealT> > y_theta = ey.Get_Component_Vector_2();
	Weight_Mat_Op_->Apply_z_Weight_Mat(y_z,x_z);
	Weight_Mat_Op_->Apply_theta_Weight_Mat(y_theta,x_theta);
      }
      
    };

  };

}

#endif
