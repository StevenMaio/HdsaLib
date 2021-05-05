#ifndef HDSA_MODEL_ERROR_KRONECKER_MATRIX_HPP
#define HDSA_MODEL_ERROR_KRONECKER_MATRIX_HPP

namespace HDSA
{

template <class RealT>
class Model_Error_Kronecker_Matrix {

private:
  int m_;
  int n_;
  int kpp_;
  HDSA::Ptr<HDSA::Processor_Distribution<RealT> > proc_dist_;
  
public:
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > uk; 	  
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > zk;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > u;	  
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > z;
  HDSA::Ptr<HDSA::Dense_Matrix<RealT> > bk; 	  
  RealT a;

  // Defines a set of kpp vectors of the form
  // [ a*u_k , kron(u_k,z) ] + [ b_k*u , kron(u,z_k) ] 

  Model_Error_Kronecker_Matrix(int m, int n, int kpp, HDSA::Ptr<HDSA::Processor_Distribution<RealT> > & proc_dist) 
  {
    m_ = m;
    n_ = n;
    kpp_ = kpp;
    proc_dist_ = proc_dist;

    uk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,kpp_); 	  
    zk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n_,kpp_);
    u = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(m_,1); 	  
    z = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(n_,1);
    bk = HDSA::makePtr<HDSA::Dense_Matrix<RealT> >(1,kpp_); 	  
    a = 0.0;
  }

  virtual ~Model_Error_Kronecker_Matrix()
  { }

  void Set_u(HDSA::Vector<RealT> & u_vec)
  {
    u->Write_Vector_to_Column(0,u_vec);
  }

  void Set_z(HDSA::Vector<RealT> & z_vec)
  {
    z->Write_Vector_to_Column(0,z_vec);
  }

  void Set_uk(int k, HDSA::Vector<RealT> & u_vec)
  {
    uk->Write_Vector_to_Column(k,u_vec);
  }

  void Set_zk(int k, HDSA::Vector<RealT> & z_vec)
  {
    zk->Write_Vector_to_Column(k,z_vec);
  }

  void Set_bk(int k, RealT val)
  {
    bk->Replace_Element(0,k,val);
  }

  void Set_a(RealT val)
  {
    a = val;
  }
 
  void Broadcast_Data(void)
  {
    proc_dist_->Broadcast_Matrix(uk);
    proc_dist_->Broadcast_Matrix(zk);
    proc_dist_->Broadcast_Matrix(bk);
  }

};

}

#endif
