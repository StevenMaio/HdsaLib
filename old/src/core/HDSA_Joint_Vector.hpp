#ifndef HDSA_JOINT_VECTOR_HPP
#define HDSA_JOINT_VECTOR_HPP

namespace HDSA
{

template <class RealT>
class Joint_Vector : public HDSA::Vector<RealT> {

  HDSA::Ptr<HDSA::Vector<RealT> > vec_1_;
  HDSA::Ptr<HDSA::Vector<RealT> > vec_2_;
  int L1_;

public:

  Joint_Vector(HDSA::Ptr<HDSA::Vector<RealT> > vec_1, HDSA::Ptr<HDSA::Vector<RealT> > vec_2) 
  {
    vec_1_ = vec_1->Clone();
    vec_2_ = vec_2->Clone();
    L1_ = vec_1_->dimension();
    HDSA::Vector<RealT>::enforce_zeros_ = vec_1_->Get_enforce_zeros();
    if(HDSA::Vector<RealT>::enforce_zeros_)
      {
        HDSA::Vector<RealT>::map_full_to_reduced_ = vec_1_->Get_map_full_to_reduced();
        HDSA::Vector<RealT>::map_reduced_to_full_ = vec_1_->Get_map_reduced_to_full();
        HDSA::Vector<RealT>::vec_zeros_ = vec_1_->Get_vec_zeros();
      }
  }

  virtual ~Joint_Vector()
  { }

  // Access the (i,j) element
  RealT operator () (int k) const
  {
    RealT val;
    if(k < L1_)
      {
	val = (*vec_1_)(k);
      }
    else
      {
	val = (*vec_2_)(k-L1_);
      }
    return val;
  }
  
  // Replace the kth element of the vector by val
  void Replace_Element(int k, RealT val) 
  {
    if(k < L1_)
      {
	vec_1_->Replace_Element(k,val);	
      }
    else
      {
	vec_2_->Replace_Element(k-L1_,val);
      }
  }

  // Get the data on this processor
  std::vector<RealT> Get_Data_on_Processor(void) const
  {
    std::vector<RealT> data;
    return data;
  }

  // Get the indices on this processor
  std::vector<int> Get_Indices_on_Processor(void) const 
  {
    std::vector<int> indices;
    return indices;
  }

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > Define_Clone() const 
  {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Joint_Vector<RealT> >(vec_1_,vec_2_);
    return vec;
  }

  // add x to this
  void plus( const HDSA::Vector<RealT> & x ) 
  {
    const Joint_Vector<RealT> &ex = dynamic_cast<const Joint_Vector<RealT>&>(x);
    vec_1_->plus(*ex.Get_Component_Vector_1());
    vec_2_->plus(*ex.Get_Component_Vector_2());
  }

  // scale this by val
  void scale( const RealT val ) 
  {
    vec_1_->scale(val);
    vec_2_->scale(val);
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const Joint_Vector<RealT> &ex = dynamic_cast<const Joint_Vector<RealT>&>(x);
    RealT val = vec_1_->dot(*ex.Get_Component_Vector_1());
    val += vec_2_->dot(*ex.Get_Component_Vector_2());
    return val;
  }
 
  // compute the norm of this
  RealT norm(void) const 
  {
    return std::sqrt(this->dot(*this));
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) 
  {
    const Joint_Vector<RealT> &ex = dynamic_cast<const Joint_Vector<RealT>&>(x);
    vec_1_->axpy(alpha,*ex.Get_Component_Vector_1());
    vec_2_->axpy(alpha,*ex.Get_Component_Vector_2());
  }
 
  // set this=0
  void zero(void) 
  {
    this->scale(0.0);
  }

  // set this= ith canonical basis vector
  void basis( const int i ) 
  {
    vec_1_->zero();
    vec_2_->zero();
    if(i < L1_)
      {
	vec_1_->basis(i);
      }
    else
      {
	vec_2_->basis(i-L1_);
      }
  }
 
  // return vector dimension
  int dimension() const
  {
    return vec_1_->dimension() + vec_2_->dimension();
  }

  // set this=x
  void set( const HDSA::Vector<RealT> &x ) 
  {
    const Joint_Vector<RealT> &ex = dynamic_cast<const Joint_Vector<RealT>&>(x);
    vec_1_->set(*ex.Get_Component_Vector_1());
    vec_2_->set(*ex.Get_Component_Vector_2());
  }

  // set this=val elementwise
  void setScalar( const RealT val ) 
  {
    vec_1_->setScalar(val);
    vec_2_->setScalar(val);
  }

  // set entries of this to random numbers in [l,u]
  void randomize( const RealT l = 0.0, const RealT u = 1.0 ) 
  {
    vec_1_->randomize(l,u);
    vec_2_->randomize(l,u);
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Get_Component_Vector_1(void) const
  {
    return vec_1_;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > Get_Component_Vector_2(void) const
  {
    return vec_2_;
  }
 
};

}

#endif
