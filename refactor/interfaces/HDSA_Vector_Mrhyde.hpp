#ifndef HDSA_VECTOR_MRHYDE_HPP
#define HDSA_VECTOR_MRHYDE_HPP

namespace HDSA
{

template <class RealT>
class Vector_Mrhyde : public Vector<RealT> {

public:
  HDSA::Ptr<MrHyDE_OptVector> mrhyde_vec;
  
  Vector_Mrhyde(MrHyDE_OptVector &mrhyde_vec_in)
  {
    mrhyde_vec = mrhyde_vec_in.clone();
  }

  virtual ~Vector_Mrhyde()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define when creating a vector interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    HDSA::Ptr<HDSA::Vector<RealT> > hdsa_vector = HDSA::makePtr<HDSA::Vector_Mrhyde<RealT> >(*mrhyde_vec);
    return hdsa_vector;
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const HDSA::Vector_Mrhyde<RealT> &ex = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(x);
    RealT val = ex.mrhyde_vec->dot(mrhyde_vec);
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const HDSA::Vector_Mrhyde<RealT> &ex = dynamic_cast<const HDSA::Vector_Mrhyde<RealT>&>(x);
    mrhyde_vec->axpy(alpha,*ex.mrhyde_vec);
  }
  // return vector dimension
  int dimension() const
  {
    return mrhyde_vec->dimension();
  }
  

  // set this=val elementwise
  void setScalar( const RealT val )
  {
  mrhyde_vec->setScalar(val);
  }
  
  void randomize_standard_normal( )
  {
    // bvbw: using uniform to get everthing mechanically working - need Normal 
    mrhyde_vec->randomize();
  }
 
 
};

}

#endif
