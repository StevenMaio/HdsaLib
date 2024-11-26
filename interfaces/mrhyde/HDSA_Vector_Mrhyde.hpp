#ifndef HDSA_VECTOR_MRHYDE_HPP
#define HDSA_VECTOR_MRHYDE_HPP

namespace HDSA
{

template <class RealT>
class Vector_MrHyDE : public Vector<RealT> {

private:
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;
public:
  HDSA::Ptr<ROL::Vector<RealT> > mrhyde_vec;
  
  Vector_MrHyDE(const MrHyDE_OptVector &mrhyde_vec_in) :  random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >())
  {
    mrhyde_vec = mrhyde_vec_in.clone();
  }

  Vector_MrHyDE(HDSA::Ptr<MrHyDE_OptVector> &mrhyde_vec_in):mrhyde_vec(mrhyde_vec_in), random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >())
  {
  }

  virtual ~Vector_MrHyDE()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define when creating a vector interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Clone the vector
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    //    const MrHyDE_OptVector &emrhyde_vec = dynamic_cast<const MrHyDE_OptVector&>(*mrhyde_vec);
    MrHyDE_OptVector* emrhyde_vec;
    emrhyde_vec = dynamic_cast<MrHyDE_OptVector*>(&const_cast<ROL::Vector<RealT> &>(*mrhyde_vec));

    HDSA::Ptr<HDSA::Vector<RealT> > hdsa_vector = HDSA::makePtr<HDSA::Vector_MrHyDE<RealT> >(*emrhyde_vec);
    return hdsa_vector;
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const HDSA::Vector_MrHyDE<RealT> &ex = dynamic_cast<const HDSA::Vector_MrHyDE<RealT>&>(x);
    RealT val = ex.mrhyde_vec->dot(*mrhyde_vec);
    return val;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const HDSA::Vector_MrHyDE<RealT> &ex = dynamic_cast<const HDSA::Vector_MrHyDE<RealT>&>(x);
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
   MrHyDE_OptVector* emrhyde_vec = dynamic_cast<MrHyDE_OptVector*>(&const_cast<ROL::Vector<RealT> &>(*mrhyde_vec));
   for (int j=0; j<emrhyde_vec->getField().size();j++) {
     auto vecT_data = emrhyde_vec->getField()[j]->getVector()->getDataNonConst(0);
     for (int i=0; i< emrhyde_vec->getField()[j]->getVector()->getLocalLength(); i++) {
       vecT_data[i] = random_number_generator_->Generate_Standard_Normal_Sample();   
     }
   }

   for (int j=0; j<emrhyde_vec->getParameter().size();j++) {
     for (int i=0; i< emrhyde_vec->getParameter()[j]->getVector()->size(); i++) {
       (*emrhyde_vec->getParameter()[j]->getVector())[i] = random_number_generator_->Generate_Standard_Normal_Sample();   
     }
   }
  }
 
 
};

}

#endif
