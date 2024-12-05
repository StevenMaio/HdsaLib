#ifndef HDSA_VECTOR_MRHYDE_HPP
#define HDSA_VECTOR_MRHYDE_HPP

template <class RealT>
class Vector_MrHyDE : public HDSA::Vector<RealT> {

private:
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;

public:
  HDSA::Ptr<ROL::Vector<RealT> > mrhyde_vec;
  
  Vector_MrHyDE(const MrHyDE_OptVector &mrhyde_vec_in,const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator) :  random_number_generator_(random_number_generator)
  {
    mrhyde_vec = mrhyde_vec_in.clone();
  }

  Vector_MrHyDE(HDSA::Ptr<MrHyDE_OptVector> &mrhyde_vec_in,const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator) : mrhyde_vec(mrhyde_vec_in), random_number_generator_(random_number_generator)
  {
  }

  virtual ~Vector_MrHyDE()
  { }

  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    MrHyDE_OptVector* emrhyde_vec = dynamic_cast<MrHyDE_OptVector*>(&const_cast<ROL::Vector<RealT> &>(*mrhyde_vec));
    HDSA::Ptr<HDSA::Vector<RealT> > hdsa_vector = HDSA::makePtr<Vector_MrHyDE<RealT> >(*emrhyde_vec, random_number_generator_);
    return hdsa_vector;
  }

  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const Vector_MrHyDE<RealT> &ex = dynamic_cast<const Vector_MrHyDE<RealT>&>(x);
    RealT val = ex.mrhyde_vec->dot(*mrhyde_vec);
    return val;
  }

  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    const Vector_MrHyDE<RealT> &ex = dynamic_cast<const Vector_MrHyDE<RealT>&>(x);
    mrhyde_vec->axpy(alpha,*ex.mrhyde_vec);
  }

  int dimension() const
  {
    return mrhyde_vec->dimension();
  }
  
  void setScalar( const RealT val )
  {
    mrhyde_vec->setScalar(val);
  }
  
  void randomize_standard_normal( )
  {
   MrHyDE_OptVector* emrhyde_vec = dynamic_cast<MrHyDE_OptVector*>(&const_cast<ROL::Vector<RealT> &>(*mrhyde_vec));
   if(emrhyde_vec->getField()[0] != ROL::nullPtr) {
     for (int j = 0; j<emrhyde_vec->getField().size();j++) {
       auto vecT_data = emrhyde_vec->getField()[j]->getVector()->getDataNonConst(0);
       for (int i = 0; i< emrhyde_vec->getField()[j]->getVector()->getLocalLength(); i++) {
	 vecT_data[i] = random_number_generator_->Generate_Standard_Normal_Sample();   
       }
     }
   }
   if(emrhyde_vec->getParameter()[0] != ROL::nullPtr) {
     for (int j=0; j<emrhyde_vec->getParameter().size();j++) {
       for (int i=0; i< emrhyde_vec->getParameter()[j]->getVector()->size(); i++) {
	 (*emrhyde_vec->getParameter()[j]->getVector())[i] = random_number_generator_->Generate_Standard_Normal_Sample();   
       }
     }
   }
  }

  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > Get_random_number_generator() const {
    return random_number_generator_;
  }
};
#endif
