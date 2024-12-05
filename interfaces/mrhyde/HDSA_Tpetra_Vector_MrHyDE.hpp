#ifndef HDSA_TPETRA_VECTOR_MRHYDE_HPP
#define HDSA_TPETRA_VECTOR_MRHYDE_HPP

template <class RealT,
          class LO=Tpetra::Map<>::local_ordinal_type,
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class Tpetra_Vector_MrHyDE : public HDSA::Vector<RealT> {
  
private:
  const HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > tpetra_vec_;
  const HDSA::Ptr<const Tpetra::Map<LO,GO,Node> > map_;
  const HDSA::Ptr<const Teuchos::Comm<int> > comm_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;

public:  
  Tpetra_Vector_MrHyDE(const HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > &tpetra_vec, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator)
    : tpetra_vec_(tpetra_vec), map_(tpetra_vec_->getMap()), comm_(map_->getComm()),
      random_number_generator_(random_number_generator) {}
  
  ~Tpetra_Vector_MrHyDE()
  { }

  HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > getVector() const {
    return tpetra_vec_;
  }

  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    int n = tpetra_vec_->getNumVectors();
    return HDSA::makePtr<Tpetra_Vector_MrHyDE>(HDSA::makePtr<Tpetra::MultiVector<RealT,LO,GO,Node>>(map_,n),random_number_generator_);
  }

  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const Tpetra_Vector_MrHyDE &ex = dynamic_cast<const Tpetra_Vector_MrHyDE&>(x);
    int n = tpetra_vec_->getNumVectors();
    Teuchos::Array<RealT> val(n,0);
    tpetra_vec_->dot( *ex.getVector(), val.view(0,n) );
    RealT xy(0);
    for (int i = 0; i < n; ++i) {
      xy += val[i];
    }
    return xy;
  }

  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    RealT one(1);
    const Tpetra_Vector_MrHyDE &ex = dynamic_cast<const Tpetra_Vector_MrHyDE&>(x);
    tpetra_vec_->update(alpha,*ex.getVector(),one);
  }

  int dimension() const
  {
    int nVecs = static_cast<int>(tpetra_vec_->getNumVectors());
    int dim   = static_cast<int>(tpetra_vec_->getGlobalLength());
    return nVecs*dim;
  }

  void setScalar( const RealT val )
  {
    tpetra_vec_->putScalar(static_cast<double>(val));
  }

  void randomize_standard_normal( ) 
  {
    int nVecs = static_cast<int>(tpetra_vec_->getNumVectors());
    for (int j=0; j<nVecs;j++) {
      auto vecT_data = tpetra_vec_->getDataNonConst(j);
      for (int i=0; i<tpetra_vec_->getLocalLength(); i++) {
	vecT_data[i] = random_number_generator_->Generate_Standard_Normal_Sample();   
      }
    }
  }
};


#endif
