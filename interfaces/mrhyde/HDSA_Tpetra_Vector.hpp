#ifndef HDSA_TPETRA_VECTOR_HPP
#define HDSA_TPETRA_VECTOR_HPP

template <class RealT,
          class LO=Tpetra::Map<>::local_ordinal_type,
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class Tpetra_Vector : public HDSA::Vector<RealT> {
  
private:
  const HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > tpetra_vec_;
  const HDSA::Ptr<const Tpetra::Map<LO,GO,Node> > map_;
  const HDSA::Ptr<const Teuchos::Comm<int> > comm_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;

public:  
  Tpetra_Vector(const HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > &tpetra_vec)
    : tpetra_vec_(tpetra_vec), map_(tpetra_vec_->getMap()), comm_(map_->getComm()),
      random_number_generator_(HDSA::makePtr<HDSA::Random_Number_Generator<RealT> >()) {}
  
  ~Tpetra_Vector()
  { }

  HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > getVector() const {
    return tpetra_vec_;
  }
  //////////////////////////////////////////////////////////////////////////////////
  // Overloading pure virtual functions in HDSA::Vector base class
  //////////////////////////////////////////////////////////////////////////////////


  //Use ROL functionality to fill these function
  HDSA::Ptr<HDSA::Vector<RealT> > clone() const
  {
    int n = tpetra_vec_->getNumVectors();
    return HDSA::makePtr<Tpetra_Vector>(HDSA::makePtr<Tpetra::MultiVector<RealT,LO,GO,Node>>(map_,n));
  }

  // compute the dot product of this and x
  RealT dot( const HDSA::Vector<RealT> &x ) const
  {
    const Tpetra_Vector &ex = dynamic_cast<const Tpetra_Vector&>(x);
    int n = tpetra_vec_->getNumVectors();
    // Perform Euclidean dot between *this and x for each vector
    Teuchos::Array<RealT> val(n,0);
    tpetra_vec_->dot( *ex.getVector(), val.view(0,n) );
    // Combine dots for each vector to get a scalar
    RealT xy(0);
    for (int i = 0; i < n; ++i) {
      xy += val[i];
    }
    return xy;
  }

  // add alpha*x to this
  void axpy( const RealT alpha, const HDSA::Vector<RealT> &x )
  {
    RealT one(1);
    const Tpetra_Vector &ex = dynamic_cast<const Tpetra_Vector&>(x);
    tpetra_vec_->update(alpha,*ex.getVector(),one);
  }

  // return vector dimension
  int dimension() const
  {
    int nVecs = static_cast<int>(tpetra_vec_->getNumVectors());
    int dim   = static_cast<int>(tpetra_vec_->getGlobalLength());
    return nVecs*dim;
  }

  // set this=val elementwise
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
