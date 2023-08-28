#ifndef HDSA_MRHYDE_STATE_VECTOR_HPP
#define HDSA_MRHYDE_STATE_VECTOR_HPP

namespace HDSA
{

template <class RealT,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class Vector_Mrhyde_State : public Vector<RealT> {

public:

  std::vector<std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > > mrhyde_state_vec;
  Vector_Mrhyde_State(Teuchos::RCP<MrHyDE::SolverManager<Node> > &solver) 
  {
    int numsets = solver->setnames.size();
    int numsteps = solver->numsteps;
    mrhyde_state_vec.resize(numsteps);
    for (int i=0; i<numsteps; i++) {
      mrhyde_state_vec[i].resize(numsets);
      for (int set=0; set<numsets; set++) {
	mrhyde_state_vec[i][set] = solver->linalg->getNewVector(set);
      }
    }
  }

  virtual ~Vector_Mrhyde_State()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define when creating a vector interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Clone the vector
  virtual HDSA::Ptr<HDSA::Vector<RealT> > clone() const {
  }

  // compute the dot product of this and x
  virtual RealT dot( const HDSA::Vector<RealT> &x ) const {
  }

  // add alpha*x to this
  virtual void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) {
  }
 
  // return vector dimension
  virtual int dimension() const {
  }

  // set this=val elementwise
  virtual void setScalar( const RealT val ) {
  }

  virtual void randomize_standard_normal( ) {
  }

};

}

#endif
