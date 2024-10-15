#ifndef HDSA_MRHYDE_STEADY_STATE_VECTOR_HPP
#define HDSA_MRHYDE_STEADY_STATE_VECTOR_HPP

namespace HDSA
{

template <class RealT,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class Vector_MrHyDE_Steady_State : public Vector<RealT> {

private:
  Teuchos::RCP<MrHyDE::SolverManager<Node> > solver_;
public:

  std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > mrhyde_steady_state_vec;
  Vector_MrHyDE_Steady_State(const Teuchos::RCP<MrHyDE::SolverManager<Node> > &solver) : solver_(solver)
  {
    int numsets = solver->setnames.size();
    mrhyde_steady_state_vec.resize(numsets);
    for (int set=0; set<numsets; set++) {
	mrhyde_steady_state_vec[set] = solver->linalg->getNewVector(set);
      }
  }

  virtual ~Vector_MrHyDE_Steady_State()
  { }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Pure virtual functions to define when creating a vector interface
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  // Clone the vector
  virtual HDSA::Ptr<HDSA::Vector<RealT> > clone() const {

    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr< HDSA::Vector_MrHyDE_Steady_State<RealT> >(solver_);
    return vec;
  }

  // compute the dot product of this and x
  virtual RealT dot( const HDSA::Vector<RealT> &x ) const {
    // const HDSA::Vector_MrHyDE_Steady_State<RealT> &ex = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(x);
    // int numsets = solver_->setnames.size();
    // RealT val = 0.0;
    // for (int set=0; set<numsets; set++) {
    //   for (int i=0; i<solver_->numsteps[set]; i++) {
    //   	val += ex.mrhyde_steady_state_vec[set][i]->dot(*mrhyde_steady_state_vec[set][i]);
    //         }
    // }
    // return val;
    // modeled after dot in ROL::pdevector.hpp
    const HDSA::Vector_MrHyDE_Steady_State<RealT> &ex = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(x);
    int numsets = solver_->setnames.size();
    RealT xy(0);    
    for (int set=0; set<numsets; set++) {
      Teuchos::Array<RealT> val(1,0);
	mrhyde_steady_state_vec[set]->dot(*ex.mrhyde_steady_state_vec[set],val.view(0,1)); 
        xy += val[0];
    }
    return xy;
  }

  // add alpha*x to this
  virtual void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) {
    const HDSA::Vector_MrHyDE_Steady_State<RealT> &ex = dynamic_cast<const HDSA::Vector_MrHyDE_Steady_State<RealT>&>(x);
    int numsets = solver_->setnames.size();
    RealT one(1);
    for (int set=0; set<numsets; set++) {
	//	mrhyde_steady_state_vec[set][i]->axpy(alpha,*ex.steady_mrhyde_steady_state_vec[set][i]);
	// from ROL_Tpetra_Multivector.hpp
	mrhyde_steady_state_vec[set]->update(alpha,*ex.mrhyde_steady_state_vec[set],one);
      }
    }
 
  // return vector dimension
  virtual int dimension() const {
    int spatialdim = 0;
    int numsets = solver_->setnames.size();
    for (int set=0; set<numsets; set++) {
	spatialdim += mrhyde_steady_state_vec[set]->getGlobalLength();
	//    return numsets*numsteps*spatialdim;
    }
    return spatialdim;
  }

  // set this=val elementwise
  virtual void setScalar( const RealT val ) {
    int numsets = solver_->setnames.size();
    for (int set=0; set<numsets; set++) {
	mrhyde_steady_state_vec[set]->putScalar(val);
      }
    }

  //  virtual void randomize_standard_normal(RealT l = 0.0, RealT u = 1.0) {
  // bvbw need normal
  virtual void randomize_standard_normal() {
    int numsets = solver_->setnames.size();
    for (int set=0; set<numsets; set++) {
	//	mrhyde_state_vec[set][i]->randomize(l,u);
	mrhyde_steady_state_vec[set]->randomize();
      }
    }
};

}

#endif
