#ifndef HDSA_MRHYDE_STATE_VECTOR_HPP
#define HDSA_MRHYDE_STATE_VECTOR_HPP

template <class RealT,
          class LO=Tpetra::Map<>::local_ordinal_type, 
          class GO=Tpetra::Map<>::global_ordinal_type,
          class Node=Tpetra::Map<>::node_type >
class Vector_MrHyDE_State : public HDSA::Vector<RealT> {

private:
  Teuchos::RCP<MrHyDE::SolverManager<Node> > solver_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;
  
public:
  std::vector<std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT,LO,GO,Node> > > > mrhyde_state_vec;

  Vector_MrHyDE_State(const Teuchos::RCP<MrHyDE::SolverManager<Node> > &solver,const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator,bool isSteady=false) : solver_(solver), random_number_generator_(random_number_generator)
  {
    int numsets = solver->setnames.size();
    mrhyde_state_vec.resize(numsets);
    for (int set=0; set<numsets; set++) {
      if(isSteady) {
	mrhyde_state_vec[set].resize(1);
	mrhyde_state_vec[set][0] = solver->linalg->getNewVector(set);
      } else {
	mrhyde_state_vec[set].resize(solver->numsteps[set]);
	for (int i = 0; i<solver->numsteps[set]; i++) {
	  mrhyde_state_vec[set][i] = solver->linalg->getNewVector(set);
	}
      }
    }
  }
  
  virtual ~Vector_MrHyDE_State()
  { }

  virtual HDSA::Ptr<HDSA::Vector<RealT> > clone() const {
    HDSA::Ptr<HDSA::Vector<RealT> > vec = HDSA::makePtr<Vector_MrHyDE_State<RealT> >(solver_,random_number_generator_);
    return vec;
  }

  virtual RealT dot( const HDSA::Vector<RealT> &x ) const {
    const Vector_MrHyDE_State<RealT> &ex = dynamic_cast<const Vector_MrHyDE_State<RealT>&>(x);
    int numsets = solver_->setnames.size();
    RealT xy(0);    
    for (int set = 0; set<numsets; set++) {
      int n = ex.mrhyde_state_vec[set].size();
      Teuchos::Array<RealT> val(1,0);
      for (int i = 0; i < n; ++i) {
	mrhyde_state_vec[set][i]->dot(*ex.mrhyde_state_vec[set][i],val.view(0,1)); 
        xy += val[0];
      }
    }
    return xy;
  }

  virtual void axpy( const RealT alpha, const HDSA::Vector<RealT> &x ) {
    const Vector_MrHyDE_State<RealT> &ex = dynamic_cast<const Vector_MrHyDE_State<RealT>&>(x);
    int numsets = solver_->setnames.size();
    RealT one(1);
    for (int set = 0; set<numsets; set++) {
      for (int i = 0; i<solver_->numsteps[set]; i++) {
	mrhyde_state_vec[set][i]->update(alpha,*ex.mrhyde_state_vec[set][i],one);
      }
    }
  }
 
  virtual int dimension() const {
    int spatialdim = 0;
    int numsets = solver_->setnames.size();
    for (int set = 0; set<numsets; set++) {
      for (int i = 0; i<solver_->numsteps[set]; i++) {
	spatialdim += mrhyde_state_vec[set][i]->getGlobalLength();
      }
    }
    return spatialdim;
  }

  virtual void setScalar( const RealT val ) {
    int numsets = solver_->setnames.size();
    for (int set = 0; set<numsets; set++) {
      for (int i = 0; i<solver_->numsteps[set]; i++) {
	mrhyde_state_vec[set][i]->putScalar(val);
      }
    }
  }

  virtual void randomize_standard_normal() {
    int numsets = solver_->setnames.size();
    for (int set = 0; set<numsets; set++) {
      for (int i = 0; i<solver_->numsteps[set]; i++) {
	auto vecT_data = mrhyde_state_vec[set][i]->getDataNonConst(0);
	for (int j = 0; j<mrhyde_state_vec[set][i]->getLocalLength(); j++) {
	  vecT_data[j] = random_number_generator_->Generate_Standard_Normal_Sample();   
	}
      }
    }
  }
  
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > Get_random_number_generator() const {
    return random_number_generator_;
  }
};
#endif
