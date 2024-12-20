#ifndef HDSA_MD_DATA_INTERFACE_MRHYDE_HPP
#define HDSA_MD_DATA_INTERFACE_MRHYDE_HPP

#include "exodusII.h"
#include "preferences.hpp"

template <class RealT,
	  class LO=Tpetra::Map<>::local_ordinal_type,
	  class GO=Tpetra::Map<>::global_ordinal_type,
	  class Node=Tpetra::Map<>::node_type >
class MD_Data_Interface_MrHyDE : public HDSA::MD_Data_Interface<RealT> {

private:
  Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
  Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solve_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > random_number_generator_;
  Teuchos::ParameterList data_load_list_;
  std::string opt_var_physics_;
  std::string opt_solution_exo_file_;
  int num_hifi_;
  std::vector<std::string> hifi_exo_files_;
  std::vector<std::string> lofi_exo_files_; 
  
public:
  MD_Data_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator, Teuchos::ParameterList &data_load_list):
    comm_(comm), solve_(solve), random_number_generator_(random_number_generator), data_load_list_(data_load_list)
  {
    opt_var_physics_ = data_load_list_.get<std::string>("OptVariablePhysics", "error");
    opt_solution_exo_file_ = data_load_list_.get<std::string>("OptimalSolutionExoFile", "error");
    num_hifi_ = data_load_list_.get<int>("NumHifi", 1);
    hifi_exo_files_.resize(num_hifi_);
    lofi_exo_files_.resize(num_hifi_);
    for (int k=0;k<num_hifi_;k++) {
      hifi_exo_files_[k] = data_load_list_.get<std::string>("HifiExoFile" + std::to_string(k+1), "error");
      lofi_exo_files_[k] = data_load_list_.get<std::string>("LofiExoFile" + std::to_string(k+1), "error");
    }
    
    if(opt_var_physics_ == "error") {
      std::cout << "Error: specify OptVariablePhysics" << std::endl;
    }
    if(opt_solution_exo_file_ == "error") {
      std::cout << "Error: specify OptimalSolutionExoFile" << std::endl;
    }
    for (int k=0;k<num_hifi_;k++) {
      if(hifi_exo_files_[k] == "error") {
	std::cout << "Error: specify HifiExoFile" + std::to_string(k+1) << std::endl;
      }
      if(lofi_exo_files_[k] == "error") {
	std::cout << "Error: specify fiExoFile" + std::to_string(k+1) << std::endl;
      }
    }
  }

  
  virtual ~MD_Data_Interface_MrHyDE()
  { }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const{

    Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > u_tpetra =  Read_Exodus_Data(opt_solution_exo_file_); 
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
    State_Vector_MrHyDE<RealT> &eu_opt = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*u_opt);
    eu_opt.mrhyde_state_vec[0][0]->update(1.0,*u_tpetra,0.0); 
    
    return u_opt;
  } 
  
  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const{
    
    Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > z_tpetra =  Read_Exodus_Data(opt_solution_exo_file_,false); 

    ROL::Ptr<ROL::Vector<RealT> > z_opt_rol = solve_->params->getCurrentVector().clone();
    HDSA::Ptr<MrHyDE_OptVector> z_opt = HDSA::dynamicPtrCast<MrHyDE_OptVector>(z_opt_rol);
    z_opt->getField()[0]->getVector()->update(1.0,*z_tpetra,0.0);

    HDSA::Ptr<HDSA::Vector<RealT> > z_opt_hdsa = HDSA::makePtr<Opt_Vector_MrHyDE<RealT> >(z_opt, random_number_generator_);

    return z_opt_hdsa; 
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const{
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > z_vecs;
    for(int k=0; k<num_hifi_; k++) {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > z_tpetra =  Read_Exodus_Data(lofi_exo_files_[k],false); 
    
      ROL::Ptr<ROL::Vector<RealT> > z_opt_rol = solve_->params->getCurrentVector().clone();
      HDSA::Ptr<MrHyDE_OptVector> z_opt = HDSA::dynamicPtrCast<MrHyDE_OptVector>(z_opt_rol);
      z_opt->getField()[0]->getVector()->update(1.0,*z_tpetra,0.0);
      
      HDSA::Ptr<HDSA::Vector<RealT> > z_opt_hdsa = HDSA::makePtr<Opt_Vector_MrHyDE<RealT> >(z_opt, random_number_generator_);
      z_vecs.push_back(z_opt_hdsa);
    }

    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(z_vecs);
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const{
    std::vector<HDSA::Ptr<HDSA::Vector<RealT> > > u_vecs; 
    for(int k=0; k<num_hifi_; k++) {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > u_tpetra_hifi =  Read_Exodus_Data(hifi_exo_files_[k]);
      Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > u_tpetra_lofi =  Read_Exodus_Data(lofi_exo_files_[k]); 
      HDSA::Ptr<HDSA::Vector<RealT> > u_vec_k = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
      State_Vector_MrHyDE<RealT> &eu_vec_k = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*u_vec_k);
      eu_vec_k.mrhyde_state_vec[0][0]->update(1.0,*u_tpetra_hifi,0.0);
      eu_vec_k.mrhyde_state_vec[0][0]->update(-1.0,*u_tpetra_lofi,1.0);
      u_vecs.push_back(u_vec_k);
    }
    HDSA::Ptr<HDSA::MultiVector<RealT> > D = HDSA::makePtr<HDSA::MultiVector<RealT> >(u_vecs);
    
    return D;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > Read_Exodus_Data(std::string exofile, bool load_state = true) const {
  
  string fname;

  std::vector<std::vector<ScalarT> > nfield_vals, efield_vals;
  std::vector<std::string> block_names, side_names, node_names, nfield_names, efield_names;
  
  if (comm_->getSize() > 1) {
    std::stringstream ssProc, ssPID;
    ssProc << comm_->getSize();
    ssPID << comm_->getRank();
    string strProc = ssProc.str();
    string strPID = ssPID.str();
    // this section may need tweaking if the input exodus mesh is
    // spread across 10's, 100's, or 1000's (etc) of processors
    fname = exofile + "." + strProc + "." + strPID;
  }
  else {
    fname = exofile;
  }
  
  int CPU_word_size, IO_word_size, exoid, exo_error;
  int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;
  char title[MAX_STR_LENGTH+1];
  float exo_version;
  CPU_word_size = sizeof(ScalarT);
  IO_word_size = 0;
  exoid = ex_open(fname.c_str(), EX_READ, &CPU_word_size,&IO_word_size,
                  &exo_version);
  exo_error = ex_get_init(exoid, title, &num_dim, &num_nods, &num_el,
                          &num_el_blk, &num_ns, &num_ss);
  
  int id = 1; 
  int step = 1;
  ex_block eblock;
  eblock.id = id;
  eblock.type = EX_ELEM_BLOCK;
  
  exo_error = ex_get_block_param(exoid, &eblock);
  
  int num_el_in_blk = eblock.num_entry;
  int num_node_per_el = eblock.num_nodes_per_entry;
    
  int *connect = new int[num_el_in_blk*num_node_per_el];
  int edgeconn, faceconn;
  exo_error = ex_get_conn(exoid, EX_ELEM_BLOCK, id, connect, &edgeconn, &faceconn);
    
  int num_node_vars = 0;
  int var_ind;
  exo_error = ex_get_variable_param(exoid, EX_NODAL, &num_node_vars);

  for (int i=0; i<num_node_vars; i++) {
    char varname[MAX_STR_LENGTH+1];
    ScalarT *var_vals = new ScalarT[num_nods];
    var_ind = i+1;
    exo_error = ex_get_variable_name(exoid, EX_NODAL, var_ind, varname);
    string vname(varname);
    nfield_names.push_back(vname);
    nfield_vals.push_back(vector<ScalarT>(num_nods));
    exo_error = ex_get_var(exoid,step,EX_NODAL,var_ind,0,num_nods,var_vals);
    for (int j=0; j<num_nods; j++) {
      nfield_vals[i][j] = var_vals[j];
    }
    delete [] var_vals;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > vec_over = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,LO,GO,SolverNode>(solve_->linalg->overlapped_map[0],1));
  auto vec_over_kv = vec_over->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
  int index, dindex;
    
  auto dev_offsets = solve_->assembler->wkset[0]->offsets;
  auto offsets = Kokkos::create_mirror_view(dev_offsets);
  Kokkos::deep_copy(offsets,dev_offsets);

  vector<string> blockNames = solve_->mesh->getBlockNames();

  
  
  for (int block=0;block<blockNames.size(); block++) {
    for (int grp=0; grp<solve_->assembler->groups[block].size(); ++grp) {
      auto LIDs = solve_->assembler->groups[block][grp]->LIDs_host;
      auto nDOF = solve_->assembler->groups[block][grp]->group_data->num_dof_host;

      int bound_nv = solve_->physics->getVarList()[0][block].size();
      if(!load_state) {
	bound_nv = solve_->params->getParamsNames(4).size();
      }
      
      for (int nv=0; nv<bound_nv; nv++) {
	std::string var;
	if (load_state) {
	  var = solve_->physics->getVarList()[0][block][nv];
	} else {
	  var = solve_->params->getParamsNames(4)[nv];
	}
	int n = -1;

	for (int j= 0; j<nfield_names.size();j++) {
	  if(nfield_names[j] == var) {
	    n = j;
	  } 
	}
	int n_var = n;
	if (!load_state) {
	  for (int j= 0; j<nfield_names.size();j++) {
	    if(nfield_names[j] == opt_var_physics_) {
	      n_var = j;
	    } 
	  }
	}
	for (int p=0; p<solve_->assembler->groups[block][grp]->numElem; p++) {
	  for( int i=0; i<nDOF(n_var); i++ ) {
	    index = LIDs[0](p,offsets(n_var,i));
	    dindex = connect[grp*num_node_per_el + i] - 1;
	    vec_over_kv(index,0) = nfield_vals[n][dindex];
	  }
	}
      }
    }
  }
  delete [] connect;

  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > vec = solve_->linalg->getNewVector(0);
  solve_->linalg->exportVectorFromOverlappedReplace(0, vec, vec_over);
  
  exo_error = ex_close(exoid);
  return vec;
  }

};
#endif
