#ifndef HDSA_DATA_READER_MRHYDE_HPP
#define HDSA_DATA_READER_MRHYDE_HPP

#include "exodusII.h"
#include "preferences.hpp"

template <class RealT,
	  class LO=Tpetra::Map<>::local_ordinal_type,
	  class GO=Tpetra::Map<>::global_ordinal_type,
	  class Node=Tpetra::Map<>::node_type >
class Data_Reader_MrHyDE {

private:
  Teuchos::RCP<Teuchos::MpiComm<int> > comm_;
  Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > solve_;
  
public:
  
  Data_Reader_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int> > &comm,Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve):comm_(comm),solve_(solve)
  { }
    
  virtual ~Data_Reader_MrHyDE()
  { }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > Read_Exodus_Data(std::string exofile) {
  
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
      
      for (int n=0; n<nDOF(0); n++) {
	for (int p=0; p<solve_->assembler->groups[block][grp]->numElem; p++) {
	  for( int i=0; i<nDOF(n); i++ ) {
	    index = LIDs[0](p,offsets(n,i));
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

