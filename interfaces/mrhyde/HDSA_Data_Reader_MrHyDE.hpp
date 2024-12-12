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

  void Read_Exodus_Data(std::string exofile) {
  
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
    //if (Comm_->MyPID() < 10)
    if (false)
      fname = exofile + "." + strProc + ".0" + strPID;
    else
      fname = exofile + "." + strProc + "." + strPID;
  }
  else {
    fname = exofile;
  }
  
  // open exodus file
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
  
  if (exo_error>0) {
    // need some debug statement
  }
  int id = 1; // only one blkid
  int step = 1; // only one time step (for now)
  ex_block eblock;
  eblock.id = id;
  eblock.type = EX_ELEM_BLOCK;
  
  exo_error = ex_get_block_param(exoid, &eblock);
  
  int num_el_in_blk = eblock.num_entry;
  int num_node_per_el = eblock.num_nodes_per_entry;
    
  // get elem vars
  /*
  bool have_element = false;
  if (have_element) {
    int num_elem_vars = 0;
    int var_ind;
    int numResponses = 1;
    exo_error = ex_get_variable_param(exoid, EX_ELEM_BLOCK, &num_elem_vars); // TMW: this is depracated
    // This turns off this feature
    for (int i=0; i<num_elem_vars; i++) {
      char varname[MAX_STR_LENGTH+1];
      ScalarT *var_vals = new ScalarT[num_el_in_blk];
      var_ind = i+1;
      exo_error = ex_get_variable_name(exoid, EX_ELEM_BLOCK, var_ind, varname);
      string vname(varname);
      efield_names.push_back(vname);
      size_t found = vname.find("Val");
      if (found != std::string::npos) {
        vector<string> results;
        std::stringstream sns, snr;
        int nr;
        results = this->breakupList(vname,"_");
        //boost::split(results, vname, [](char u){return u == '_';});
        snr << results[3];
        snr >> nr;
        numResponses = std::max(numResponses,nr);
      }
      efield_vals.push_back(vector<ScalarT>(num_el_in_blk));
      exo_error = ex_get_var(exoid,step,EX_ELEM_BLOCK,var_ind,id,num_el_in_blk,var_vals);
      for (int j=0; j<num_el_in_blk; j++) {
        efield_vals[i][j] = var_vals[j];
      }
      delete [] var_vals;
    }
  }
  */
  
  int *connect = new int[num_el_in_blk*num_node_per_el];
  int edgeconn, faceconn;
  //exo_error = ex_get_elem_conn(exoid, id, connect);
  exo_error = ex_get_conn(exoid, EX_ELEM_BLOCK, id, connect, &edgeconn, &faceconn);
    
  // get nodal vars
  int num_node_vars = 0;
  int var_ind;
  exo_error = ex_get_variable_param(exoid, EX_NODAL, &num_node_vars);
  // This turns off this feature
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
  std::cout << " proc  "  << comm_->getRank() <<  "  first value " << nfield_vals[0][0] << std::endl;
  std::cout << " proc  "  << comm_->getRank() <<  "  last value " << nfield_vals[0].back() << std::endl;

  ScalarT vecnorm = 0.0;
  for (int i=0;i<nfield_vals[0].size();i++) {
    vecnorm += std::pow(nfield_vals[0][i],2.0); 
  }

  std::cout << "vecnorm   " << std::sqrt(vecnorm) << std::endl;
  
  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > meas = Teuchos::rcp(new Tpetra::MultiVector<ScalarT,LO,GO,SolverNode>(solve_->linalg->overlapped_map[0],1)); // empty solution
  int b = 0;
  //meas->sync<HostDevice>();
  auto meas_kv = meas->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
    
  //meas.modify_host();
  int index, dindex;
    
  //  auto dev_offsets = solve_->assembler->groups[b][0]->wkset->offsets;
  auto dev_offsets = solve_->assembler->wkset[0]->offsets;
  auto offsets = Kokkos::create_mirror_view(dev_offsets);
  Kokkos::deep_copy(offsets,dev_offsets);

  vector<string> blockNames = solve_->mesh->getBlockNames();
  for (int block=0;block<blockNames.size(); block++) {
    for (int grp=0; grp<solve_->assembler->groups[block].size(); ++grp) {
      //cindex = groups[block][grp]->index;
      auto LIDs = solve_->assembler->groups[block][grp]->LIDs_host;
      auto nDOF = solve_->assembler->groups[block][grp]->group_data->num_dof_host;
      
      for (int n=0; n<nDOF(0); n++) {
	//Kokkos::View<GO**,HostDevice> GIDs = assembler->groups[block][grp]->GIDs;
	for (int p=0; p<solve_->assembler->groups[block][grp]->numElem; p++) {
	  for( int i=0; i<nDOF(n); i++ ) {
	    index = LIDs[0](p,offsets(n,i));//cindex(p,n,i);//LA_overlapped_map->getLocalElement(GIDs(p,curroffsets[n][i]));
	    dindex = connect[grp*num_node_per_el + i] - 1;
	    meas_kv(index,0) = nfield_vals[n][dindex];  
	    if(comm_->getRank() == 1) {
	      std::cout << "grp " << grp << " n " << n << " p " << p << " i " << i << " index " << index << "  dindex " << dindex << std::endl;
	    }
	    //(*meas)[0][index] = nfield_vals[n][dindex];
	  }
	}
      }
    }
  }
  //meas.sync<>();
  delete [] connect;
 
  exo_error = ex_close(exoid);
  Teuchos::Array<RealT> val(1,0);
  meas->norm2(val.view(0,1));
  std::cout << "val[0]  " << val[0] << std::endl;
  }
  std::vector<string> breakupList(const string & list, const string & delimiter) {
  // Script to break delimited list into pieces
  string tmplist = list;
  vector<string> terms;
  size_t pos = 0;
  if (tmplist.find(delimiter) == string::npos) {
    terms.push_back(tmplist);
  }
  else {
    string token;
    while ((pos = tmplist.find(delimiter)) != string::npos) {
      token = tmplist.substr(0, pos);
      terms.push_back(token);
      tmplist.erase(0, pos + delimiter.length());
    }
    terms.push_back(tmplist);
  }
  return terms;
  }
};
#endif

