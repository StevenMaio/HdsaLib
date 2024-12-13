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
  
public:
  MD_Data_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int> > &comm, Teuchos::RCP<MrHyDE::SolverManager<SolverNode> > &solve, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT> > &random_number_generator):
    comm_(comm), solve_(solve), random_number_generator_(random_number_generator)
  { }

  virtual ~MD_Data_Interface_MrHyDE()
  { }

  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_u(void) const{
    std::string exofile = "output.exo";
    Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > u_tpetra =  Read_Exodus_Data(exofile); 
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt_tmp = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
    State_Vector_MrHyDE<RealT> &eu_opt_tmp = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*u_opt_tmp);
    eu_opt_tmp.mrhyde_state_vec[0][0]->update(1.0,*u_tpetra,0.0); 
    
    HDSA::Ptr<HDSA::Vector<RealT> > u_opt = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
    int num_coeff_load = 200;
    std::vector<RealT> opt_u_coeff = std::vector<RealT>(num_coeff_load);
    std::ifstream in("u_opt.txt");
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> opt_u_coeff[j];
        }
      }
    else
     {
        std::cout << "Error loading the data from u_opt.txt" << std::endl;
     }

    State_Vector_MrHyDE<RealT> &eu_opt = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*u_opt);
    for(int k = 0; k < num_coeff_load; k++)
        {
          eu_opt.mrhyde_state_vec[0][0]->replaceGlobalValue(k,0,opt_u_coeff[k]);
        }
    u_opt_tmp->axpy(-1.0,*u_opt);
    std::cout << "u_opt_tmp norm   " << u_opt_tmp->norm() << std::endl;
    
    return u_opt;
  } 
  
  HDSA::Ptr<HDSA::Vector<RealT> > Load_Optimal_z(void) const{
    int num_coeff_load = 200;
    ROL::Ptr<ROL::Vector<RealT> > z_opt_rol = solve_->params->getCurrentVector().clone();
    MrHyDE_OptVector &z_opt = dynamic_cast<MrHyDE_OptVector&>(*z_opt_rol);
    
    std::vector<RealT> opt_z_coeff = std::vector<RealT>(num_coeff_load);
    std::ifstream in("z_opt.txt");
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> opt_z_coeff[j];
        }
      }
    else
     {
        std::cout << "Error loading the data from z_opt.txt" << std::endl;
     }

    for(int k = 0; k < num_coeff_load; k++)
        {
	  z_opt.getField()[0]->getVector()->replaceGlobalValue(k,0,opt_z_coeff[k]);
        }
    HDSA::Ptr<HDSA::Vector<RealT> > z_opt_hdsa = HDSA::makePtr<Opt_Vector_MrHyDE<RealT> >(z_opt, random_number_generator_);
    Opt_Vector_MrHyDE<RealT> &ez_opt_hdsa = dynamic_cast<Opt_Vector_MrHyDE<RealT>&>(*z_opt_hdsa);
    ez_opt_hdsa.mrhyde_vec->set(z_opt);	
    return z_opt_hdsa; 
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_Z_Data(void) const{
    int num_coeff_load = 200;

    ROL::Ptr<ROL::Vector<RealT> > z1_rol = solve_->params->getCurrentVector().clone();
    ROL::Ptr<ROL::Vector<RealT> > z2_rol = solve_->params->getCurrentVector().clone();

    HDSA::Ptr<MrHyDE_OptVector> ez1 = HDSA::dynamicPtrCast<MrHyDE_OptVector> (z1_rol);
    HDSA::Ptr<MrHyDE_OptVector> ez2 = HDSA::dynamicPtrCast<MrHyDE_OptVector> (z2_rol);

    std::ifstream in("Z.txt");
    RealT val = 0.0;
    if (in)
      {
        for(int j = 0; j < num_coeff_load; j++)
        {
          in >> val;
	  ez1->getField()[0]->getVector()->replaceGlobalValue(j,0,val);
          in >> val;
          ez2->getField()[0]->getVector()->replaceGlobalValue(j,0,val);
        }
      }
    else
      {
        std::cout << "Error loading the data from Z.txt" << std::endl;
      }

    HDSA::Ptr<HDSA::Vector<RealT> > z1 = HDSA::makePtr<Opt_Vector_MrHyDE<RealT> >(ez1, random_number_generator_);
    HDSA::Ptr<HDSA::Vector<RealT> > z2 = HDSA::makePtr<Opt_Vector_MrHyDE<RealT> >(ez2, random_number_generator_);

    HDSA::Ptr<HDSA::MultiVector<RealT> > Z = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*z1);
    (*Z)[0]->set(*z1);
    (*Z)[1]->set(*z2);
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT> > Load_D_Data(void) const{
    int num_coeff_load = 200;
    HDSA::Ptr<HDSA::Vector<RealT> > d1 = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
    HDSA::Ptr<HDSA::Vector<RealT> > d2 = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
    State_Vector_MrHyDE<RealT> &ed1 = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*d1);
    State_Vector_MrHyDE<RealT> &ed2 = dynamic_cast<State_Vector_MrHyDE<RealT>&>(*d2);
    std::ifstream in("D.txt");
    RealT val = 0.0;
    if (in)
      {   
        for(int j = 0; j < num_coeff_load; j++)
        {   
          in >> val;
          ed1.mrhyde_state_vec[0][0]->replaceGlobalValue(j,0,val);
          in >> val;
          ed2.mrhyde_state_vec[0][0]->replaceGlobalValue(j,0,val);
        }   
      }   
    else
     {   
        std::cout << "Error loading the data from D.txt" << std::endl;
     }   

    HDSA::Ptr<HDSA::MultiVector<RealT> > D = HDSA::makePtr<HDSA::MultiVector<RealT> >(2,*d1);
    (*D)[0]->set(*d1);
    (*D)[1]->set(*d2);
    return D;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT,LO,GO,SolverNode> > Read_Exodus_Data(std::string exofile) const {
  
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
