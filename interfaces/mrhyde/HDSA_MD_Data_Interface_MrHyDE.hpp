#ifndef HDSA_MD_DATA_INTERFACE_MRHYDE_HPP
#define HDSA_MD_DATA_INTERFACE_MRHYDE_HPP

#include "exodusII.h"
#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_Tpetra_Vector.hpp"
#include "HDSA_Std_Vector.hpp"

template <class RealT,
          class LO = Tpetra::Map<>::local_ordinal_type,
          class GO = Tpetra::Map<>::global_ordinal_type,
          class Node = Tpetra::Map<>::node_type>
class MD_Data_Interface_MrHyDE : public HDSA::MD_Data_Interface<RealT>
{

private:
  Teuchos::RCP<Teuchos::MpiComm<int>> comm_;
  Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solve_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;
  Teuchos::ParameterList data_load_list_;
  int num_hifi_;

  std::string opt_solution_exo_file_;
  std::vector<std::string> hifi_exo_files_;

  std::string opt_solution_txt_file_u_;
  std::string opt_solution_txt_file_z_;
  std::vector<std::string> hifi_txt_files_u_;
  std::vector<std::string> txt_files_z_;

public:
  MD_Data_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int>> &comm, Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solve, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator, Teuchos::ParameterList &data_load_list) : comm_(comm), solve_(solve), params_(params), random_number_generator_(random_number_generator), data_load_list_(data_load_list)
  {
    num_hifi_ = data_load_list_.get<int>("NumHifi", 1);
    opt_solution_exo_file_ = data_load_list_.get<std::string>("OptimalSolutionExoFile", "error");

    opt_solution_txt_file_u_ = data_load_list_.get<std::string>("OptimalSolutionTxtFileU", "error");
    opt_solution_txt_file_z_ = data_load_list_.get<std::string>("OptimalSolutionTxtFileZ", "error");

    hifi_exo_files_.resize(num_hifi_);

    hifi_txt_files_u_.resize(num_hifi_);
    txt_files_z_.resize(num_hifi_);
    for (int k = 0; k < num_hifi_; k++)
    {
      hifi_exo_files_[k] = data_load_list_.get<std::string>("HifiExoFile" + std::to_string(k + 1), "error");

      hifi_txt_files_u_[k] = data_load_list_.get<std::string>("HifiTxtFileU" + std::to_string(k + 1), "error");
      txt_files_z_[k] = data_load_list_.get<std::string>("TxtFileZ" + std::to_string(k + 1), "error");
    }
  }

  virtual ~MD_Data_Interface_MrHyDE()
  {
  }

  Teuchos::RCP<Teuchos::MpiComm<int>> Get_Communicator(void) const
  {
    return comm_;
  }

  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> Get_Parameter_Manager(void) const
  {
    return params_;
  }

  std::string Get_Opt_Solution_Exo_File(void) const
  {
    return opt_solution_exo_file_;
  }

  void Overwrite_Opt_Solution_Exo_File(std::string &exo_file)
  {
    opt_solution_exo_file_ = exo_file;
  }

  std::vector<std::string> Get_HiFi_Exo_Files(void) const
  {
    return hifi_exo_files_;
  }

  void Overwrite_HiFi_Exo_Files(std::vector<std::string> &exo_files)
  {
    for (int k = 0; k < exo_files.size(); k++)
    {
      hifi_exo_files_[k] = exo_files[k];
    }
  }

  void State_Solve(HDSA::Vector<RealT> &u, const HDSA::Vector<RealT> &z) const
  {

    ROL::Ptr<MrHyDE_OptVector> z_rol;
    if (const HDSA_Tpetra_Vector<RealT> *ez = dynamic_cast<const HDSA_Tpetra_Vector<RealT> *>(&z))
    {
      HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> fvec = ez->getVector();
      ROL::Ptr<std::vector<ScalarT>> svec = ROL::makePtr<std::vector<RealT>>(0);
      z_rol = ROL::makePtr<MrHyDE_OptVector>(fvec, svec);
    }
    else if (const Std_Vector<RealT> *ez = dynamic_cast<const Std_Vector<RealT> *>(&z))
    {
      ROL::Ptr<std::vector<ScalarT>> svec = ez->get_std_vec();
      z_rol = ROL::makePtr<MrHyDE_OptVector>(svec);
    }

    params_->updateParams(*z_rol);
    ScalarT val = 0.0;
    solve_->forwardModel(val);

    if (solve_->isTransient)
    {
      Transient_Vector<RealT> &u_trans = dynamic_cast<Transient_Vector<RealT> &>(u);
      int n_t = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
      for (int i = 0; i < n_t; i++)
      {
        HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_vec;
        solve_->postproc->soln[0]->extract(u_vec, i);
        HDSA_Tpetra_Vector<RealT> &eu_i = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(*u_trans[i]);
        eu_i.getVector()->update(1.0, *u_vec, 0.0);
      }
    }
    else
    {
      HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_vec;
      solve_->postproc->soln[0]->extract(u_vec, 0);
      HDSA_Tpetra_Vector<RealT> &eu = dynamic_cast<HDSA_Tpetra_Vector<RealT> &>(u);
      eu.getVector()->update(1.0, *u_vec, 0.0);
    }
  }

  HDSA::Ptr<const HDSA::Vector<RealT>> Extract_State_Component(const HDSA::Vector<RealT> & u, int component_id) const override
  { 

    HDSA::Ptr<const HDSA::Vector<RealT>> u_component;
    int num_states = solve_->varlist[0][0].size();
    if(num_states == 1) 
    {
      u_component = HDSA::makePtrFromRef(u);
    }
    else
    {
      const HDSA_Tpetra_Vector<RealT> &eu = dynamic_cast<const HDSA_Tpetra_Vector<RealT> &>(u);
      HDSA::Ptr<Tpetra::MultiVector<RealT>> eu_tpetra = eu.getVector();
      Teuchos::ArrayRCP<const RealT> u_view = eu_tpetra->get1dView();
      Teuchos::RCP<const Tpetra::Map<LO,GO>> map = eu_tpetra->getMap();

      int num_local_elements = map->getLocalNumElements()/num_states;
      int init_index = map->getMinLocalIndex()/num_states;
      int num_global_element = map->getGlobalNumElements()/num_states;
      Teuchos::Array<GO> component_ids(num_local_elements);
      for (int i = 0; i < num_local_elements; ++i) 
      {
        component_ids[i] = init_index + i;
      }
      Teuchos::RCP<const Tpetra::Map<LO,GO>> component_map = HDSA::makePtr<Tpetra::Map<LO,GO>>(num_global_element, component_ids, 0, solve_->Comm);

      HDSA::Ptr<Tpetra::MultiVector<ScalarT, LO, GO, Node>> tpetra_vec = HDSA::makePtr<Tpetra::MultiVector<ScalarT, LO, GO, Node>>(component_map, 1);
      for (int k = 0; k < num_local_elements; k++)
      {
        tpetra_vec->replaceLocalValue(k, 0, u_view[num_states*k + component_id]);
      }
      u_component = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(tpetra_vec, random_number_generator_);
    }
    return u_component;
  }

  void Set_State_Component(HDSA::Vector<RealT> & u, const HDSA::Vector<RealT> & u_component, int component_id) const override
  { 
    int num_states = solve_->varlist[0][0].size();
    if(num_states == 1) 
    {
      u.set(u_component);
    }
    else
    {
      const HDSA_Tpetra_Vector<RealT> u_tpetra = dynamic_cast<const HDSA_Tpetra_Vector<RealT> &>(u);
      const HDSA_Tpetra_Vector<RealT> u_component_tpetra = dynamic_cast<const HDSA_Tpetra_Vector<RealT> &>(u_component);
      Teuchos::ArrayRCP<const RealT> u_component_view = u_component_tpetra.getVector()->get1dView();
      int local_dim = u_component_tpetra.getVector()->getLocalLength();
      for(int k = 0; k < local_dim; k++)
      {
        u_tpetra.getVector()->replaceLocalValue(num_states*k + component_id, 0, u_component_view[k]);
      }
    }
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_u(void) const override
  {

    int num_time_nodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;

    HDSA::Ptr<HDSA::Vector<RealT>> u_opt;
    if (num_time_nodes > 1)
    {
      std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>> u_tpetra;
      std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> u_hdsa;
      u_tpetra.resize(num_time_nodes);
      u_hdsa.resize(num_time_nodes);
      for (int i = 0; i < num_time_nodes; i++)
      {
        if (opt_solution_exo_file_ != "error")
        {
          u_tpetra[i] = Read_Exodus_Data(opt_solution_exo_file_, true, i + 1);
        }
        else if (opt_solution_txt_file_u_ != "error")
        {
          u_tpetra[i] = Read_Text_Data(opt_solution_txt_file_u_, true, i + 1);
        }
        else
        {
          std::cout << "no valid input file given for Load_Optimal_u" << std::endl;
        }
        u_hdsa[i] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra[i], random_number_generator_);
      }
      u_opt = HDSA::makePtr<Transient_Vector<ScalarT>>(u_hdsa);
    }
    else
    {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_tpetra;
      if (opt_solution_exo_file_.substr(opt_solution_exo_file_.length() - 3, opt_solution_exo_file_.length()) == "exo")
      {
        u_tpetra = Read_Exodus_Data(opt_solution_exo_file_);
      }
      else if (opt_solution_txt_file_u_ != "error")
      {
        u_tpetra = Read_Text_Data(opt_solution_txt_file_u_);
      }
      else
      {
        std::cout << "no valid input file given for Load_Optimal_u" << std::endl;
      }

      u_opt = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra, random_number_generator_);
    }
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_z(void) const override
  {
    HDSA::Ptr<HDSA::Vector<RealT>> z_opt_hdsa;

    if (opt_solution_txt_file_z_ != "error")
    {
      if (params_->getNumParams("discretized") > 0)
      {
        Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra = Read_Text_Data(opt_solution_txt_file_z_, false);
        z_opt_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
      }
      else
      {
        std::vector<RealT> z_vec = Read_Text_Data_std(opt_solution_txt_file_z_);
        z_opt_hdsa = HDSA::makePtr<Std_Vector<RealT>>(z_vec, random_number_generator_);
      }
    }
    else if (opt_solution_exo_file_ != "error")
    {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra = Read_Exodus_Data(opt_solution_exo_file_, false);
      z_opt_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
    }
    else
    {
      std::cout << "no valid input file given for Load_Optimal_z" << std::endl;
    }

    return z_opt_hdsa;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Z_Data(void) const override
  {
    HDSA::Ptr<HDSA::MultiVector<RealT>> Z = HDSA::makePtr<HDSA::MultiVector<RealT>>();

    for (int k = 0; k < num_hifi_; k++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> z_hdsa;
      if (txt_files_z_[k] != "error")
      {
        if (params_->getNumParams("discretized") > 0)
        {
          Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra = Read_Text_Data(txt_files_z_[k], false);
          z_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
        }
        else
        {
          std::vector<RealT> z_vec = Read_Text_Data_std(txt_files_z_[k]);
          z_hdsa = HDSA::makePtr<Std_Vector<RealT>>(z_vec, random_number_generator_);
        }
      }
      else if (hifi_exo_files_[k] != "error")
      {
        Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra = Read_Exodus_Data(hifi_exo_files_[k], false);
        z_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
      }
      else
      {
        std::cout << "no valid input file given for Load_Optimal_Z" << std::endl;
      }
      Z->push_back(z_hdsa);
    }
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_D_Data(void) const override
  {
    int num_time_nodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> u_vecs;
    HDSA::Ptr<HDSA::MultiVector<RealT>> Z = Load_Z_Data();
    for (int k = 0; k < num_hifi_; k++)
    {
      HDSA::Ptr<HDSA::Vector<RealT>> u_k_lofi = Load_Optimal_u()->clone();
      State_Solve(*u_k_lofi, *(*Z)[k]);

      if (num_time_nodes > 1)
      {
        std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>> u_tpetra_hifi;
        std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> u_hdsa_hifi;
        u_tpetra_hifi.resize(num_time_nodes);
        u_hdsa_hifi.resize(num_time_nodes);

        for (int i = 0; i < num_time_nodes; i++)
        {
          if (hifi_exo_files_[k] != "error")
          {
            u_tpetra_hifi[i] = Read_Exodus_Data(hifi_exo_files_[k], true, i + 1);
          }
          else if (hifi_txt_files_u_[k] != "error")
          {
            u_tpetra_hifi[i] = Read_Text_Data(hifi_txt_files_u_[k], true, i + 1);
          }
          else
          {
            std::cout << "no valid input file given for Load_Optimal_D" << std::endl;
          }
          u_hdsa_hifi[i] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra_hifi[i], random_number_generator_);
        }
        HDSA::Ptr<HDSA::Vector<RealT>> u_k_hifi = HDSA::makePtr<Transient_Vector<ScalarT>>(u_hdsa_hifi);
        u_k_hifi->axpy(-1.0, *u_k_lofi);
        u_vecs.push_back(u_k_hifi);
      }
      else
      {
        Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_tpetra_hifi;
        if (hifi_exo_files_[k] != "error")
        {
          u_tpetra_hifi = Read_Exodus_Data(hifi_exo_files_[k]);
        }
        else if (hifi_txt_files_u_[k] != "error")
        {
          u_tpetra_hifi = Read_Text_Data(hifi_txt_files_u_[k]);
        }
        else
        {
          std::cout << "no valid input file given for Load_Optimal_D" << std::endl;
        }
        HDSA::Ptr<HDSA::Vector<RealT>> u_k = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra_hifi, random_number_generator_);
        u_k->axpy(-1.0, *u_k_lofi);
        u_vecs.push_back(u_k);
      }
    }
    HDSA::Ptr<HDSA::MultiVector<RealT>> D = HDSA::makePtr<HDSA::MultiVector<RealT>>(u_vecs);
    return D;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> Read_Exodus_Data(std::string exofile, bool load_state = true, int step = 1) const
  {

    string fname;

    std::vector<std::vector<ScalarT>> nfield_vals, efield_vals;
    std::vector<std::string> block_names, side_names, node_names, nfield_names, efield_names;

    if (comm_->getSize() > 1)
    {
      std::stringstream ssProc, ssPID;
      ssProc << comm_->getSize();
      ssPID << comm_->getRank();
      string strProc = ssProc.str();
      string strPID = ssPID.str();
      // this section may need tweaking if the input exodus mesh is
      // spread across 10's, 100's, or 1000's (etc) of processors
      fname = exofile + "." + strProc + "." + strPID;
    }
    else
    {
      fname = exofile;
    }

    int CPU_word_size, IO_word_size, exoid, exo_error;
    int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;
    char title[MAX_STR_LENGTH + 1];
    float exo_version;
    CPU_word_size = sizeof(ScalarT);
    IO_word_size = 0;
    exoid = ex_open(fname.c_str(), EX_READ, &CPU_word_size, &IO_word_size,
                    &exo_version);
    exo_error = ex_get_init(exoid, title, &num_dim, &num_nods, &num_el,
                            &num_el_blk, &num_ns, &num_ss);
    int id = 1;
    // int step = 1;
    ex_block eblock;
    eblock.id = id;
    eblock.type = EX_ELEM_BLOCK;

    exo_error = ex_get_block_param(exoid, &eblock);

    int num_el_in_blk = eblock.num_entry;
    int num_node_per_el = eblock.num_nodes_per_entry;

    int *connect = new int[num_el_in_blk * num_node_per_el];
    int edgeconn, faceconn;
    exo_error = ex_get_conn(exoid, EX_ELEM_BLOCK, id, connect, &edgeconn, &faceconn);

    int num_node_vars = 0;
    int var_ind;
    exo_error = ex_get_variable_param(exoid, EX_NODAL, &num_node_vars);

    for (int i = 0; i < num_node_vars; i++)
    {
      char varname[MAX_STR_LENGTH + 1];
      ScalarT *var_vals = new ScalarT[num_nods];
      var_ind = i + 1;
      exo_error = ex_get_variable_name(exoid, EX_NODAL, var_ind, varname);
      string vname(varname);
      nfield_names.push_back(vname);
      nfield_vals.push_back(vector<ScalarT>(num_nods));
      exo_error = ex_get_var(exoid, step, EX_NODAL, var_ind, 0, num_nods, var_vals);
      for (int j = 0; j < num_nods; j++)
      {
        nfield_vals[i][j] = var_vals[j];
      }
      delete[] var_vals;
    }

    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec_over;
    if(load_state)
    {
      vec_over = solve_->linalg->getNewOverlappedVector(0);
    } 
    else
    {
      vec_over = solve_->linalg->getNewOverlappedParamVector();
    }
    auto vec_over_kv = vec_over->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
    int index, dindex;

    auto dev_offsets = solve_->assembler->wkset[0]->offsets;
    auto offsets = Kokkos::create_mirror_view(dev_offsets);
    Kokkos::deep_copy(offsets, dev_offsets);

    vector<string> blockNames = solve_->mesh->getBlockNames();

    std::vector<std::string> state_vars = solve_->varlist[0][0]; // [0][0] accesses set=0 and block=0
    std::vector<std::string> opt_vars = solve_->params->getParamsNames(4);

    std::vector<int> n_list;
    if (load_state)
    {
      for (int j = 0; j < nfield_names.size(); j++)
      {
        auto it = std::find(state_vars.begin(),state_vars.end(),nfield_names[j]);
        if (it != state_vars.end())
        {
          n_list.push_back(j);
        }
      }
    }
    else
    {
      for (int j = 0; j < nfield_names.size(); j++)
      {
        auto it = std::find(opt_vars.begin(),opt_vars.end(),nfield_names[j]);
        if (it != opt_vars.end())
        {
          n_list.push_back(j);
        }
      }
    }

    int index_normalization = 1;
    if (!load_state)
    {
      index_normalization = state_vars.size(); 
      if(opt_vars.size() > 1)
      {
        std::cout << "Data loading with multiple optimization variables defined on the mesh is currently not supported" << std::endl;
      }
    }

    for(int k = 0; k < n_list.size(); k++)
    {
      int n = n_list[k];
      int n_var = n;
      if(!load_state)
      {
        n_var = 0;
      }
      for (int block = 0; block < blockNames.size(); block++)
      {
        int e = 0;
        for (int grp = 0; grp < solve_->assembler->groups[block].size(); ++grp)
        {
          auto LIDs = solve_->assembler->groups[block][grp]->LIDs_host;
          auto nDOF = solve_->assembler->groups[block][grp]->group_data->num_dof_host;
          for (int p = 0; p < solve_->assembler->groups[block][grp]->numElem; p++)
          {
            for (int i = 0; i < nDOF(n_var); i++)
            {
              index = LIDs[0](p, offsets(n_var, i))/index_normalization;
              dindex = connect[e * num_node_per_el + i] - 1;
              vec_over_kv(index, 0) = nfield_vals[n][dindex];
            }
            e += 1;
          }
        }
      }
    }
    delete[] connect;

    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec;
    if(load_state)
    {
      vec = solve_->linalg->getNewVector(0);
      solve_->linalg->exportVectorFromOverlappedReplace(0, vec, vec_over);
    } 
    else
    {
      vec = solve_->linalg->getNewParamVector();
      solve_->linalg->exportParamVectorFromOverlappedReplace(vec, vec_over);
    }

    exo_error = ex_close(exoid);
    if (exo_error != 0)
    {
      std::cout << "Exodus reader error" << std::endl;
    }
    return vec;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Read_Spatial_Node_Data() const
  {
    HDSA::Ptr<HDSA::MultiVector<RealT>> spatial_nodes;
    if (opt_solution_exo_file_ != "error")
    {
      spatial_nodes = Read_Exodus_Spatial_Node_Data(opt_solution_exo_file_);
    }
    else
    {
      std::cout << "Error: Read_Spatial_Node_Data is currently only supported for Exodus file reading" << std::endl;

      // Default to assume 1D on the domain [0,1] so that test problems will run
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
      int dim = vec->getGlobalLength();
      for (int k = 0; k < dim; k++)
      {
        RealT val = double(k) / double(dim - 1);
        vec->replaceGlobalValue(k, 0, val);
      }
      std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> coord_vecs;
      coord_vecs.resize(1);
      coord_vecs[0] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(vec, random_number_generator_);
      spatial_nodes = HDSA::makePtr<HDSA::MultiVector<RealT>>(coord_vecs);
    }
    return spatial_nodes;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Read_Exodus_Spatial_Node_Data(std::string exofile) const
  {

    string fname;
    std::vector<std::string> block_names, side_names, node_names;

    if (comm_->getSize() > 1)
    {
      std::stringstream ssProc, ssPID;
      ssProc << comm_->getSize();
      ssPID << comm_->getRank();
      string strProc = ssProc.str();
      string strPID = ssPID.str();
      // this section may need tweaking if the input exodus mesh is
      // spread across 10's, 100's, or 1000's (etc) of processors
      fname = exofile + "." + strProc + "." + strPID;
    }
    else
    {
      fname = exofile;
    }

    int CPU_word_size, IO_word_size, exoid, exo_error;
    int num_dim, num_nods, num_el, num_el_blk, num_ns, num_ss;
    char title[MAX_STR_LENGTH + 1];
    float exo_version;
    CPU_word_size = sizeof(ScalarT);
    IO_word_size = 0;
    exoid = ex_open(fname.c_str(), EX_READ, &CPU_word_size, &IO_word_size,
                    &exo_version);
    exo_error = ex_get_init(exoid, title, &num_dim, &num_nods, &num_el,
                            &num_el_blk, &num_ns, &num_ss);
    int id = 1;
    // int step = 1;
    ex_block eblock;
    eblock.id = id;
    eblock.type = EX_ELEM_BLOCK;

    exo_error = ex_get_block_param(exoid, &eblock);

    int num_el_in_blk = eblock.num_entry;
    int num_node_per_el = eblock.num_nodes_per_entry;

    int *connect = new int[num_el_in_blk * num_node_per_el];
    int edgeconn, faceconn;
    exo_error = ex_get_conn(exoid, EX_ELEM_BLOCK, id, connect, &edgeconn, &faceconn);

    int num_node_vars = 0;
    exo_error = ex_get_variable_param(exoid, EX_NODAL, &num_node_vars);

    ScalarT *x_coords = new double[num_nods];
    ScalarT *y_coords = new double[num_nods];
    ScalarT *z_coords = new double[num_nods];
    ex_get_coord(exoid, x_coords, y_coords, z_coords);
    std::vector<std::vector<ScalarT>> spatial_coords;
    spatial_coords.resize(num_dim);
    for (int i = 0; i < num_dim; i++)
    {
      spatial_coords[i].resize(num_nods);
      for (int j = 0; j < num_nods; j++)
      {
        if (i == 0)
        {
          spatial_coords[i][j] = x_coords[j];
        }
        if (i == 1)
        {
          spatial_coords[i][j] = y_coords[j];
        }
        if (i == 2)
        {
          spatial_coords[i][j] = z_coords[j];
        }
      }
    }

    std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> coord_vecs;
    coord_vecs.resize(num_dim);

    for (int spatial_id = 0; spatial_id < num_dim; spatial_id++)
    {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec_over = Teuchos::rcp(new Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>(solve_->linalg->overlapped_map[0], 1));
      auto vec_over_kv = vec_over->template getLocalView<HostDevice>(Tpetra::Access::ReadWrite);
      int index, dindex;

      auto dev_offsets = solve_->assembler->wkset[0]->offsets;
      auto offsets = Kokkos::create_mirror_view(dev_offsets);
      Kokkos::deep_copy(offsets, dev_offsets);

      vector<string> blockNames = solve_->mesh->getBlockNames();
      for (int block = 0; block < blockNames.size(); block++)
      {
        for (int grp = 0; grp < solve_->assembler->groups[block].size(); ++grp)
        {
          auto LIDs = solve_->assembler->groups[block][grp]->LIDs_host;
          auto nDOF = solve_->assembler->groups[block][grp]->group_data->num_dof_host;
          for (int p = 0; p < solve_->assembler->groups[block][grp]->numElem; p++)
          {
            for (int i = 0; i < nDOF(0); i++)
            {
              index = LIDs[0](p, offsets(0, i));
              dindex = connect[grp * num_node_per_el + i] - 1;
              vec_over_kv(index, 0) = spatial_coords[spatial_id][dindex];
            }
          }
        }
      }
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
      solve_->linalg->exportVectorFromOverlappedReplace(0, vec, vec_over);
      coord_vecs[spatial_id] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(vec, random_number_generator_);
    }

    delete[] connect;
    exo_error = ex_close(exoid);
    if (exo_error != 0)
    {
      std::cout << "Exodus reader error" << std::endl;
    }

    HDSA::Ptr<HDSA::MultiVector<RealT>> coords = HDSA::makePtr<HDSA::MultiVector<RealT>>(coord_vecs);
    return coords;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> Read_Text_Data(std::string txtfile, bool load_state = true, int step = 1) const
  {
    RealT val = 0.0;
    int num_time_nodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
    int spatial_dim = vec->getGlobalLength();

    // read in data
    std::ifstream in(txtfile);
    if (in)
    {
      if (load_state)
      {
        for (int i = 0; i < num_time_nodes; i++)
        {
          for (int j = 0; j < spatial_dim; j++)
          {
            in >> val;
            if (i == step - 1)
            {
              vec->replaceGlobalValue(j, 0, val);
            }
          }
        }
      }
      else
      {
        for (int j = 0; j < spatial_dim; j++)
        {
          in >> val;
          vec->replaceGlobalValue(j, 0, val);
        }
      }
    }
    else
    {
      std::cout << "Error loading the data from " << txtfile << std::endl;
    }
    return vec;
  }

  std::vector<RealT> Read_Text_Data_std(std::string txtfile) const
  {
    RealT val = 0.0;
    int dim = params_->getNumParams("active");
    std::vector<RealT> vec = std::vector<RealT>(dim, 0.0);

    // read in data
    std::ifstream in(txtfile);
    if (in)
    {
      for (int j = 0; j < dim; j++)
      {
        in >> val;
        vec[j] = val;
      }
    }
    else
    {
      std::cout << "Error loading the data from " << txtfile << std::endl;
    }
    return vec;
  }
};
#endif
