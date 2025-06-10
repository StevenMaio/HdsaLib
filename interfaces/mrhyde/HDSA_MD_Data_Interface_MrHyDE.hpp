#ifndef HDSA_MD_DATA_INTERFACE_MRHYDE_HPP
#define HDSA_MD_DATA_INTERFACE_MRHYDE_HPP

#include "exodusII.h"
#include "preferences.hpp"

template <class RealT,
          class LO = Tpetra::Map<>::local_ordinal_type,
          class GO = Tpetra::Map<>::global_ordinal_type,
          class Node = Tpetra::Map<>::node_type>
class MD_Data_Interface_MrHyDE : public HDSA::MD_Data_Interface<RealT>
{

private:
  Teuchos::RCP<Teuchos::MpiComm<int>> comm_;
  Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solve_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;
  Teuchos::ParameterList data_load_list_;
  int num_hifi_;

  std::string opt_var_physics_;
  std::string opt_solution_exo_file_;
  std::vector<std::string> hifi_exo_files_;
  std::vector<std::string> lofi_exo_files_;

  std::string opt_solution_txt_file_u_;
  std::string opt_solution_txt_file_z_;
  std::vector<std::string> hifi_txt_files_u_;
  std::vector<std::string> lofi_txt_files_u_;
  std::vector<std::string> txt_files_z_;

public:
  MD_Data_Interface_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int>> &comm, Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solve, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator, Teuchos::ParameterList &data_load_list) : comm_(comm), solve_(solve), random_number_generator_(random_number_generator), data_load_list_(data_load_list)
  {
    num_hifi_ = data_load_list_.get<int>("NumHifi", 1);

    opt_var_physics_ = data_load_list_.get<std::string>("OptVariablePhysics", "error");
    opt_solution_exo_file_ = data_load_list_.get<std::string>("OptimalSolutionExoFile", "error");

    opt_solution_txt_file_u_ = data_load_list_.get<std::string>("OptimalSolutionTxtFileU", "error");
    opt_solution_txt_file_z_ = data_load_list_.get<std::string>("OptimalSolutionTxtFileZ", "error");

    hifi_exo_files_.resize(num_hifi_);
    lofi_exo_files_.resize(num_hifi_);

    hifi_txt_files_u_.resize(num_hifi_);
    lofi_txt_files_u_.resize(num_hifi_);
    txt_files_z_.resize(num_hifi_);
    for (int k = 0; k < num_hifi_; k++)
    {
      hifi_exo_files_[k] = data_load_list_.get<std::string>("HifiExoFile" + std::to_string(k + 1), "error");
      lofi_exo_files_[k] = data_load_list_.get<std::string>("LofiExoFile" + std::to_string(k + 1), "error");

      hifi_txt_files_u_[k] = data_load_list_.get<std::string>("HifiTxtFileU" + std::to_string(k + 1), "error");
      lofi_txt_files_u_[k] = data_load_list_.get<std::string>("LofiTxtFileU" + std::to_string(k + 1), "error");
      txt_files_z_[k] = data_load_list_.get<std::string>("TxtFileZ" + std::to_string(k + 1), "error");
    }

    if (opt_var_physics_ == "error")
    {
      std::cout << "Error: specify OptVariablePhysics" << std::endl;
    }
    if (opt_solution_exo_file_ == "error")
    {
      std::cout << "Error: specify OptimalSolutionExoFile" << std::endl;
    }
    for (int k = 0; k < num_hifi_; k++)
    {
      if (hifi_exo_files_[k] == "error")
      {
        std::cout << "Error: specify HifiExoFile" + std::to_string(k + 1) << std::endl;
      }
      if (lofi_exo_files_[k] == "error")
      {
        std::cout << "Error: specify fiExoFile" + std::to_string(k + 1) << std::endl;
      }
    }
  }

  virtual ~MD_Data_Interface_MrHyDE()
  {
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_u(void) const
  {

    int numtimenodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;

    HDSA::Ptr<HDSA::Vector<RealT>> u_opt;
    if (numtimenodes > 1)
    {
      std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>> u_tpetra;
      std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> u_hdsa;
      u_tpetra.resize(numtimenodes);
      u_hdsa.resize(numtimenodes);
      for (int i = 0; i < numtimenodes; i++)
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
      else
      {
        u_tpetra = Read_Text_Data(opt_solution_txt_file_u_);
      }

      HDSA::Ptr<State_Vector_MrHyDE<RealT>> u_opt_mrhyde = HDSA::makePtr<State_Vector_MrHyDE<RealT>>(solve_, random_number_generator_);
      u_opt_mrhyde->mrhyde_state_vec[0][0]->update(1.0, *u_tpetra, 0.0);
      u_opt = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_opt_mrhyde->mrhyde_state_vec[0][0], random_number_generator_);
    }
    return u_opt;
  }

  HDSA::Ptr<HDSA::Vector<RealT>> Load_Optimal_z(void) const
  {
    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra;
    if (opt_solution_exo_file_ != "error")
    {
      z_tpetra = Read_Exodus_Data(opt_solution_exo_file_, false);
    }
    else if (opt_solution_txt_file_z_ != "error")
    {
      z_tpetra = Read_Text_Data(opt_solution_txt_file_z_, false);
    }
    else
    {
      std::cout << "no valid input file given for Load_Optimal_z" << std::endl;
    }
    HDSA::Ptr<HDSA::Vector<RealT>> z_opt_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
    return z_opt_hdsa;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_Z_Data(void) const
  {
    HDSA::Ptr<HDSA::MultiVector<RealT>> Z = HDSA::makePtr<HDSA::MultiVector<RealT>>();

    for (int k = 0; k < num_hifi_; k++)
    {
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> z_tpetra;
      if (hifi_exo_files_[k] != "error")
      {
        z_tpetra = Read_Exodus_Data(hifi_exo_files_[k], false);
      }
      else if (txt_files_z_[k] != "error")
      {
        z_tpetra = Read_Text_Data(txt_files_z_[k], false);
      }
      else
      {
        std::cout << "no valid input file given for Load_Optimal_Z" << std::endl;
      }
      HDSA::Ptr<HDSA::Vector<RealT>> z_hdsa = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(z_tpetra, random_number_generator_);
      Z->push_back(z_hdsa);
    }
    return Z;
  }

  HDSA::Ptr<HDSA::MultiVector<RealT>> Load_D_Data(void) const
  {
    int numtimenodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
    std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> u_vecs;
    for (int k = 0; k < num_hifi_; k++)
    {
      if (numtimenodes > 1)
      {
        std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>> u_tpetra_lofi;
        std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> u_hdsa_lofi;
        u_tpetra_lofi.resize(numtimenodes);
        u_hdsa_lofi.resize(numtimenodes);

        std::vector<Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>>> u_tpetra_hifi;
        std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> u_hdsa_hifi;
        u_tpetra_hifi.resize(numtimenodes);
        u_hdsa_hifi.resize(numtimenodes);

        for (int i = 0; i < numtimenodes; i++)
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

          if (lofi_exo_files_[k] != "error")
          {
            u_tpetra_lofi[i] = Read_Exodus_Data(lofi_exo_files_[k], true, i + 1);
          }
          else if (lofi_txt_files_u_[k] != "error")
          {
            u_tpetra_lofi[i] = Read_Text_Data(lofi_txt_files_u_[k], true, i + 1);
          }
          else
          {
            std::cout << "no valid input file given for Load_Optimal_D" << std::endl;
          }

          u_hdsa_hifi[i] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra_hifi[i], random_number_generator_);
          u_hdsa_lofi[i] = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra_lofi[i], random_number_generator_);
        }
        HDSA::Ptr<HDSA::Vector<RealT>> u_k_hifi = HDSA::makePtr<Transient_Vector<ScalarT>>(u_hdsa_hifi);
        HDSA::Ptr<HDSA::Vector<RealT>> u_k_lofi = HDSA::makePtr<Transient_Vector<ScalarT>>(u_hdsa_lofi);
        u_k_hifi->axpy(-1.0, *u_k_lofi);
        u_vecs.push_back(u_k_hifi);
      }
      else
      {
        Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_tpetra_hifi;
        Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> u_tpetra_lofi;
        if (hifi_exo_files_[k] != "error")
        {
          u_tpetra_hifi = Read_Exodus_Data(hifi_exo_files_[k]);
          u_tpetra_lofi = Read_Exodus_Data(lofi_exo_files_[k]);
        }
        else if (hifi_txt_files_u_[k] != "error")
        {
          u_tpetra_hifi = Read_Text_Data(hifi_txt_files_u_[k]);
          u_tpetra_lofi = Read_Text_Data(lofi_txt_files_u_[k]);
        }
        else
        {
          std::cout << "no valid input file given for Load_Optimal_D" << std::endl;
        }

        // HDSA::Ptr<State_Vector_MrHyDE<RealT> > u_vec_k = HDSA::makePtr<State_Vector_MrHyDE<RealT> >(solve_,random_number_generator_);
        // u_vec_k->mrhyde_state_vec[0][0]->update(1.0,*u_tpetra_hifi,0.0);
        // u_vec_k->mrhyde_state_vec[0][0]->update(-1.0,*u_tpetra_lofi,1.0);
        // HDSA::Ptr<HDSA::Vector<RealT> > u_k = HDSA::makePtr<HDSA_Tpetra_Vector<RealT> >(u_vec_k->mrhyde_state_vec[0][0],random_number_generator_);

        u_tpetra_hifi->update(-1.0, *u_tpetra_lofi, 1.0);
        HDSA::Ptr<HDSA::Vector<RealT>> u_k = HDSA::makePtr<HDSA_Tpetra_Vector<RealT>>(u_tpetra_hifi, random_number_generator_);
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

        int bound_nv = solve_->physics->getVarList()[0][block].size();
        if (!load_state)
        {
          bound_nv = solve_->params->getParamsNames(4).size();
        }

        for (int nv = 0; nv < bound_nv; nv++)
        {
          std::string var;
          if (load_state)
          {
            var = solve_->physics->getVarList()[0][block][nv];
          }
          else
          {
            var = solve_->params->getParamsNames(4)[nv];
          }
          int n = -1;

          for (int j = 0; j < nfield_names.size(); j++)
          {
            if (nfield_names[j] == var)
            {
              n = j;
            }
          }
          int n_var = n;
          if (!load_state)
          {
            for (int j = 0; j < nfield_names.size(); j++)
            {
              if (nfield_names[j] == opt_var_physics_)
              {
                n_var = j;
              }
            }
          }
          for (int p = 0; p < solve_->assembler->groups[block][grp]->numElem; p++)
          {
            for (int i = 0; i < nDOF(n_var); i++)
            {
              index = LIDs[0](p, offsets(n_var, i));
              dindex = connect[grp * num_node_per_el + i] - 1;
              vec_over_kv(index, 0) = nfield_vals[n][dindex];
            }
          }
        }
      }
    }
    delete[] connect;

    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
    solve_->linalg->exportVectorFromOverlappedReplace(0, vec, vec_over);

    exo_error = ex_close(exoid);
    return vec;
  }

  Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> Read_Text_Data(std::string txtfile, bool load_state = true, int step = 1) const
  {
    RealT val = 0.0;
    int numtimenodes = solve_->settings->sublist("Solver").get<int>("number of steps", 0) + 1;
    Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
    int spatialdim = vec->getGlobalLength();

    // read in data
    std::ifstream in(txtfile);
    if (in)
    {
      if (load_state)
      {
        for (int i = 0; i < numtimenodes; i++)
        {
          for (int j = 0; j < spatialdim; j++)
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
        for (int j = 0; j < spatialdim; j++)
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
};
#endif
