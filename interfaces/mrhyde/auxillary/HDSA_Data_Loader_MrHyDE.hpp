/***********************************************************************
 HdsaLib - A library for Hyper-differential Sensitivity Analysis
 
 Questions? Contact Joseph Hart (joshart@sandia.gov)
************************************************************************/

#ifndef HDSA_DATA_LOADER_MRHYDE_HPP
#define HDSA_DATA_LOADER_MRHYDE_HPP

#include "exodusII.h"
#include "HDSA_MD_Data_Interface.hpp"
#include "HDSA_Tpetra_Vector.hpp"
#include "HDSA_Std_Vector.hpp"

template <class RealT,
          class LO = Tpetra::Map<>::local_ordinal_type,
          class GO = Tpetra::Map<>::global_ordinal_type,
          class Node = Tpetra::Map<>::node_type>
class Data_Loader_MrHyDE
{

private:
  Teuchos::RCP<Teuchos::MpiComm<int>> comm_;
  Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> solve_;
  HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> params_;
  const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> random_number_generator_;

public:
  Data_Loader_MrHyDE(Teuchos::RCP<Teuchos::MpiComm<int>> &comm, Teuchos::RCP<MrHyDE::SolverManager<SolverNode>> &solve, HDSA::Ptr<MrHyDE::ParameterManager<SolverNode>> &params, const HDSA::Ptr<HDSA::Random_Number_Generator<RealT>> &random_number_generator) : comm_(comm), solve_(solve), params_(params), random_number_generator_(random_number_generator)
  {
  }

  virtual ~Data_Loader_MrHyDE()
  {
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
    if (load_state)
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
        auto it = std::find(state_vars.begin(), state_vars.end(), nfield_names[j]);
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
        auto it = std::find(opt_vars.begin(), opt_vars.end(), nfield_names[j]);
        if (it != opt_vars.end())
        {
          n_list.push_back(j);
        }
      }
    }

    int index_Normalization = 1;
    if (!load_state)
    {
      index_Normalization = state_vars.size();
      if (opt_vars.size() > 1)
      {
        HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                                "Error in HDSA_Data_Loader_MrHyDE: Data loading with multiple optimization variables defined on the mesh is currently not supported" << std::endl);
      }
    }

    for (int k = 0; k < n_list.size(); k++)
    {
      int n = n_list[k];
      int n_var = n;
      if (!load_state)
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
              index = LIDs[0](p, offsets(n_var, i)) / index_Normalization;
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
    if (load_state)
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
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA_Data_Loader_MrHyDE: Exodus reader failure" << std::endl);
    }
    return vec;
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
        int e = 0;
        for (int grp = 0; grp < solve_->assembler->groups[block].size(); ++grp)
        {
          auto LIDs = solve_->assembler->groups[block][grp]->LIDs_host;
          auto nDOF = solve_->assembler->groups[block][grp]->group_data->num_dof_host;
          for (int p = 0; p < solve_->assembler->groups[block][grp]->numElem; p++)
          {
            for (int i = 0; i < nDOF(0); i++)
            {
              index = LIDs[0](p, offsets(0, i));
              dindex = connect[e * num_node_per_el + i] - 1;
              vec_over_kv(index, 0) = spatial_coords[spatial_id][dindex];
            }
            e += 1;
          }
        }
      }
      Teuchos::RCP<Tpetra::MultiVector<ScalarT, LO, GO, SolverNode>> vec = solve_->linalg->getNewVector(0);
      solve_->linalg->exportVectorFromOverlappedReplace(0, vec, vec_over);
      coord_vecs[spatial_id] = HDSA::makePtr<HDSA::Tpetra_Vector<RealT>>(vec, random_number_generator_);
    }

    delete[] connect;
    exo_error = ex_close(exoid);
    if (exo_error != 0)
    {
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA_Data_Loader_MrHyDE: Exodus reader failure" << std::endl);
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
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA_Data_Loader_MrHyDE: Failed to load data from " << txtfile << std::endl);
    }
    return vec;
  }

  std::vector<std::vector<RealT>> Read_Text_Data_Dynamic_std(std::string txtfile) const
  {
    int num_params = params_->getNumParams("active");
    int num_time_steps = params_->getCurrentVector().dimension() / num_params; // This assumes that all active parameters are dynamic
    ScalarT val = 0.0;
    std::ifstream in(txtfile);
    std::vector<std::vector<ScalarT>> vec;
    vec.resize(num_time_steps);
    if (in)
    {
      for (int i = 0; i < num_time_steps; i++)
      {
        vec[i].resize(num_params);
        // read the elements in the file into a vector
        for (int j = 0; j < num_params; j++)
        {
          in >> val;
          vec[i][j] = val;
        }
        params_->dynamic_timeindex = i;
        params_->updateParams(vec[i], "active");
      }
      params_->dynamic_timeindex = 0;
    }
    else
    {
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA_Data_Loader_MrHyDE: Failed to load data from " << txtfile << std::endl);
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
      HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
                              "Error in HDSA_Data_Loader_MrHyDE: Failed to load data from " << txtfile << std::endl);
    }
    return vec;
  }
};
#endif
