#ifndef HDSA_WRITE_OUTPUT_MRHYDE_HPP
#define HDSA_WRITE_OUTPUT_MRHYDE_HPP

#include <filesystem>
#include "HDSA_MD_u_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Multi_State_u_Hyperparameter_Interface.hpp"
#include "HDSA_MD_z_Hyperparameter_Interface.hpp"
#include "HDSA_MD_Posterior_Vectors.hpp"
#include "HDSA_Ensemble_Vector.hpp"

template <class RealT>
class Write_Output_MrHyDE
{

private:
    HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> postproc_;
    HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> solver_;
    bool write_exo_;
    std::string output_dir_name_;

public:
    Write_Output_MrHyDE(const HDSA::Ptr<MrHyDE::PostprocessManager<SolverNode>> &postproc, const HDSA::Ptr<MrHyDE::SolverManager<SolverNode>> &solver, bool write_exo) : postproc_(postproc), solver_(solver), write_exo_(write_exo)
    {
        if (write_exo_)
        {
            postproc_->write_optimization_solution = false;
        }

        output_dir_name_ = "hdsa_output";
        if (solver_->Comm->getRank() == 0)
        {
            try
            {
                std::filesystem::create_directory(output_dir_name_);
            }
            catch (const std::exception &e)
            {
            }
        }
        solver_->Comm->barrier();
    }

    virtual ~Write_Output_MrHyDE()
    {
    }

    void Write_Hyperparameters(const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface) const
    {
        std::ofstream outfile(output_dir_name_ + "/hyperparameters.txt");
        if (outfile.is_open())
        {
            // u hyperparameters
            if (u_hyperparam_interface->Is_Multi_State_Interface())
            {
                HDSA::MD_Multi_State_u_Hyperparameter_Interface<RealT> *u_hyperparam_interface_multistate = dynamic_cast<HDSA::MD_Multi_State_u_Hyperparameter_Interface<RealT> *>(&(*u_hyperparam_interface));
                int num_states = u_hyperparam_interface_multistate->Get_Number_of_States();
                for (int k = 0; k < num_states; k++)
                {
                    outfile << "Component " << k << std::endl;
                    outfile << "alpha_u: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_alpha_u() << std::endl;
                    outfile << "beta_u: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_beta_u() << std::endl;
                    if (u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Is_Transient())
                    {
                        outfile << "alpha_t:";
                        for (int k = 0; k < u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_alpha_t().size(); k++)
                        {
                            outfile << " " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_alpha_t()[k];
                        }
                        outfile << " " << std::endl;
                        outfile << "beta_t: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_beta_t() << std::endl;
                    }
                    outfile << "alpha_d: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_alpha_d() << std::endl;
                    outfile << "gsvd_num_sing_vals: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_gsvd_num_sing_vals() << std::endl;
                    outfile << "gsvd_oversampling: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_gsvd_oversampling() << std::endl;
                    outfile << "gsvd_num_subspace_iter: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Get_gsvd_num_subspace_iter() << std::endl;
                    outfile << "Center Data: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Center_Data() << std::endl;
                    if (u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Is_Transient())
                    {
                        outfile << "Adapt_Time_Variance: " << u_hyperparam_interface_multistate->Get_Hyperparameter_Interface_k(k)->Adapt_Time_Variance() << std::endl;
                    }
                    outfile << " " << std::endl;
                }
            }
            else
            {
                outfile << "alpha_u: " << u_hyperparam_interface->Get_alpha_u() << std::endl;
                outfile << "beta_u: " << u_hyperparam_interface->Get_beta_u() << std::endl;
                if (u_hyperparam_interface->Is_Transient())
                {
                    outfile << "alpha_t:";
                    for (int k = 0; k < u_hyperparam_interface->Get_alpha_t().size(); k++)
                    {
                        outfile << " " << u_hyperparam_interface->Get_alpha_t()[k];
                    }
                    outfile << " " << std::endl;
                    outfile << "beta_t: " << u_hyperparam_interface->Get_beta_t() << std::endl;
                }
                outfile << "alpha_d: " << u_hyperparam_interface->Get_alpha_d() << std::endl;
                outfile << "gsvd_num_sing_vals: " << u_hyperparam_interface->Get_gsvd_num_sing_vals() << std::endl;
                outfile << "gsvd_oversampling: " << u_hyperparam_interface->Get_gsvd_oversampling() << std::endl;
                outfile << "gsvd_num_subspace_iter: " << u_hyperparam_interface->Get_gsvd_num_subspace_iter() << std::endl;
                outfile << "Center Data: " << u_hyperparam_interface->Center_Data() << std::endl;
                if (u_hyperparam_interface->Is_Transient())
                {
                    outfile << "Adapt_Time_Variance: " << u_hyperparam_interface->Adapt_Time_Variance() << std::endl;
                }
            }

            outfile << " " << std::endl;
            // z hyperparameters
            outfile << "alpha_z: " << z_hyperparam_interface->Get_alpha_z() << std::endl;
            outfile << "beta_z: " << z_hyperparam_interface->Get_beta_z() << std::endl;
            outfile << "beta_t: " << z_hyperparam_interface->Get_beta_t() << std::endl;

            // Close the file
            outfile.close();
        }
        else
        {
            std::cout << "Error: Unable to write hyperparameters to a file" << std::endl;
        }
    }

    void Write_Prior_Discrepancy_Samples(HDSA::Ptr<HDSA::MultiVector<ScalarT>> &prior_delta_z_opt, std::vector<HDSA::Ptr<HDSA::Vector<ScalarT>>> &prior_z_pert, std::vector<HDSA::Ptr<HDSA::MultiVector<ScalarT>>> &prior_delta_z_pert) const
    {
        std::filesystem::create_directory(output_dir_name_ + "/prior");

        std::string name = output_dir_name_ + "/prior/prior_delta_z_opt";
        std::filesystem::create_directory(name);
        Write_to_File(prior_delta_z_opt, name, true);

        int N = prior_z_pert.size();
        for (int k = 0; k < N; k++)
        {
            name = output_dir_name_ + "/prior/prior_z_pert_" + std::to_string(k + 1);
            Write_to_File(prior_z_pert[k], name, false);

            name = output_dir_name_ + "/prior/prior_delta_z_pert_" + std::to_string(k + 1);
            std::filesystem::create_directory(name);
            Write_to_File(prior_delta_z_pert[k], name, true);
        }
    }

    void Write_Prior_Discrepancy_Time_Evolution(std::vector<std::vector<std::vector<RealT>>> &prior_delta_z_opt_time_evol, std::vector<std::vector<RealT>> &prior_discrep_data_time_evol) const
    {
        std::string name = output_dir_name_ + "/prior/summary_statistics";
        std::filesystem::create_directory(name);

        int num_samps = prior_delta_z_opt_time_evol.size();
        int num_time = prior_delta_z_opt_time_evol[0].size();
        int num_states = prior_delta_z_opt_time_evol[0][0].size();

        for (int j = 0; j < num_states; j++)
        {
            std::string name_j = name + "/prior_delta_z_opt_time_evol_" + std::to_string(j + 1) + ".txt";
            std::ofstream fout;
            fout.open(name_j);
            for (int i = 0; i < num_samps; i++)
            {
                for (int k = 0; k < num_time; k++)
                {
                    fout << std::setprecision(8) << prior_delta_z_opt_time_evol[i][k][j] << "  ";
                }
                fout << " " << std::endl;
            }
            fout.close();
        }

        for (int j = 0; j < num_states; j++)
        {
            std::string name_j = name + "/prior_discrep_data_time_evol_" + std::to_string(j + 1) + ".txt";
            std::ofstream fout;
            fout.open(name_j);
            for (int k = 0; k < num_time; k++)
            {
                fout << std::setprecision(8) << prior_discrep_data_time_evol[k][j] << "  ";
            }
            fout << " " << std::endl;
            fout.close();
        }
    }

    void Write_Posterior_Discrepancy_Samples(std::vector<HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT>>> post_delta) const
    {
        std::filesystem::create_directory(output_dir_name_ + "/posterior");
        int N = post_delta.size();
        for (int k = 0; k < N; k++)
        {
            std::string name = output_dir_name_ + "/posterior/posterior_delta_z_" + std::to_string(k + 1);
            std::filesystem::create_directory(name);
            std::string filename = name + "/posterior_mean";
            Write_to_File(post_delta[k]->mean, filename, true);
            name = name + "/posterior_samples";
            std::filesystem::create_directory(name);
            Write_to_File(post_delta[k]->samples, name, true);
        }
    }

    void Write_Hessian_Eigenvalues(HDSA::Ptr<HDSA::Dense_Matrix<RealT>> &evals) const
    {
        std::string name = output_dir_name_ + "/hessian_evals.txt";
        std::ofstream fout;
        fout.open(name);
        int num_evals = evals->numRows();
        for (int i = 0; i < num_evals; i++)
        {
            fout << std::setprecision(8) << (*evals)(i, 0) << "  ";
        }
        fout.close();
    }

    void Write_Optimal_Solution_Update(HDSA::Ptr<HDSA::MD_Posterior_Vectors<ScalarT>> &posterior_update_samples) const
    {
        std::string name = output_dir_name_ + "/posterior";
        std::filesystem::create_directory(name);
        name = name + "/z_update";
        std::filesystem::create_directory(name);
        std::string filename = name + "/mean";
        Write_to_File(posterior_update_samples->mean, filename, false);
        name = name + "/posterior_samples";
        std::filesystem::create_directory(name);
        Write_to_File(posterior_update_samples->samples, name, false);
    }

    void Write_Optimal_Solution_Update(HDSA::Ptr<HDSA::Vector<ScalarT>> &posterior_update_mean) const
    {
        std::string name = output_dir_name_ + "/posterior";
        std::filesystem::create_directory(name);
        name = name + "/z_update";
        std::filesystem::create_directory(name);
        name = name + "/mean";
        Write_to_File(posterior_update_mean, name, false);
    }

    void Write_to_File(const HDSA::Ptr<HDSA::MultiVector<RealT>> &vec, std::string &name, bool is_state) const
    {
        int num_vecs = vec->Number_of_Vectors();
        for (int k = 0; k < num_vecs; k++)
        {
            HDSA::Ptr<HDSA::Vector<RealT>> vec_k = (*vec)[k];
            std::string filename = name + "/Vector_" + std::to_string(k + 1);
            Write_to_File(vec_k, filename, is_state);
        }
    }

    void Write_to_File(const HDSA::Ptr<HDSA::Vector<RealT>> &vec, std::string &filename, bool is_state) const
    {
        if (write_exo_)
        {
            std::string name = filename + ".exo";
            postproc_->write_solution = true;
            postproc_->exodus_filename = name;

            if (is_state)
            {
                if (HDSA::Ensemble_Vector<RealT> *evec = dynamic_cast<HDSA::Ensemble_Vector<RealT> *>(&(*vec)))
                {
                    for (int s = 0; s < evec->Number_of_Vectors(); s++)
                    {
                        std::string name_s = filename + "_ens_" + std::to_string(s + 1) + ".exo";
                        Write_to_File((*evec)[s], name_s, is_state);
                    }
                }
                else
                {
                    std::vector<std::string> discretized_param_names = postproc_->params->discretized_param_names;
                    postproc_->params->discretized_param_names.clear();
                    postproc_->setNewExodusFile(name);
                    std::vector<HDSA::Ptr<Tpetra::MultiVector<RealT>>> sol;
                    sol.resize(1);
                    RealT current_time = 0.0;
                    if (HDSA::Transient_Vector<RealT> *evec = dynamic_cast<HDSA::Transient_Vector<RealT> *>(&(*vec)))
                    {
                        int n_t = evec->Get_n_t();
                        for (int i = 0; i < n_t; i++)
                        {
                            HDSA::Tpetra_Vector<RealT> &evec_i = dynamic_cast<HDSA::Tpetra_Vector<RealT> &>(*(*evec)[i]);
                            sol[0] = evec_i.getVector();
                            postproc_->writeSolution(sol, current_time);
                            current_time = current_time + solver_->deltat;
                        }
                    }
                    else if (HDSA::Tpetra_Vector<RealT> *evec = dynamic_cast<HDSA::Tpetra_Vector<RealT> *>(&(*vec)))
                    {
                        sol[0] = evec->getVector();
                        postproc_->writeSolution(sol, current_time);
                    }
                    postproc_->params->discretized_param_names = discretized_param_names;
                }
            }
            else
            {
                if (HDSA::Tpetra_Vector<RealT> *evec = dynamic_cast<HDSA::Tpetra_Vector<RealT> *>(&(*vec)))
                {
                    postproc_->mesh->setupOptimizationExodusFile(name);
                    postproc_->params->updateParams(evec->getVector());
                    postproc_->writeOptimizationSolution(name);
                }
                else if (HDSA::Std_Vector<RealT> *evec = dynamic_cast<HDSA::Std_Vector<RealT> *>(&(*vec)))
                {
                    vec->Write_to_File(filename + ".txt");
                }
            }
            postproc_->write_solution = false;
        }
        else
        {
            vec->Write_to_File(filename + ".txt");
        }
    }
};
#endif
