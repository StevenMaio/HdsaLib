#ifndef HDSA_WRITE_OUTPUT_MRHYDE_HPP
#define HDSA_WRITE_OUTPUT_MRHYDE_HPP

template <class RealT>
class Write_Output_MrHyDE
{

private:
    bool write_exo_;

public:
    Write_Output_MrHyDE(const HDSA::Ptr<MD_Data_Interface_MrHyDE<RealT>> data_interface)
    {
        std::string opt_solution_exo_file = data_interface->Get_Opt_Solution_Exo_File();
        write_exo_ = false;
        if (opt_solution_exo_file != "error")
        {
            write_exo_ = true;
        }
    }

    virtual ~Write_Output_MrHyDE()
    {
    }

    void Write_Hyperparameters(const HDSA::Ptr<HDSA::MD_u_Hyperparameter_Interface<RealT>> &u_hyperparam_interface, const HDSA::Ptr<HDSA::MD_z_Hyperparameter_Interface<RealT>> &z_hyperparam_interface) const
    {
        std::ofstream outfile("hyperparameters.txt");
        if (outfile.is_open())
        {
            // u hyperparameters
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

            // z hyperparameters
            outfile << "alpha_z: " << z_hyperparam_interface->Get_alpha_z() << std::endl;
            outfile << "beta_z: " << z_hyperparam_interface->Get_beta_z() << std::endl;

            // Close the file
            outfile.close();
        }
        else
        {
            std::cout << "Error: Unable to write hyperparameters to a file" << std::endl;
        }
    }
};
#endif
