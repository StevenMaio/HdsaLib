#ifndef HDSA_MD_POSTERIOR_SAMPLING_HPP
#define HDSA_MD_POSTERIOR_SAMPLING_HPP

namespace HDSA
{

	template <class RealT>
	class MD_Posterior_Sampling
	{

	private:
		HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> data_interface_;
		HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> u_prior_interface_;
		HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT>> z_prior_interface_;

	public:
		HDSA::Ptr<HDSA::MD_Posterior_Data<RealT>> post_data;

		MD_Posterior_Sampling(const HDSA::Ptr<HDSA::MD_Data_Interface<RealT>> &data_interface, const HDSA::Ptr<HDSA::MD_u_Prior_Interface<RealT>> &u_prior_interface, const HDSA::Ptr<HDSA::MD_z_Prior_Interface<RealT>> &z_prior_interface)
		{
			data_interface_ = data_interface;
			u_prior_interface_ = u_prior_interface;
			z_prior_interface_ = z_prior_interface;
			post_data = HDSA::makePtr<HDSA::MD_Posterior_Data<RealT>>();
		}

		~MD_Posterior_Sampling(void)
		{
		}

		void Compute_Posterior_Data(RealT &alpha_d, int &num_samples)
		{
			post_data->Compute_Posterior_Data(*data_interface_, *u_prior_interface_, *z_prior_interface_, alpha_d, num_samples);
		}

		std::vector<HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT>>> Posterior_Discrepancy_Samples(std::vector<HDSA::Ptr<HDSA::Vector<RealT>>> &z) const
		{
			int p = z.size();
			std::vector<HDSA::Ptr<HDSA::MD_Posterior_Vectors<RealT>>> delta;
			delta.resize(p);

			for (int k = 0; k < p; k++)
			{
				delta[k] = HDSA::makePtr<HDSA::MD_Posterior_Vectors<RealT>>(post_data->num_samples, *z[k]);
				HDSA::Ptr<HDSA::Vector<RealT>> delta_mean_k = delta[k]->mean;
				HDSA::Ptr<HDSA::MultiVector<RealT>> delta_samples_k = delta[k]->samples;

				HDSA::Ptr<HDSA::Vector<RealT>> dz_k = z[k]->clone();
				dz_k->set(*z[k]);
				dz_k->axpy(-1.0, *data_interface_->get_z_opt());
				HDSA::Ptr<HDSA::Vector<RealT>> M_z_dz_k = z[k]->clone();
				z_prior_interface_->Apply_M_z(*M_z_dz_k,*dz_k);

				for (int ell = 0; ell < post_data->N; ell++)
				{
					RealT coeff = 1.0 + (*post_data->W_z_inv_M_z_Z)[ell]->dot(*M_z_dz_k) - post_data->W_z_inv_M_z_z_opt->dot(*M_z_dz_k);
					delta_mean_k->axpy(coeff, *(*post_data->u_ell)[ell]);
				}

				for (int i = 0; i < post_data->N; i++)
				{
					// Compute W_z_inv_M_z_yi
					HDSA::Ptr<HDSA::Vector<RealT>> W_z_inv_M_z_yi = dz_k->clone();
					W_z_inv_M_z_yi->axpy(-post_data->sum_g_vecs[i], *post_data->W_z_inv_M_z_z_opt);
					for (int j = 0; j < post_data->N; j++)
					{
						W_z_inv_M_z_yi->axpy((*post_data->g_vecs)(j, i), *(*post_data->W_z_inv_M_z_Z)[j]);
					}

					// Add terms to delta_mean_k
					RealT dz_k_M_z_W_z_inv_M_z_yi = W_z_inv_M_z_yi->dot(*M_z_dz_k);
					for (int ell = 0; ell < post_data->N; ell++)
					{
						RealT coeff = (*post_data->b_i_ell)(i, ell) * (post_data->sum_g_vecs[i] + dz_k_M_z_W_z_inv_M_z_yi);
						delta_mean_k->axpy(-coeff, *(*post_data->u_i_ell[i])[ell]);
					}

					// Add delta_hat terms to delta_samples_k
					RealT coeff = (1.0 / std::sqrt((*post_data->Mu)(i, 0))) * (post_data->sum_g_vecs[i] + dz_k_M_z_W_z_inv_M_z_yi);
					delta_samples_k->axpy(coeff, *post_data->u_i_hat[i]);
				}
				delta_mean_k->scale(1.0 / post_data->alpha_d);
				delta_samples_k->scale(std::sqrt(post_data->alpha_d));

				// Add delta_breve terms to delta_samples_k
				HDSA::Ptr<HDSA::Vector<RealT>> W_z_inv_M_z_dz_k = dz_k->clone();
				z_prior_interface_->Apply_W_z_Inverse(*W_z_inv_M_z_dz_k, *M_z_dz_k);

				// Compute z_tmp = M_z_dz_k - Zc * linsolve(Zc_M_z_W_z_inv_M_z_Zc, M_z_Zc' * W_z_inv_dz_k)
				HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = dz_k->clone();
				z_tmp->set(*M_z_dz_k);

				HDSA::Ptr<HDSA::Dense_Matrix<RealT>> b = post_data->M_z_Zc->MatVec(*W_z_inv_M_z_dz_k);
				HDSA::Ptr<HDSA::Dense_Matrix<RealT>> x = HDSA::makePtr<HDSA::Dense_Matrix<RealT>>(post_data->N - 1, 1);
				HDSA::Linear_Algebra::Symmetric_Direct_Linear_Solve<RealT>(*post_data->Zc_M_z_W_z_inv_M_z_Zc, *x, *b);
				for (int j = 0; j < post_data->N - 1; j++)
				{
					z_tmp->axpy(-(*x)(j, 0), *(*post_data->M_z_Zc)[j]);
				}

				RealT tmp = W_z_inv_M_z_dz_k->dot(*z_tmp);
				if (tmp < -1.e-12)
				{
					HDSA_TEST_FOR_EXCEPTION(true, std::logic_error,
											"Error in Posterior Discrepancy Samples: delta breve coeff < 0" << std::endl);
				}
				RealT breve_coeff = std::sqrt(std::abs(tmp));
				delta_samples_k->axpy(breve_coeff, *post_data->u_breve);

				// Add delta_mean term to delta_samples_k
				delta_samples_k->axpy(1.0, *delta_mean_k);
			}

			return delta;
		}
	};

}

#endif
