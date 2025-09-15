#ifndef HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP
#define HDSA_PC_PSEUDO_TIME_CONTINUATION_HPP

#include "HDSA_Vector.hpp"
#include "HDSA_PC_Sensitivity_Operator_Interface.hpp"
#include "HDSA_PC_Quasi_Newton_Preconditioner.hpp"

namespace HDSA
{

	template <class RealT>
	class PC_Pseudo_Time_Continuation
	{

	private:
		HDSA::Ptr<HDSA::Vector<RealT>> z_bar_;
		HDSA::Ptr<HDSA::Vector<RealT>> theta_bar_;
		HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT>> sen_op_interface_;
		HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner<RealT>> qn_prec_;
		RealT grad_tol_;
		bool use_qn_prec_;
		bool print_cg_output_;
		bool print_cg_iter_;
		RealT cg_tol_;
		int max_cg_iter_;

	protected:
		virtual int Apply_Inverse_Hessian(HDSA::Vector<RealT> &z_out, const HDSA::Vector<RealT> &z_in, const HDSA::Vector<RealT> &z, const HDSA::Vector<RealT> &theta, RealT &cg_tol) const
		{

			z_out.zeros();
			RealT z_in_norm = z_in.norm();
			HDSA::Ptr<HDSA::Vector<RealT>> r = z_in.clone();
			r->set(z_in);
			r->scale(1.0 / z_in_norm);
			HDSA::Ptr<HDSA::Vector<RealT>> v = z_in.clone();
			qn_prec_->Apply_Inverse_Hessian_Approximation(*v, *r);
			HDSA::Ptr<HDSA::Vector<RealT>> p = z_in.clone();
			p->set(*v);
			RealT scalar = r->dot(*p);

			RealT r_norm = r->norm();
			int iter = 0;
			HDSA::Ptr<HDSA::Vector<RealT>> w = z_in.clone();

			while ((std::sqrt(scalar) > cg_tol) && (r_norm > cg_tol) && (iter < max_cg_iter_))
			{
				iter += 1;
				sen_op_interface_->Apply_Hessian(*w, *p, z, theta);
				qn_prec_->Add_Block_Quasi_Newton_Step(p, w);

				RealT scalar_tmp = w->dot(*p);
				RealT alpha = scalar / scalar_tmp;
				z_out.axpy(alpha, *p);
				r->axpy(-alpha, *w);
				qn_prec_->Apply_Inverse_Hessian_Approximation(*v, *r);
				RealT scalar_old = scalar;
				scalar = v->dot(*r);
				p->scale(scalar / scalar_old);
				p->plus(*v);
				r_norm = r->norm();

				if (print_cg_iter_)
				{
					std::cout << "Iteration = " << iter << " with relative residual = " << r_norm << std::endl;
				}
			}

			z_out.scale(z_in_norm);

			if (print_cg_output_)
			{
				std::cout << "Total iterations = " << iter << std::endl;
				std::cout << "Relative residual = " << r_norm << std::endl;
			}

			return iter;
		}

	public:
		PC_Pseudo_Time_Continuation(const HDSA::Ptr<HDSA::Vector<RealT>> &z_bar, const HDSA::Ptr<HDSA::Vector<RealT>> &theta_bar,
									const HDSA::Ptr<HDSA::PC_Sensitivity_Operator_Interface<RealT>> &sen_op_interface, const HDSA::Ptr<HDSA::PC_Quasi_Newton_Preconditioner<RealT>> &qn_prec,
									const RealT grad_tol = 1.e-7, const bool use_qn_prec = true, const bool print_cg_output = true, const bool print_cg_iter = false, const RealT cg_tol = 1.e-5, const int max_cg_iter = 100) : z_bar_(z_bar), theta_bar_(theta_bar), sen_op_interface_(sen_op_interface), qn_prec_(qn_prec), grad_tol_(grad_tol), use_qn_prec_(use_qn_prec), print_cg_output_(print_cg_output), print_cg_iter_(print_cg_iter), cg_tol_(cg_tol), max_cg_iter_(max_cg_iter)
		{
		}

		virtual ~PC_Pseudo_Time_Continuation()
		{
		}

		void Pseudo_Time_Continuation_Forward_Euler(HDSA::Vector<RealT> &z_star, HDSA::Vector<RealT> &grad_star, const HDSA::Vector<RealT> &theta_star, const int &N)
		{
			HDSA::Ptr<HDSA::Vector<RealT>> z_current = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> z_new = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> grad_current = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> theta_current = theta_star.clone();

			int num_Hvecs = 0;
			int num_Bvecs = 0;
			int num_grads = 0;

			if (use_qn_prec_)
			{
				num_Hvecs += qn_prec_->Get_Number_HessVecs();
			}

			RealT dt = 1.0 / static_cast<RealT>(N);
			HDSA::Ptr<HDSA::Vector<RealT>> d_theta = theta_star.clone();
			d_theta->set(theta_star);
			d_theta->axpy(-1.0, *theta_bar_);

			z_current->set(*z_bar_);
			theta_current->set(*theta_bar_);
			sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
			num_grads += 1;

			HDSA::Ptr<HDSA::Vector<RealT>> s, y;
			if (use_qn_prec_)
			{
				qn_prec_->Set_N(2 * N + 1);
				s = z_new->clone();
				y = z_new->clone();
			}

			for (int k = 0; k < N; k++)
			{
				std::cout << "-----------------------------------------------------" << std::endl;
				std::cout << "The gradient norm after step: " << k << " is " << grad_current->norm() << std::endl;

				std::cout << "Beginning B matvec at time step " << k + 1 << std::endl;
				sen_op_interface_->Apply_B(*z_tmp, *d_theta, *z_current, *theta_current);
				num_Bvecs += 1;

				std::cout << "Beginning inverse Hessian matvec at Euler step " << k + 1 << std::endl;
				int iters = Apply_Inverse_Hessian(*z_new, *z_tmp, *z_current, *theta_current, cg_tol_);
				num_Hvecs += iters;

				z_new->scale(-dt);
				z_new->plus(*z_current);

				if (use_qn_prec_)
				{
					s->set(*z_new);
					s->axpy(-1.0, *z_current);

					y->set(*grad_current);
					y->axpy(dt, *z_tmp);
					y->scale(-1.0);
				}
				z_current->set(*z_new);
				theta_current->axpy(dt, *d_theta);

				sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
				num_grads += 1;
				if (use_qn_prec_ & (k < N - 1))
				{
					qn_prec_->Add_Block_Quasi_Newton_Data();
					y->plus(*grad_current);
					qn_prec_->Add_Parametric_Quasi_Newton_Data(*s, *y);
				}
				RealT sol_grad_norm = grad_current->norm();

				if (sol_grad_norm > grad_tol_)
				{
					std::cout << "Beginning inverse Hessian matvec at Newton step " << k + 1 << std::endl;
					std::cout << "The current gradient norm = " << sol_grad_norm << std::endl;
					iters = Apply_Inverse_Hessian(*z_new, *grad_current, *z_current, *theta_current, cg_tol_);
					num_Hvecs += iters;

					z_current->axpy(-1.0, *z_new);
					sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
					num_grads += 1;
					if (use_qn_prec_ & (k < N - 1))
					{
						qn_prec_->Add_Block_Quasi_Newton_Data();
					}
					sol_grad_norm = grad_current->norm();
				}

				RealT variable_cg_tol = cg_tol_;
				while (sol_grad_norm > grad_tol_)
				{
					std::cout << "Taking an extra Newton iteration at step " << k + 1 << std::endl;
					std::cout << "The current gradient norm = " << sol_grad_norm << std::endl;
					iters = Apply_Inverse_Hessian(*z_new, *grad_current, *z_current, *theta_current, variable_cg_tol);
					num_Hvecs += iters;
					z_current->axpy(-1.0, *z_new);
					sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
					num_grads += 1;
					sol_grad_norm = grad_current->norm();
					variable_cg_tol = (1.e-1) * variable_cg_tol;
				}
			}

			z_star.set(*z_current);
			grad_star.set(*grad_current);

			RealT sol_grad_norm = grad_current->norm();
			std::cout << " " << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << "Solution gradient norm = " << sol_grad_norm << std::endl;
			std::cout << " " << std::endl;

			std::string name = "Forward_Euler_Cost_Report.txt";
			std::ofstream fout;
			fout.open(name);
			fout << "Number of B-vector products: " << num_Bvecs << std::endl;
			fout << "Number of H-vector products: " << num_Hvecs << std::endl;
			fout << "Number of gradient evaluations: " << num_grads << std::endl;
			fout << "Number of vectors stored for preconditioner: " << qn_prec_->Get_Number_Vecs_Stored() << std::endl;
			fout << "Solution gradient norm = " << sol_grad_norm << std::endl;
			fout.close();
		}

		void Pseudo_Time_Continuation_Modified_Euler(HDSA::Vector<RealT> &z_star, HDSA::Vector<RealT> &grad_star, const HDSA::Vector<RealT> &theta_star, const int &N)
		{
			HDSA::Ptr<HDSA::Vector<RealT>> z_current = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> z_new = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> z_tmp = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> z_store = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> grad_current = z_star.clone();
			HDSA::Ptr<HDSA::Vector<RealT>> theta_current = theta_star.clone();

			int num_Hvecs = 0;
			int num_Bvecs = 0;
			int num_grads = 0;

			if (use_qn_prec_)
			{
				num_Hvecs += qn_prec_->Get_Number_HessVecs();
			}

			RealT dt = 1.0 / static_cast<RealT>(N);
			HDSA::Ptr<HDSA::Vector<RealT>> d_theta = theta_star.clone();
			d_theta->set(theta_star);
			d_theta->axpy(-1.0, *theta_bar_);

			z_current->set(*z_bar_);
			theta_current->set(*theta_bar_);
			sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
			num_grads += 1;

			HDSA::Ptr<HDSA::Vector<RealT>> s, y;
			if (use_qn_prec_)
			{
				qn_prec_->Set_N(3 * N + 1);
				s = z_new->clone();
				y = z_new->clone();
			}

			for (int k = 0; k < N; k++)
			{
				std::cout << "-----------------------------------------------------" << std::endl;
				std::cout << "The gradient norm after step: " << k << " is " << grad_current->norm() << std::endl;

				z_store->set(*z_current);

				std::cout << "Beginning B matvec at Euler step " << k + 1 << std::endl;
				sen_op_interface_->Apply_B(*z_tmp, *d_theta, *z_current, *theta_current);
				num_Bvecs += 1;

				std::cout << "Beginning inverse Hessian matvec at Euler step " << k + 1 << std::endl;
				int iters = Apply_Inverse_Hessian(*z_new, *z_tmp, *z_current, *theta_current, cg_tol_);
				num_Hvecs += iters;

				z_new->scale(-0.5 * dt);
				z_new->plus(*z_current);

				if (use_qn_prec_)
				{
					s->set(*z_new);
					s->axpy(-1.0, *z_current);

					y->set(*grad_current);
					y->axpy(0.5 * dt, *z_tmp);
					y->scale(-1.0);
				}

				z_current->set(*z_new);
				z_new->zeros();
				theta_current->axpy(0.5 * dt, *d_theta);

				sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
				num_grads += 1;
				if (use_qn_prec_)
				{
					qn_prec_->Add_Block_Quasi_Newton_Data();
					y->plus(*grad_current);
					qn_prec_->Add_Parametric_Quasi_Newton_Data(*s, *y);
				}

				std::cout << "Beginning B matvec at Euler step " << k + 1 << " 1/2" << std::endl;
				sen_op_interface_->Apply_B(*z_tmp, *d_theta, *z_current, *theta_current);
				num_Bvecs += 1;

				std::cout << "Beginning inverse Hessian matvec at Euler step " << k + 1 << " 1/2" << std::endl;
				iters = Apply_Inverse_Hessian(*z_new, *z_tmp, *z_current, *theta_current, cg_tol_);
				num_Hvecs += iters;

				z_new->scale(-dt);
				z_new->plus(*z_store);

				if (use_qn_prec_)
				{
					s->set(*z_new);
					s->axpy(-1.0, *z_current);

					y->set(*grad_current);
					y->axpy(0.5 * dt, *z_tmp);
					y->scale(-1.0);
				}

				z_current->set(*z_new);
				theta_current->axpy(0.5 * dt, *d_theta);
				sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
				num_grads += 1;

				if (use_qn_prec_)
				{
					qn_prec_->Add_Block_Quasi_Newton_Data();
					y->plus(*grad_current);
					qn_prec_->Add_Parametric_Quasi_Newton_Data(*s, *y);
				}
				RealT sol_grad_norm = grad_current->norm();

				if (sol_grad_norm > grad_tol_)
				{
					std::cout << "Beginning inverse Hessian matvec at Newton step " << k + 1 << std::endl;
					std::cout << "The current gradient norm = " << sol_grad_norm << std::endl;
					iters = Apply_Inverse_Hessian(*z_new, *grad_current, *z_current, *theta_current, cg_tol_);
					num_Hvecs += iters;

					z_current->axpy(-1.0, *z_new);
					sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
					num_grads += 1;

					if (use_qn_prec_)
					{
						qn_prec_->Add_Block_Quasi_Newton_Data();
					}

					sol_grad_norm = grad_current->norm();
				}

				RealT variable_cg_tol = cg_tol_;
				while (sol_grad_norm > grad_tol_)
				{
					std::cout << "Taking an extra Newton iteration at step " << k + 1 << std::endl;
					std::cout << "The current gradient norm = " << sol_grad_norm << std::endl;
					iters = Apply_Inverse_Hessian(*z_new, *grad_current, *z_current, *theta_current, variable_cg_tol);
					num_Hvecs += iters;
					z_current->axpy(-1.0, *z_new);
					sen_op_interface_->Gradient(*grad_current, *z_current, *theta_current);
					num_grads += 1;
					sol_grad_norm = grad_current->norm();
					variable_cg_tol = (1.e-1) * variable_cg_tol;
				}
			}
			z_star.set(*z_current);
			grad_star.set(*grad_current);

			RealT sol_grad_norm = grad_current->norm();
			std::cout << " " << std::endl;
			std::cout << "-----------------------------------------------------" << std::endl;
			std::cout << "Solution gradient norm = " << sol_grad_norm << std::endl;
			std::cout << " " << std::endl;

			std::string name = "Modified_Euler_Cost_Report.txt";
			std::ofstream fout;
			fout.open(name);
			fout << "Number of B-vector products: " << num_Bvecs << std::endl;
			fout << "Number of H-vector products: " << num_Hvecs << std::endl;
			fout << "Number of gradient evaluations: " << num_grads << std::endl;
			fout << "Number of vectors stored for preconditioner: " << qn_prec_->Get_Number_Vecs_Stored() << std::endl;
			fout << "Solution gradient norm = " << sol_grad_norm << std::endl;
			fout.close();
		}
	};

}

#endif
