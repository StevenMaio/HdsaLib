classdef Bayesian_Model_Discrepancy_HDSA < handle
    
    properties
        hdsa_obj;
        num_prior_samps;
        num_post_samps;
        post_data;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        
        % Compute (L+beta*I)^{-1}*v
        [Linv_v] = Apply_L_plus_shift_inv(obj,v,beta);
        
        % Compute L^{-1/2}*v
        [Linv_v] = Apply_Sqrt_L_inv(obj,v);
        
        % Compute (L+beta*I)^{-1/2}*v
        [Linv_v] = Apply_Sqrt_L_plus_shift_inv(obj,v,beta);

        % Compute Gamma^{1/2}*v
        [G_v] = Apply_Sqrt_Gamma_Mat(obj,v);
        
        % Compute Gamma^{-1/2}*v
        [Ginv_v] = Apply_Sqrt_Gamma_Mat_Inv(obj,v);
        
    end
    
    methods
        function obj = Bayesian_Model_Discrepancy_HDSA(hdsa_obj,num_prior_samps,num_post_samps)
            obj.hdsa_obj = hdsa_obj;
            obj.num_prior_samps = num_prior_samps;
            obj.num_post_samps = num_post_samps;
        end
        
        function [z_prior_samps,delta_prior_samps] = Compute_Prior_Samples(obj)
            omega_n = randn(obj.hdsa_obj.n,obj.num_prior_samps);
            Z0_samps = obj.Apply_Sqrt_Gamma_Mat(omega_n);
            U0_samps = obj.Apply_Sqrt_L_inv(randn(obj.hdsa_obj.m,obj.num_prior_samps));
            z_prior_samps = Z0_samps + obj.hdsa_obj.z_star;
            delta_prior_samps = U0_samps*diag(sqrt(1 + diag(omega_n'*omega_n)));
        end
        
        % Compute data to define Bayesian posterior for model discrepancy
        function [] = Compute_Posterior_Data(obj,Z,Y,alpha)
            N = size(Y,2);
            Gamma_inv_Z = obj.hdsa_obj.Apply_Gamma_Mat_Inv(Z);
            G = (1+obj.hdsa_obj.beta) - Z'*obj.hdsa_obj.gamma_inv_z_star - obj.hdsa_obj.gamma_inv_z_star'*Z + Z'*Gamma_inv_Z;
            [g_vecs,Lambda] = eig(G);
            
            u_ell = obj.hdsa_obj.Apply_L_Mat_Inv(Y);
            u_i_ell = cell(N,1);
            for i = 1:N
                u_i_ell{i} = (1/alpha)*obj.Apply_L_plus_shift_inv(u_ell,Lambda(i,i)/alpha);
            end
            
            a = zeros(N,1);
            b = zeros(N,N);
            for ell = 1:N
                a(ell) = 1 - obj.hdsa_obj.gamma_inv_z_star'*(Z(:,ell)-obj.hdsa_obj.z_star);
                for i = 1:N
                    b(i,ell) = (Z*g_vecs(:,i))'*(Gamma_inv_Z(:,ell)-obj.hdsa_obj.gamma_inv_z_star) + sum(g_vecs(:,i))*a(ell);
                end
            end
            
            u_hat = zeros(obj.hdsa_obj.m,obj.num_post_samps,N);
            V = randn(obj.hdsa_obj.m,obj.num_post_samps,N);
            for i=1:N
                u_hat(:,:,i) = (1/sqrt(alpha))*obj.Apply_Sqrt_L_plus_shift_inv(V(:,:,i),Lambda(i,i)/alpha);
            end
            
            obj.post_data = struct;
            obj.post_data.Z = Z;
            obj.post_data.Y = Y;
            obj.post_data.alpha = alpha;
            obj.post_data.Gamma_inv_Z = Gamma_inv_Z;
            obj.post_data.G = G;
            obj.post_data.g_vecs = g_vecs;
            obj.post_data.Lambda = Lambda;
            obj.post_data.u_ell = u_ell;
            obj.post_data.u_i_ell = u_i_ell;
            obj.post_data.a = a;
            obj.post_data.b = b;
            obj.post_data.u_hat = u_hat;
        end
        
        % Compute model discrepancy evaluated at theta MAP point
        function [delta_map] = Compute_Mean_Discrepancy(obj,z)
            N = size(obj.post_data.Y,2);
            dz = z - obj.hdsa_obj.z_star;
            
            coeff = 1 + obj.post_data.Gamma_inv_Z'*dz - obj.hdsa_obj.gamma_inv_z_star'*dz;
            delta_map = (1/obj.post_data.alpha)*obj.post_data.u_ell*coeff;
            
            for i = 1:N
                Gamma_inv_wi = obj.post_data.Gamma_inv_Z*obj.post_data.g_vecs(:,i) - sum(obj.post_data.g_vecs(:,i))*obj.hdsa_obj.gamma_inv_z_star;
                coeff_i = sum(obj.post_data.g_vecs(:,i)) + Gamma_inv_wi'*dz;
                delta_map = delta_map - (1/obj.post_data.alpha)*(obj.post_data.u_i_ell{i}*obj.post_data.b(i,:)')*coeff_i;
            end
            
        end
        
        % Compute model discrepancy samples
        function [delta_map_Z,delta_samps_Z] = Compute_Discrepancy_Posterior_Samples(obj)
            N = size(obj.post_data.Y,2);
            
            delta_map_Z = obj.Compute_Mean_Discrepancy(obj.post_data.Z);
            
            delta_samps_Z = zeros(obj.hdsa_obj.m,obj.num_post_samps,size(delta_map_Z,2));
            Gamma_inv_w_vecs = obj.post_data.Gamma_inv_Z*obj.post_data.g_vecs - obj.hdsa_obj.gamma_inv_z_star*sum(obj.post_data.g_vecs,1);
            for k = 1:N
                coeffs = (obj.post_data.Z(:,k)-obj.hdsa_obj.z_star)'*Gamma_inv_w_vecs + sum(obj.post_data.g_vecs,1);
                coeffs = coeffs./sqrt(diag(obj.post_data.Lambda))';
                coeffs = sqrt(obj.post_data.alpha)*coeffs;
                S = coeffs(1)*obj.post_data.u_hat(:,:,1);
                for j = 2:N
                    S = S + coeffs(j)*obj.post_data.u_hat(:,:,j);
                end
                delta_samps_Z(:,:,k) = delta_map_Z(:,k) + S;
            end
            
        end
        
        % Sample from optimal solution posterior distribution
        function [z_map,Z_samps] = Compute_Optimal_Solution_Update_Samples(obj)
            
            z_map = obj.Mean_Optimal_Solution_Update();
            
            N = size(obj.post_data.Y,2);
            
            coeffs = sum(obj.post_data.g_vecs,1)'./sqrt(diag(obj.post_data.Lambda));
            u_tmp = zeros(obj.hdsa_obj.m,obj.num_post_samps);
            for k = 1:obj.num_post_samps
                u_tmp(:,k) = reshape(obj.post_data.u_hat(:,k,:),obj.hdsa_obj.m,N)*coeffs;
            end
            
            tmp1 = obj.hdsa_obj.Apply_u_u_Hessian_FS(u_tmp,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star);
            z_tmp1 = sqrt(obj.post_data.alpha)*obj.hdsa_obj.Apply_Solution_Operator_Jacobian_Transpose(tmp1,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star);
            
            z_tmp2 = zeros(obj.hdsa_obj.n,obj.num_post_samps);
            for k = 1:obj.num_post_samps
                coeffs = sqrt(obj.post_data.alpha)*(obj.hdsa_obj.g'*reshape(obj.post_data.u_hat(:,k,:),obj.hdsa_obj.m,N))'./sqrt(diag(obj.post_data.Lambda));
                z_tmp2(:,k) = ( obj.post_data.Gamma_inv_Z*obj.post_data.g_vecs-obj.hdsa_obj.gamma_inv_z_star*sum(obj.post_data.g_vecs,1) )*coeffs;
            end
            
            coeff = sqrt(obj.hdsa_obj.Linv_g'*obj.hdsa_obj.g);
            Zhat = obj.post_data.Z(:,2:end)-obj.post_data.Z(:,1)*ones(1,N-1);
            Gamma_inv_Zhat = obj.post_data.Gamma_inv_Z(:,2:end) - obj.post_data.Gamma_inv_Z(:,1)*ones(1,N-1);
            z_hat = obj.Apply_Sqrt_Gamma_Mat_Inv(randn(obj.hdsa_obj.n,obj.num_prior_samps));
            z_tmp3 = coeff*( z_hat - Gamma_inv_Zhat*linsolve(Zhat'*Gamma_inv_Zhat,Zhat'*z_hat) );
            
            Z_samps = -obj.hdsa_obj.Apply_Inv_Hessian_RS(z_tmp1+z_tmp2+z_tmp3,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star) + z_map;
        end
        
        % Compute optimal solution upate at theta MAP point
        function [z_pert] = Mean_Optimal_Solution_Update(obj)
            N = size(obj.post_data.Y,2);
            u = zeros(obj.hdsa_obj.m,1);
            for ell = 1:N
                u = u + obj.post_data.u_ell(:,ell);
                for i = 1:N
                    u = u - obj.post_data.b(i,ell)*sum(obj.post_data.g_vecs(:,i))*obj.post_data.u_i_ell{i}(:,ell);
                end
            end
            tmp1 = obj.hdsa_obj.Apply_u_u_Hessian_FS(u,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star);
            z_tmp = obj.hdsa_obj.Apply_Solution_Operator_Jacobian_Transpose(tmp1,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star);
            
            for ell = 1:N
                z_tmp = z_tmp + (obj.hdsa_obj.g'*obj.post_data.u_ell(:,ell))*(obj.post_data.Gamma_inv_Z(:,ell)-obj.hdsa_obj.gamma_inv_z_star);
                for i = 1:N
                    z_tmp = z_tmp - obj.post_data.b(i,ell)*(obj.hdsa_obj.g'*obj.post_data.u_i_ell{i}(:,ell))*(obj.post_data.Gamma_inv_Z*obj.post_data.g_vecs(:,i)-sum(obj.post_data.g_vecs(:,i))*obj.hdsa_obj.gamma_inv_z_star);
                end
            end
            
            z_tmp = (1/obj.post_data.alpha)*z_tmp;
            z_pert = -obj.hdsa_obj.Apply_Inv_Hessian_RS(z_tmp,obj.hdsa_obj.u_star,obj.hdsa_obj.z_star);
        end
        
    end
    
end