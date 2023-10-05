classdef Constrained_ODE_Optimization_Model_Form_Error_HDSA < Model_Form_Error_HDSA
    
    properties
        con_opt; % Object of type "Constrained_Optimization"
        hessian_data; % Data needed for hessian vector products
        M_time; % Time integration weights
        u_star_0;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        [Mz_v] = Apply_z_Mass_Mat(obj,v);
        
        [Mz_v] = Apply_z_Mass_Mat_Inv(obj,v);
        
        [Linv_v] = Apply_Stationary_L_Mat_Inv(obj,v);
        
        [Ginv_v] = Apply_Gamma_Mat_Inv(obj,v);
        
    end
    
    methods
        function obj = Constrained_ODE_Optimization_Model_Form_Error_HDSA(u_star,z_star,con_opt,alpha)
            obj = obj@Model_Form_Error_HDSA(u_star((con_opt.m+1):end),z_star);
            
            obj.u_star_0 = u_star(1:con_opt.m);
            obj.con_opt = con_opt;
            
            M_t = diag((1/3)*ones(1,con_opt.N)) + diag((1/6)*ones(1,con_opt.N-1),1);
            M_t = M_t + M_t';
            M_t(1,1) = 1/3;
            M_t(end,end) = 1/3;
            M_t = (con_opt.T/(con_opt.N-1))*M_t;
            
            M_dt = diag(ones(1,con_opt.N)) + diag(-ones(1,con_opt.N-1),1);
            M_dt = M_dt + M_dt';
            M_dt(1,1) = 1;
            M_dt(end,end) = 1;
            M_dt = ((con_opt.N-1)/(con_opt.T))*M_dt;

            obj.M_time = M_t(2:end,2:end) + alpha*M_dt(2:end,2:end);  
        end
        
        %% Implementation of pure virtual functions from base class
        function [Linv_v] = Apply_L_Mat_Inv(obj,v)
            Linv_v = zeros(obj.con_opt.m*(obj.con_opt.N-1),size(v,2));
            for k = 1:size(v,2)
                vr = reshape(v(:,k),obj.con_opt.m,obj.con_opt.N-1);
                Linv_v_tmp = obj.Apply_Stationary_L_Mat_Inv(vr);
                Linv_v_tmp = linsolve(obj.M_time,Linv_v_tmp')';
                Linv_v(:,k) = Linv_v_tmp(:);
            end
        end
        
        function [Hinv_v] = Apply_Inv_Hessian_RS(obj,v,u,z)
            if isempty(obj.hessian_data)
                [~,~,obj.hessian_data] = obj.con_opt.Jhat(z);
            end
            Hinv_v = zeros(size(v));
            for k = 1:size(v,2)
                [Hinv_v(:,k),flag,relres,iter,resvec] = pcg(@(w)obj.con_opt.Jhat_hessVec(obj.hessian_data,w),v(:,k),10^-6,size(v,1));
                if flag~=0
                    disp('CG Solver Error')
                end
            end
        end
        
        function [g] = Compute_u_Gradient_FS(obj,u,z)
            [~, g] = obj.con_opt.Objective([obj.u_star_0;u],z);
            g = g((obj.con_opt.m+1):end);
        end
        
        function [H_v] = Apply_u_u_Hessian_FS(obj,v,u,z)
            H_v = obj.con_opt.Objective_uu_Apply([zeros(obj.con_opt.m,size(v,2)) ; v],[obj.u_star ; u],z);
            H_v = H_v((obj.con_opt.m+1):end,:);
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian(obj,v,u,z)
            J_v = -obj.con_opt.c_u_Inverse_Apply(obj.con_opt.c_z_Apply(v,[obj.u_star ; u],z),[obj.u_star ; u],z);
            J_v = J_v((obj.con_opt.m+1):end,:);
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian_Transpose(obj,v,u,z)
            J_v = -obj.con_opt.c_z_Transpose_Apply(obj.con_opt.c_u_Transpose_Inverse_Apply([zeros(obj.con_opt.m,size(v,2)) ; v],[obj.u_star ; u],z),[obj.u_star ; u],z);
        end
        
        %% Add zeros to model form error for IC
        function [U,Sigma,V] = Compute_HDSA_GSVD(obj,k,p,q)
            [U,Sigma,V] = Compute_HDSA_GSVD@Model_Form_Error_HDSA(obj,k,p,q);
            V.u_k = [zeros(obj.con_opt.m,V.K) ; V.u_k ];
            V.u = [zeros(obj.con_opt.m,1) ; V.u ];
            V.m = V.m + obj.con_opt.m;
        end
        
        function [z_pert] = Update_Solution_with_Fit_delta(obj,theta_fit)
            coeff_1 = theta_fit.b_k + theta_fit.z_k'*obj.M_z;
            tmp1 = obj.Apply_u_u_Hessian_FS(theta_fit.u_k((obj.con_opt.m+1):end,:),obj.u_star((obj.con_opt.m+1):end),obj.z_star);
            tmp1 = obj.Apply_Solution_Operator_Jacobian_Transpose(tmp1,obj.u_star((obj.con_opt.m+1):end),obj.z_star);
            
            coeff_2 = theta_fit.u_k((obj.con_opt.m+1):end,:)'*obj.g;
            tmp2 = obj.Apply_z_Mass_Mat(theta_fit.z_k);
            
            z_pert = -obj.Apply_Inv_Hessian_RS(tmp1*coeff_1 + tmp2*coeff_2,obj.u_star((obj.con_opt.m+1):end),obj.z_star);
        end
    end
        
end

