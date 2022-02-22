classdef Constrained_Optimization_Model_Form_Error_HDSA < Model_Form_Error_HDSA
    
    properties
        con_opt; % Object of type "Constrained_Optimization"
        hessian_data; % Data needed for hessian vector products
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        [Mz_v] = Apply_z_Mass_Mat(obj,v);
        
        [Mz_v] = Apply_z_Mass_Mat_Inv(obj,v);
        
        [Linv_v] = Apply_L_Mat_Inv(obj,v);
        
        [Ginv_v] = Apply_Gamma_Mat_Inv(obj,v);
        
    end
    
    methods
        function obj = Constrained_Optimization_Model_Form_Error_HDSA(u_star,z_star,con_opt)
            obj = obj@Model_Form_Error_HDSA(u_star,z_star);
            obj.con_opt = con_opt;
        end
        
        %% Implementation of pure virtual functions from base class
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
            [~, g] = obj.con_opt.Objective(u,z);
        end
        
        function [H_v] = Apply_u_u_Hessian_FS(obj,v,u,z)
            H_v = obj.con_opt.Objective_uu_Apply(v,u,z);
        end

        function [J_v] = Apply_Solution_Operator_Jacobian(obj,v,u,z)
            J_v = -obj.con_opt.c_u_Inverse_Apply(obj.con_opt.c_z_Apply(v,u,z),u,z);
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian_Transpose(obj,v,u,z)
            J_v = -obj.con_opt.c_z_Transpose_Apply(obj.con_opt.c_u_Transpose_Inverse_Apply(v,u,z),u,z);
        end
        
    end
        
end

