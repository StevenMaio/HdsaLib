classdef Example_HDSA < Constrained_Optimization_Model_Form_Error_HDSA
    
    properties
        M;
        L;
        Gamma;
    end
    
    methods (Access = public)
        
        %% Implementation of pure virtual functions
        
        function [Hinv_v] = Apply_Inv_Hessian_RS(obj,v,u,z)
            w = 9*z.^4 + 6*z.*(z.^3-obj.con_opt.d); 
            Hr = diag(w) + obj.con_opt.reg_beta*(obj.Gamma/sqrt(40));
            Hinv_v = linsolve(Hr,v);
        end
        
        function [Mz_v] = Apply_z_Mass_Mat(obj,v)
            Mz_v = obj.M*v;
        end
        
        function [Mz_v] = Apply_z_Mass_Mat_Inv(obj,v)
            Mz_v = linsolve(obj.M,v);
        end
        
        function [Linv_v] = Apply_L_Mat_Inv(obj,v)
            Linv_v = linsolve(obj.L,v);
        end
        
        function [Ginv_v] = Apply_Gamma_Mat_Inv(obj,v)
            Ginv_v = linsolve(obj.Gamma,v);
        end
                
    end
    
    methods
        
        function obj = Example_HDSA(u_star,z_star,con_opt,M,L,Gamma)
            obj = obj@Constrained_Optimization_Model_Form_Error_HDSA(u_star,z_star,con_opt);
            obj.M = M;
            obj.L = L;
            obj.Gamma = Gamma;
        end
        
    end
        
end

