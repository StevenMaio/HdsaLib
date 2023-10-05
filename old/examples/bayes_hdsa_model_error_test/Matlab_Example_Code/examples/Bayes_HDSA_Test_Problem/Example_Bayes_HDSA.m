classdef Example_Bayes_HDSA < Bayesian_Model_Discrepancy_HDSA
    
    properties
        M;
        L;
        Gamma;
    end
    
    methods (Access = public)
        
        %% Implementation of pure virtual functions
        
        % Compute (L+beta*I)^{-1}*v
        function [Linv_v] = Apply_L_plus_shift_inv(obj,v,beta)
            Linv_v = linsolve(obj.L + beta*eye(size(v,1)),v);
        end
                
        % Compute L^{-1/2}*v
        function [Linv_v] = Apply_Sqrt_L_inv(obj,v)
            Linv_v = linsolve(sqrtm(obj.L),v);
        end
        
        % Compute (L+beta*I)^{-1/2}*v
        function [Linv_v] = Apply_Sqrt_L_plus_shift_inv(obj,v,beta)
            Linv_v = linsolve(sqrtm(obj.L + beta*eye(size(v,1))),v);
        end
        
        % Compute Gamma^{1/2}*v
        function [G_v] = Apply_Sqrt_Gamma_Mat(obj,v)
           G_v = sqrtm(obj.Gamma)*v; 
        end
        
        % Compute Gamma^{-1/2}*v
        function [Ginv_v] = Apply_Sqrt_Gamma_Mat_Inv(obj,v)
            Ginv_v = linsolve(sqrtm(obj.Gamma),v);
        end
        
    end
    
    methods
        
        function obj = Example_Bayes_HDSA(hdsa_obj,num_prior_samps,num_post_samps,M,L,Gamma)
            obj = obj@Bayesian_Model_Discrepancy_HDSA(hdsa_obj,num_prior_samps,num_post_samps);
            obj.M = M;
            obj.L = L;
            obj.Gamma = Gamma;
        end
        
    end
        
end