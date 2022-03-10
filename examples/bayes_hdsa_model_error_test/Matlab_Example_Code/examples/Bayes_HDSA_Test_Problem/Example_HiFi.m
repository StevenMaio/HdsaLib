classdef Example_HiFi < Constrained_Optimization 
    
    % Implement the optimization problem
    % min_{z} .5*(S(z)-d)^T*(S(z)-d) + .5*beta*z^T*D*z
    % s.t. u = S(z) solves 
    % u = z^3+.3*z^2 on the mesh linspace(0,1,n_mesh)
    
    properties
        example_obj;
        n_mesh;
        reg_beta;
        D;
        L;
        Gamma;
        M;
        d;
    end
    
    
    methods (Access = public)
        function obj = Example_HiFi(example_obj)
            obj = obj@Constrained_Optimization();
            obj.example_obj = example_obj;
        end

    end
      
   
    methods (Access = public)
        
        %% Instantiation of base class pure virtual functions for gradient computation
        function [val,grad_u, grad_z] = Objective(obj,u,z)
            val = .5*(u-obj.example_obj.d)'*(u-obj.example_obj.d) + .5*obj.example_obj.reg_beta*z'*(obj.example_obj.Gamma_bc_inv/sqrt(40))*z;
            grad_u = u-obj.example_obj.d;
            grad_z = obj.example_obj.reg_beta*(obj.example_obj.Gamma_bc_inv/sqrt(40))*z;
        end
        
        function [u] = State_Solve(obj,z) % Input z, evaluate u=S(z)
            u = z.^3 + .3*z.^2;
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z) 
            Mv = v;
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z) 
            Mv = -3*(z.^2).*v - .6*(z.*v);
        end
        
        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z) 
            Mv = v;
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z) 
            Mv = -3*(z.^2).*v - .6*(z.*v);
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda) 
            Mv = 0*u;
        end
        
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda)
            Mv = 0*u;
        end
        
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda) 
            Mv = 0*z;
        end
        
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda) 
            Mv = -6*(lambda.*z).*v - .6*(lambda.*v);
        end
        
        function [Mv] = Objective_uu_Apply(obj,v,u,z) 
            Mv = v;
        end
        
        function [Mv] = Objective_uz_Apply(obj,v,u,z) 
            Mv = 0*u;
        end
        
        function [Mv] = Objective_zu_Apply(obj,v,u,z) 
            Mv = 0*z;
        end
        
        function [Mv] = Objective_zz_Apply(obj,v,u,z)
            Mv = obj.example_obj.reg_beta*(obj.example_obj.Gamma_bc_inv/sqrt(40))*v;
        end

        
    end
end
