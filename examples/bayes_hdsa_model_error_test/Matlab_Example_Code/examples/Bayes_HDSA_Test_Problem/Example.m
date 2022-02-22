classdef Example < Constrained_Optimization 
    
    % Implement the optimization problem
    % min_{z} .5*(S(z)-d)^T*(S(z)-d) + .5*beta*z^T*D*z
    % s.t. u = S(z) solves 
    % u = z^3 on the mesh linspace(0,1,n_mesh)
    
    properties
        n_mesh;
        reg_beta;
        D;
        L;
        Gamma_bc;
        Gamma;
        M;
        d;
    end
    
    
    methods (Access = public)
        function obj = Example( )
            obj = obj@Constrained_Optimization();
            n_mesh = 51;
            reg_beta = 10^(-3);
            d = linspace(0,1,n_mesh)';
            
            h = 1/(n_mesh-1);
            D = (diag(2*ones(1,n_mesh)) + diag(-1*ones(1,n_mesh-1),1) + diag(-1*ones(1,n_mesh-1),-1));
            D(1,1) = 1;
            D(1,2) = -1;
            D(end,end) = 1;
            D(end,end-1) = -1;
            D = (1/h)*D;
            
            M = diag(4*ones(1,n_mesh)) + diag(ones(1,n_mesh-1),1) + diag(ones(1,n_mesh-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            
            L = ( 10^1*M + D )'*inv(M)*( 10^1*M + D );
            Linv = linsolve(L,eye(n_mesh));
            L = 5*10^2*diag(sqrt(diag(Linv)))*L*diag(sqrt(diag(Linv)));
            
            Gamma_bc = ( 10^1*M + D )'*inv(M)*( 10^1*M + D );
            Gamma_bc_inv = linsolve(Gamma_bc,eye(n_mesh));
            Gamma = 10^1*diag(sqrt(diag(Gamma_bc_inv)))*Gamma_bc*diag(sqrt(diag(Gamma_bc_inv)));
            
            obj.n_mesh = n_mesh;
            obj.reg_beta = reg_beta;
            obj.D = D;
            obj.Gamma_bc = Gamma_bc;
            obj.Gamma = Gamma;
            obj.M = M;
            obj.L = L;
            obj.d = d;
        end

    end
      
   
    methods (Access = public)
        
        %% Instantiation of base class pure virtual functions for gradient computation
        function [val,grad_u, grad_z] = Objective(obj,u,z)
            val = .5*(u-obj.d)'*(u-obj.d) + .5*obj.reg_beta*z'*(obj.Gamma_bc/sqrt(40))*z;
            grad_u = u-obj.d;
            grad_z = obj.reg_beta*(obj.Gamma_bc/sqrt(40))*z;
        end
        
        function [u] = State_Solve(obj,z) % Input z, evaluate u=S(z)
            u = z.^3;
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z) 
            Mv = v;
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z) 
            Mv = -3*(z.^2).*v;
        end
        
        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z) 
            Mv = v;
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z) 
            Mv = -3*(z.^2).*v;
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
            Mv = -6*(lambda.*z).*v;
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
            Mv = obj.reg_beta*(obj.Gamma_bc/sqrt(40))*v;
        end

        
    end
end

