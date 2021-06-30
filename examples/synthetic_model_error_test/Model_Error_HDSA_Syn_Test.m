classdef Model_Error_HDSA_Syn_Test < Model_Error_HDSA_Kronecker
    
    properties
        n_mesh;
        reg_beta;
        I;
        D;
        L;
        Gamma;
        M;
        d;
        x;
        ztrue;
        utrue;
        
    end
    
    methods (Access = public)
        
        %% Pure virtual functions
        
        function [u,z] = Solve_Inv_Prob(this)
            options = optimoptions('fminunc','SpecifyObjectiveGradient',true,'Display','iter','MaxIterations',10^4,'OptimalityTolerance',10^-7);
            z0 = randn(this.n_mesh,1);
            z = fminunc(@(z)this.J(z),z0,options);
            u = this.Model_Fun(z);
        end
        
        function [Hinv_v] = Apply_Inv_Hessian_RS(this,v_z,u,z)
            c = zeros(this.n_mesh,1);
            [u,u_diff,u_diff_diff] = this.Model_Fun(z);
            c(this.I) = u_diff(this.I).^2 + u_diff_diff(this.I).*(u(this.I) - this.d);
            H = this.reg_beta*this.D + diag(c);
            Hinv_v = linsolve(H,v_z);
        end
        
        function [Mz_v] = Apply_z_Mass(this,v_z)
            Mz_v = this.M*v_z;
        end
        
        function [g] = Compute_u_Gradient_FS(this,u,z)
            g = zeros(length(u),1);
            g(this.I) = (u(this.I)-this.d);
        end
        
        function [H_v] = Apply_u_u_Hessian_FS(this,v_u,u,z)
            H_v = zeros(length(v_u),1);
            H_v(this.I) = v_u(this.I);
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian(this,v_z,u,z)
            [~,diff] = this.Model_Fun(z);
            J_v = diff.*v_z;
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian_Transpose(this,v_u,u,z)
            [~,diff] = this.Model_Fun(z);
            J_v = diff.*v_u;
        end
        
        function [L_v] = Apply_L_Operator(this,v_u)
            L_v = this.L*v_u;
        end
        
        function [Linv_v] = Apply_L_Operator_Inv(this,v_u)
            Linv_v = linsolve(this.L,v_u);
        end
        
        function [Gamma_v] = Apply_Gamma_Operator(this,v_z)
            Gamma_v = this.Gamma*v_z;
        end
        
        function [Gammainv_v] = Apply_Gamma_Operator_Inv(this,v_z)
            Gammainv_v = linsolve(this.Gamma,v_z);
        end
        
    end
    
    methods
        function obj = Model_Error_HDSA_Syn_Test(val_oper)
            obj = obj@Model_Error_HDSA_Kronecker();
            Set_Up(obj);
        end
                
        function [] = Set_Up(this)
          
            n_mesh = 51;
            S = 10;
            reg_beta = 10^-2;
            
            x = linspace(0,1,n_mesh)';
            ztrue = x.^3;
            utrue = this.Model_Fun(ztrue) + this.Model_Error(ztrue);
            
            I = round(linspace(1,n_mesh,S));
            d = utrue(I) + .02*randn(S,1);
            
            h = 1/(n_mesh-1);
            D = (diag(2*ones(1,n_mesh)) + diag(-1*ones(1,n_mesh-1),1) + diag(-1*ones(1,n_mesh-1),-1));
            D(1,1) = 1; D(end,end) = 1;
            D = (1/h)*D;
            
            M = diag(4*ones(1,n_mesh)) + diag(ones(1,n_mesh-1),1) + diag(ones(1,n_mesh-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            
            L = M + 10^-2*D;
            
            Gamma = exp(-1000*(x-x').^2);
            
            this.n_mesh = n_mesh;
            this.reg_beta = reg_beta;
            this.I = I;
            this.D = D;
            this.M = M;
            this.L = L;
            this.Gamma = Gamma;
            this.d = d;
            this.x = x;
            this.ztrue = ztrue;
            this.utrue = utrue;
        end

        function [val,grad] = J(this,z)
            [u,u_diff] = this.Model_Fun(z);
            val =  (1/2)*(this.d-u(this.I))'*(this.d-u(this.I)) + (this.reg_beta/2)*z'*this.D*z;
            grad = this.reg_beta*this.D*z;
            grad(this.I) = grad(this.I) + (this.d-u(this.I)).*(-u_diff(this.I));
        end
        
        function [u,u_diff,u_diff_diff] = Model_Fun(this,z)
           u = z.^3; 
           u_diff = 3*z.^2;
           u_diff_diff = 6*z;
        end
        
        function [err] = Model_Error(this,z)
            err = .01*cos(2*pi*z);
        end
        
    end
end

