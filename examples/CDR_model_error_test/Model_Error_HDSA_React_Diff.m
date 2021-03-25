classdef Model_Error_HDSA_React_Diff < Model_Error_HDSA_Kronecker
    
    properties
        num_data;
        n_mesh;
        beta;
        noise;
        Pe;
        lambda;
        source_nodes;
        B_source_nodes;
        smoothing_coeff;
        
        A_inv;
        A_hat;
        A_hat_inv;
        W_misfit;
        d;
        R;
        R_z;
        L;
        
    end
    
    methods (Access = public)
        
        %% Pure virtual functions
        
        function [u,z] = Solve_Inv_Prob(this)
            C = this.B_source_nodes'*this.A_hat_inv'*this.W_misfit*this.A_hat_inv*this.B_source_nodes + this.beta*(this.R_z'*this.R_z);
            r = zeros(this.n_mesh,1); r(end) = 1;
            b = this.B_source_nodes'*( this.A_hat_inv'*this.W_misfit*this.d - this.A_hat_inv'*this.W_misfit*this.A_hat_inv*r );
            z = linsolve(C,b);
            u = this.A_hat_inv*(this.B_source_nodes*z + r);
        end
        
        function [Hinv_v] = Apply_Inv_Hessian_RS(this,v_z,u,z)
            H = this.B_source_nodes'*this.A_hat_inv'*this.W_misfit*this.A_hat_inv*this.B_source_nodes + this.beta*(this.R_z'*this.R_z);
            Hinv_v = linsolve(H,v_z);
        end
        
        function [Mu_v] = Apply_u_Mass(this,v_u)
            w = ones(this.m,1); w(1) = .5; w(end) = .5; w = w/sum(w);
            Mu_v = w.*v_u;
        end
        
        function [Mz_v] = Apply_z_Mass(this,v_z)
            w = ones(length(v_z),1); w(1) = .5; w(end) = .5; w = w/sum(w);
            Mz_v = w.*v_z;
        end
        
        function [g] = Compute_u_Gradient_FS(this,u,z)
            g = this.W_misfit*(u-this.d);
        end
        
        function [H_v] = Apply_u_u_Hessian_FS(this,v_u,u,z)
            H_v = this.W_misfit*v_u;
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian(this,v_z,u,z)
            J_v = this.A_hat_inv*this.B_source_nodes*v_z;
        end
        
        function [J_v] = Apply_Solution_Operator_Jacobian_Transpose(this,v_u,u,z)
            J_v = this.B_source_nodes'*this.A_hat_inv'*v_u;
        end
        
        function [L_v] = Apply_L_Operator(this,v_u)
            L_v = this.L*v_u;
        end
        
        function [Linv_v] = Apply_L_Operator_Inv(this,v_u)
            Linv_v = linsolve(this.L,v_u);
        end
        
    end
    
    methods
        function obj = Model_Error_HDSA_React_Diff(num_data,n_mesh,beta,noise,Pe,lambda,source_nodes,smoothing_coeff)
            obj = obj@Model_Error_HDSA_Kronecker();
            obj.num_data = num_data;
            obj.n_mesh = n_mesh;
            obj.beta = beta;
            obj.noise = noise;
            obj.Pe = Pe;
            obj.lambda = lambda;
            obj.source_nodes = source_nodes;
            obj.B_source_nodes = eye(n_mesh);
            obj.B_source_nodes = obj.B_source_nodes(:,obj.source_nodes);
            obj.smoothing_coeff = smoothing_coeff;
            
            Set_Up(obj);
        end
        
        function [u_hat,z_hat] = Solve_Perturbed_Inverse_Problem(this,theta)
           
            % Model error function (linear function)
            delta = this.Construct_delta(theta);
            r = zeros(this.n_mesh,1); r(end) = 1;
            C = (this.A_hat_inv*this.B_source_nodes+delta)'*this.W_misfit*(this.A_hat_inv*this.B_source_nodes+delta) + this.beta*(this.R_z'*this.R_z);
            b = (this.A_hat_inv*this.B_source_nodes+delta)'*this.W_misfit*(this.d-this.A_hat_inv*r);
            z_hat = linsolve(C,b);  
            u_hat = this.A_hat_inv*(this.B_source_nodes*z_hat+r) + delta*z_hat;
        end
        
        function [u,z] = Solve_HiFi_Inv_Prob(this)
            C = this.B_source_nodes'*this.A_inv'*this.W_misfit*this.A_inv*this.B_source_nodes + this.beta*(this.R_z'*this.R_z);
            r = zeros(this.n_mesh,1); r(end) = 1;
            b = this.B_source_nodes'*( this.A_inv'*this.W_misfit*this.d - this.A_inv'*this.W_misfit*this.A_inv*r );
            z = linsolve(C,b);
            u = this.A_inv*(this.B_source_nodes*z + r);
        end
        
        function [] = Set_Up(this)
          
            % Need to define A as the discretization of diffusive reaction pde
            % Laplacian
            A = zeros(this.n_mesh,this.n_mesh);
            h = 1/(this.n_mesh-1);
            for k = 2:(this.n_mesh-1)
                A(k,(k-1):(k+1)) = (1/h^2)*[-1, 2,-1] - this.Pe*(1/h)*[-1,1,0] + [0,this.lambda,0];
            end
            % Dirichlet boundary conditions
            A(1,1) = 1;
            A(this.n_mesh,this.n_mesh) = 1;
            % Computing the inverse of "A"
            this.A_inv = zeros(this.n_mesh,this.n_mesh);
            [L,U] = lu(A);
            Ainv = @(b)linsolve(U,linsolve(L,b));
            for k = 1:this.n_mesh
                e = zeros(this.n_mesh,1);
                e(k) = 1;
                this.A_inv(:,k) = Ainv(e);
            end

            % True source
            t = linspace(0,1,this.n_mesh)';
            z_true = exp(-10*(t-.5).^2);
            z_true = z_true(this.source_nodes);
            
            % State solve to generate data
            r = zeros(this.n_mesh,1); r(end)=1;
            this.d = linsolve(A,this.B_source_nodes*z_true+r);
            
            % Derivative norm regularization
            this.R = zeros(this.n_mesh,this.n_mesh);
            for k = 2:(this.n_mesh-1)
                this.R(k,(k-1):(k+1)) = (1/h)*[0,-1,1];
            end
            this.R(1,1:2) = (1/h)*[-1,1];
            
            this.R_z = zeros(length(this.source_nodes),length(this.source_nodes));
            for k = 2:(length(this.source_nodes)-1)
                this.R_z(k,(k-1):(k+1)) = (1/h)*[0,-1,1];
            end
            this.R_z(1,1:2) = (1/h)*[-1,1];
            this.R_z(length(this.source_nodes),(length(this.source_nodes)-1):length(this.source_nodes)) = (1/h)*[-1,1];
            
            % Generating approximate operator
            % Laplacian
            this.A_hat = zeros(this.n_mesh,this.n_mesh);
            h = 1/(this.n_mesh-1);
            for k = 2:(this.n_mesh-1)
                this.A_hat(k,(k-1):(k+1)) = (1/h^2)*[-1,2,-1] - this.Pe*(1/h)*[-1,1,0];
            end
            % Dirichlet boundary conditions
            this.A_hat(1,1) = 1;
            this.A_hat(this.n_mesh,this.n_mesh) = 1;
            
            % Computing the inverse of "A_hat"
            this.A_hat_inv = zeros(this.n_mesh,this.n_mesh);
            [L,U] = lu(this.A_hat);
            A_hatinv = @(b)linsolve(U,linsolve(L,b));
            for k = 1:this.n_mesh
                e = zeros(this.n_mesh,1);
                e(k) = 1;
                this.A_hat_inv(:,k) = A_hatinv(e);
            end
            
            % Generate observation weighting matrix
            this.W_misfit = zeros(this.n_mesh,this.n_mesh);
            obs = randperm(this.n_mesh-2)+1;
            obs = obs(1:this.num_data);
            for k = 1:this.num_data
                this.W_misfit(obs(k),obs(k)) = 1;
            end
            
            this.L = eye(this.n_mesh,this.n_mesh);
            this.L(1,1) = 1/2;
            this.L(end,end) = 1/2;
            this.L = this.L/sum(diag(this.L));
            this.L = this.L + this.smoothing_coeff*this.R'*this.L*this.R;
        end

    end
end

