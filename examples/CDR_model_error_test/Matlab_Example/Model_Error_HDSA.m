classdef Model_Error_HDSA < handle
    
    properties
        alpha;
        m;
        n;
        z_cov;
        u_star;
        z_star;
        M_z;
        gamma_inv_z_star;
        coeff;
        g;
        Mtheta;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        
        [u,z] = Solve_Inv_Prob(this);
        
        [Hinv_v] = Apply_Inv_Hessian_RS(this,v_z,u,z);
        
        [Mu_v] = Apply_u_Mass(this,v_u);
        
        [Mz_v] = Apply_z_Mass(this,v_z);

        [g] = Compute_u_Gradient_FS(this,u,z);
        
        [H_v] = Apply_u_u_Hessian_FS(this,v_u,u,z);

        [J_v] = Apply_Solution_Operator_Jacobian(this,v_z,u,z);
        
        [J_v] = Apply_Solution_Operator_Jacobian_Transpose(this,v_u,u,z);
        
        [L_v] = Apply_L_Operator(this,v_u);
        
        [Linv_v] = Apply_L_Operator_Inv(this,v_u);
        
    end
    
    methods
        function obj = Model_Error_HDSA()

        end
        
        function [] = HDSA_Setup(this,u_star,z_star,alpha,z_cov)
            disp('Starting HDSA setup')
            this.alpha = alpha;
            this.m = length(u_star);
            this.n = length(z_star);
            this.u_star = u_star;
            this.z_star = z_star;
            
            this.coeff = (1/alpha^2)*(1/(u_star'*this.Apply_L_Operator(u_star)));
            
            this.M_z = this.Apply_z_Mass(z_star);
            this.g = this.Compute_u_Gradient_FS(u_star,z_star);
            this.z_cov = z_cov;    
            this.gamma_inv_z_star = this.z_star./this.z_cov;
            
            disp('Finished HDSA setup')
        end
        
        function [] = Set_Mtheta(this,Mtheta)
           this.Mtheta = Mtheta; 
        end
          
        %% GSVD Functions
        
        function [U,Sigma,V] = Compute_HDSA_GSVD(this,k,p,q)
           kpp = k + p;
           Omega = randn(this.m*this.n,kpp);
           Q_RS = this.RandSubspace(Omega,q);
           
           disp('Projecting onto sampled subspace')
           B = zeros(this.m*this.n,kpp);
           for j = 1:kpp
               B(:,j) = this.Apply_B_Transpose(this.Apply_Inv_Hessian_RS(this.Apply_z_Mass(Q_RS(:,j))));
           end
           
           disp('Applying theta mass matrix inverse to projected vectors')
           TinvB = zeros(this.m*this.n,kpp);
           for j = 1:kpp
              TinvB(:,j) = this.Apply_theta_Mass_Inv(B(:,j)); 
           end

           disp('Orthogonalizing projected vectors')
           [QB,RB] = this.CholQR(TinvB,@(theta)this.Apply_theta_Mass(theta));
           
           [U,Sigma,Vt] = svd(RB','econ');
           U = Q_RS*U;
           V = QB*Vt;
           
           U = U(:,1:k);
           V = V(:,1:k);
           Sigma = Sigma(1:k,1:k);
        end
        
        function [Q] = RandSubspace(this,Omega,q)
            kpp = size(Omega,2);
            Y = zeros(this.n,kpp);
            disp('Starting first round of subspace iteration matvecs')
            for k = 1:kpp
                Y(:,k) = this.Apply_Inv_Hessian_RS(this.Apply_B(Omega(:,k)));
            end
            
            disp('Orthogonalizing range space samples')
            [Q,~] = this.CholQR(Y,@(z)this.Apply_z_Mass(z));
            
            for j = 1:q
            
                disp('Starting first subspace iteration round of matvecs')
                Y = zeros(this.m*this.n,kpp);
                for k = 1:kpp
                    Y(:,k) = this.Apply_B_Transpose(this.Apply_Inv_Hessian_RS(this.Apply_z_Mass(Q(:,k))));
                end
                
                disp('Orthogonalizing row space samples')
                [Q,~] = this.CholQR(Y,@(b)this.Apply_theta_Mass_Inv(b));
                
                disp('Starting second subspace iteration round of matvecs')
                Y = zeros(this.n,kpp);
                for k = 1:kpp
                    Y(:,k) = this.Apply_Inv_Hessian_RS(this.Apply_B(this.Apply_theta_Mass_Inv(Q(:,k))));
                end
                
                disp('Orthogonalizing range space samples')
                [Q,~] = this.CholQR(Y,@(z)this.Apply_z_Mass(z));
            
            end
        end
        
        function [Q,R] = CholQR(this,Z,weighting_mat)
            [QZ,RZ] = qr(Z,0);
            QW = zeros(size(QZ));
            for k = 1:size(QZ,2)
               QW(:,k) = weighting_mat(QZ(:,k)); 
            end
            RW = chol(QZ'*QW);
            R = RW*RZ;
            Q = zeros(size(QZ));
            for k = 1:size(RW,2)
               e = zeros(size(RW,2),1);
               e(k) = 1;
               Q(:,k) = QZ*linsolve(RW,e); 
            end
        end
        
        %% HDSA Operators
        
        function Bv = Apply_B(this,v_theta)
           v1 = this.Apply_A(v_theta);
           v2 = this.Apply_u_u_Hessian_FS(v1,this.u_star,this.z_star);
           v3 = this.Apply_Solution_Operator_Jacobian_Transpose(v2,this.u_star,this.z_star);
           v4 = this.Apply_X(v_theta);
           Bv = v3 + v4;
           Bv = -Bv;
        end
        
        function Btv = Apply_B_Transpose(this,v_z)
            v1 = this.Apply_Solution_Operator_Jacobian(v_z,this.u_star,this.z_star);
            v2 = this.Apply_u_u_Hessian_FS(v1,this.u_star,this.z_star);
            v3 = this.Apply_A_Transpose(v2);
            v4 = this.Apply_X_Transpose(v_z);
            Btv = v3 + v4;
            Btv = -Btv;
        end
        
        function Mtheta_v = Apply_theta_Mass(this,v_theta) 
            
            Mtheta_v = zeros(this.m*this.n,1);
            
            V_tmp = zeros(this.m,this.n);
            for i = 1:this.m
               Ii = ((i-1)*this.n+1):(i*this.n);
               V_tmp(i,:) = this.Apply_E(v_theta(Ii));
            end
            
            U_tmp = zeros(this.m,this.n);
            for k = 1:this.n
                U_tmp(:,k) = this.Apply_L_Operator(V_tmp(:,k));
            end
            
            for i = 1:this.m
                Ii = ((i-1)*this.n+1):(i*this.n);
                Mtheta_v(Ii) = this.coeff*U_tmp(i,:);
            end
        end
        
              
        function Mthetainv_v = Apply_theta_Mass_Inv(this,b)
                        
            Mthetainv_v = zeros(this.m*this.n,1);
            V_tmp = zeros(this.m,this.n);
            for i = 1:this.m
                Ii = ((i-1)*this.n+1):(i*this.n);
                V_tmp(i,:) = this.Apply_E_Inv(b(Ii));
            end
            
            V_tmp = (1/this.coeff)*V_tmp;
            
            U_tmp = zeros(this.m,this.n);
            for k = 1:this.n
                U_tmp(:,k) = this.Apply_L_Operator_Inv(V_tmp(:,k));
            end
            
            for i = 1:this.m
                Ii = ((i-1)*this.n+1):(i*this.n);
                Mthetainv_v(Ii) = U_tmp(i,:);
            end
        end
                
        function E_v = Apply_E(this,v)
            vi = this.Apply_z_Mass(v);
            vi = (this.z_star'*vi)*this.z_star + vi.*this.z_cov;
            E_v =  this.Apply_z_Mass(vi);
        end
        
        function Einv_v = Apply_E_Inv(this,b)
            tmp = this.Apply_z_Mass_Inv(b);
            tmp = tmp./this.z_cov - this.gamma_inv_z_star*(this.gamma_inv_z_star'*tmp)/(1+this.gamma_inv_z_star'*this.z_star);
            Einv_v = this.Apply_z_Mass_Inv(tmp);
        end
        
        function Mzinv_v = Apply_z_Mass_Inv(this,b)
            [Mzinv_v,flag] = pcg(@(z)this.Apply_z_Mass(z),b,10^-8,length(b));
            if flag~=0
                disp('CG Solver Error')
            end
        end
        
        %% Operators needed to compute HDSA operators
        
        % The operator A here corresponds to A*Y in the paper
        function Av = Apply_A(this,v_theta)
            Av = zeros(this.m,1);
            for i = 1:this.m
                Ii = ((i-1)*this.n+1):(i*this.n);
                Av(i) = this.M_z'*v_theta(Ii);
            end
        end
        
        function Atv = Apply_A_Transpose(this,v_z)
            Atv = zeros(this.m*this.n,1);
            for i = 1:this.m
               Ii = ((i-1)*this.n+1):(i*this.n);
               Atv(Ii) = v_z(i)*this.M_z;
            end
        end
        
        function Xv = Apply_X(this,v_theta)
            Xv = zeros(this.n,1);
            for k = 1:this.m
               Ik = ((k-1)*this.n+1):(k*this.n);
               Xv = Xv + this.g(k)*this.Apply_z_Mass(v_theta(Ik));
            end
        end
        
        function Xtv = Apply_X_Transpose(this,v_z)
            Xtv = zeros(this.m*this.n,1);
            Mv = this.Apply_z_Mass(v_z);
            for k = 1:this.m
                Ik = ((k-1)*this.n+1):(k*this.n);
                Xtv(Ik) = this.g(k)*Mv;
            end
        end
        
        %% Testing codes
        function B = Construct_B(this)
           B = zeros(this.n,this.m*this.n);
           for k = 1:this.m*this.n
              disp(['Constructing the ',num2str(k),'th out of ',num2str(this.m*this.n),' column of B'])
              theta = zeros(this.m*this.n,1);
              theta(k) = 1;
              B(:,k) = this.Apply_B(theta);
           end
        end
        
        function Bt = Construct_B_Transpose(this)
            Bt = zeros(this.m*this.n,this.n);
            for k = 1:this.n
                disp(['Constructing the ',num2str(k),'th out of ',num2str(this.n),' column of B Transpose'])
                z = zeros(this.n,1);
                z(k) = 1;
                Bt(:,k) = this.Apply_B_Transpose(z);
            end
        end
        
        function Mu = Construct_u_Mass(this)
            Mu = zeros(this.m,this.m);
            for k = 1:this.m
                u = zeros(this.m,1);
                u(k) = 1;
                Mu(:,k) = this.Apply_u_Mass(u);
            end
        end
        
        function Mz = Construct_z_Mass(this)
            Mz = zeros(this.n,this.n);
            for k = 1:this.n
                z = zeros(this.n,1);
                z(k) = 1;
                Mz(:,k) = this.Apply_z_Mass(z);
            end
        end
        
        function Mtheta = Construct_theta_Mass(this)
            Mtheta = zeros(this.m*this.n,this.m*this.n);
            for k = 1:this.m*this.n
                disp(['Constructing the ',num2str(k),'th out of ',num2str(this.m*this.n),' column of M_\theta'])
                theta = zeros(this.m*this.n,1);
                theta(k) = 1;
                Mtheta(:,k) = this.Apply_theta_Mass(theta);
            end
        end

        function delta = Construct_delta(this,theta)
            delta = zeros(this.m,this.n);
            for k = 1:this.n
                z1 = zeros(this.n,1);
                z1(k) = 1;
                z2 = this.Apply_z_Mass(z1);
                delta(:,k) = reshape(theta,this.n,this.m)'*z2;
            end
        end
        
    end
end

