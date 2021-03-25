classdef Model_Error_HDSA_Kronecker < handle
    
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
        
        [Mz_v] = Apply_z_Mass(this,v_z);

        [g] = Compute_u_Gradient_FS(this,u,z);
        
        [H_v] = Apply_u_u_Hessian_FS(this,v_u,u,z);

        [J_v] = Apply_Solution_Operator_Jacobian(this,v_z,u,z);
        
        [J_v] = Apply_Solution_Operator_Jacobian_Transpose(this,v_u,u,z);
        
        [L_v] = Apply_L_Operator(this,v_u);
        
        [Linv_v] = Apply_L_Operator_Inv(this,v_u);
        
    end
    
    methods
        function obj = Model_Error_HDSA_Kronecker()

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
        
        function [U,Sigma,Vu,Vz,z1,u2] = Compute_HDSA_GSVD(this,k,p,q)
           kpp = k + p;
           [Q_RS,WQ_RS] = this.RandSubspace(kpp,q);
           
           disp('Projecting onto sampled subspace')
           Bu = zeros(this.m,kpp);
           Bz = zeros(this.n,kpp);
           for j = 1:kpp
               vec = this.Apply_Inv_Hessian_RS(WQ_RS(:,j),this.u_star,this.z_star);
               Bu(:,j) = -this.Apply_u_u_Hessian_FS(this.Apply_Solution_Operator_Jacobian(vec,this.u_star,this.z_star),this.u_star,this.z_star);
               Bz(:,j) = -this.Apply_z_Mass(vec);
           end
           % B(:,k) = kron(Bu(:,k),this.M_z) + kron(this.g,Bz(:,k));
           
           disp('Applying theta mass matrix inverse to projected vectors')
           TinvBu = zeros(this.m,kpp);
           TinvBz = zeros(this.n,kpp);
           for j = 1:kpp
              TinvBu(:,j) = (1/this.coeff)*this.Apply_L_Operator_Inv(Bu(:,j));
              TinvBz(:,j) = this.Apply_E_Inv(Bz(:,j));
           end
           % TinvB(:,k) = kron(TinvBu(:,k),this.Apply_E_Inv(this.M_z)) + kron((1/this.coeff)*this.Apply_L_Operator_Inv(this.g),TinvBz(:,k))

           disp('Orthogonalizing projected vectors')      
           C1 = TinvBu'*Bu*(this.M_z'*this.Apply_E_Inv(this.M_z));
           C2 = (TinvBu'*this.g)*(this.Apply_E_Inv(this.M_z)'*Bz);
           C3 = (1/this.coeff)*(TinvBz'*this.M_z)*(this.Apply_L_Operator_Inv(this.g)'*Bu);
           C4 = ((1/this.coeff)*this.g'*this.Apply_L_Operator_Inv(this.g))*TinvBz'*Bz;
           C = C1 + C2 + C3 + C4;
           RB_test = chol(C);
           QBu = zeros(this.m,kpp);
           QBz = zeros(this.n,kpp);
           for j = 1:size(RB_test,2)
               e = zeros(size(RB_test,2),1);
               e(j) = 1;
               x = linsolve(RB_test,e);
               QBu(:,j) = TinvBu*x;
               QBz(:,j) = TinvBz*x;
           end
           % RB = RB_test and QB(:,k) = kron(QBu(:,k),this.Apply_E_Inv(this.M_z)) + kron((1/this.coeff)*this.Apply_L_Operator_Inv(this.g),QBz(:,k))
           
           [U,Sigma,Vt] = svd(RB_test','econ');
           U = Q_RS*U;
           Vu = QBu*Vt;
           Vz = QBz*Vt;
           
           U = U(:,1:k);
           Vu = Vu(:,1:k);
           Vz = Vz(:,1:k);
           Sigma = Sigma(1:k,1:k);
           z1 = this.Apply_E_Inv(this.M_z);
           u2 = (1/this.coeff)*this.Apply_L_Operator_Inv(this.g);
           % V = kron(Vu,z1) + kron(u2,Vz);
        end
        
        function [Q,WQ] = RandSubspace(this,kpp,q)
            Y = zeros(this.n,kpp);
            disp('Starting first round of subspace iteration matvecs')
            Omega_u = sqrt(this.M_z'*this.M_z)*randn(this.m,kpp);
            Omega_z = sqrt(this.g'*this.g)*randn(this.n,kpp);
            
            for k = 1:kpp
                Bk1 = this.Apply_Solution_Operator_Jacobian_Transpose(this.Apply_u_u_Hessian_FS(Omega_u(:,k),this.u_star,this.z_star),this.u_star,this.z_star);
                Bk2 = this.Apply_z_Mass(Omega_z(:,k));
                Y(:,k) = this.Apply_Inv_Hessian_RS(Bk1+Bk2,this.u_star,this.z_star);
            end
            
            disp('Orthogonalizing range space samples')
            [Q,~,WQ] = this.CholQR(Y,@(z)this.Apply_z_Mass(z));
            
            for j = 1:q
            
                disp('Starting first subspace iteration round of matvecs')                
                Yu = zeros(this.m,kpp);
                Yz = zeros(this.n,kpp);
                WYu = zeros(this.m,kpp);
                WYz = zeros(this.n,kpp);
                for k = 1:kpp
                    v = this.Apply_Inv_Hessian_RS(WQ(:,k),this.u_star,this.z_star);
                    Yu(:,k) = -this.Apply_u_u_Hessian_FS(this.Apply_Solution_Operator_Jacobian(v,this.u_star,this.z_star),this.u_star,this.z_star);
                    Yz(:,k) = -this.Apply_z_Mass(v);
                    WYu(:,k) = (1/this.coeff)*this.Apply_L_Operator_Inv(Yu(:,k));
                    WYz(:,k) = this.Apply_E_Inv(Yz(:,k));
                end
                
                disp('Orthogonalizing row space samples')
                C1 = Yu'*WYu*(this.M_z'*this.Apply_E_Inv(this.M_z));
                C2 = (Yu'*(1/this.coeff)*this.Apply_L_Operator_Inv(this.g))*(this.M_z'*WYz);
                C3 = (Yz'*this.Apply_E_Inv(this.M_z))*(this.g'*WYu);
                C4 = ((1/this.coeff)*this.g'*this.Apply_L_Operator_Inv(this.g))*Yz'*WYz;
                C = C1 + C2 + C3 + C4;
                R = chol(C);
                WQu = zeros(this.m,kpp);
                WQz = zeros(this.n,kpp);
                Qu = zeros(this.m,kpp);
                Qz = zeros(this.n,kpp);
                for k = 1:size(R,2)
                    e = zeros(size(R,2),1);
                    e(k) = 1;
                    x = linsolve(R,e);
                    Qu(:,k) = Yu*x;
                    Qz(:,k) = Yz*x;
                    WQu(:,k) = WYu*x;
                    WQz(:,k) = WYz*x;
                end
                
                disp('Starting second subspace iteration round of matvecs')
                Y = zeros(this.n,kpp);
                for k = 1:kpp
                    vec1 = (this.M_z'*this.Apply_E_Inv(this.M_z))*this.Apply_Solution_Operator_Jacobian_Transpose(this.Apply_u_u_Hessian_FS(WQu(:,k),this.u_star,this.z_star),this.u_star,this.z_star);
                    vec2 = (this.g'*WQu(:,k))*this.Apply_z_Mass(this.Apply_E_Inv(this.M_z));
                    vec3 = (1/this.coeff)*(WQz(:,k)'*this.M_z)*this.Apply_Solution_Operator_Jacobian_Transpose(this.Apply_u_u_Hessian_FS(this.Apply_L_Operator_Inv(this.g),this.u_star,this.z_star),this.u_star,this.z_star);
                    vec4 = (1/this.coeff)*(this.g'*this.Apply_L_Operator_Inv(this.g))*this.Apply_z_Mass(WQz(:,k));
                    vec = vec1 + vec2 + vec3 + vec4;
                    Y(:,k) = this.Apply_Inv_Hessian_RS(vec,this.u_star,this.z_star);
                end
                
                disp('Orthogonalizing range space samples')
                [Q,~,WQ] = this.CholQR(Y,@(z)this.Apply_z_Mass(z));
            
            end
        end
        
        function [Q,R,WQ] = CholQR(this,A,weighting_mat)
            [Z,RA] = qr(A,0);
            X = zeros(size(Z));
            for k = 1:size(Z,2)
               X(:,k) = weighting_mat(Z(:,k)); 
            end
            C = Z'*X;
            RC = chol(C);
            R = RC*RA;
            WQ = zeros(size(Z));
            Q = zeros(size(Z));
            for k = 1:size(RC,2)
               e = zeros(size(RC,2),1);
               e(k) = 1;
               x = linsolve(RC,e);
               Q(:,k) = Z*x;
               WQ(:,k) = X*x;
            end
        end
        
        %% HDSA Operators  
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
                
    end
end

