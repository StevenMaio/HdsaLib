classdef Model_Error_HDSA_Kronecker < handle
    
    properties
        m;
        n;
        z_cov;
        u_star;
        z_star;
        M_z;
        gamma_inv_z_star;
        Einv_M_z;
        beta;
        g;
        Linv_g;
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
        
        function [] = HDSA_Setup(this,u_star,z_star,z_cov)
            disp('Starting HDSA setup')
            this.m = length(u_star);
            this.n = length(z_star);
            this.u_star = u_star;
            this.z_star = z_star;
            
            this.M_z = this.Apply_z_Mass(z_star);
            this.g = this.Compute_u_Gradient_FS(u_star,z_star);
            this.Linv_g = this.Apply_L_Operator_Inv(this.g);
            this.z_cov = z_cov;    
            this.gamma_inv_z_star = this.z_star./this.z_cov;
            this.beta = this.z_star'*this.gamma_inv_z_star;
            this.Einv_M_z = this.Apply_E_Inv(this.M_z);
            
            disp('Finished HDSA setup')
        end
        
        function [] = Set_Mtheta(this,Mtheta)
           this.Mtheta = Mtheta; 
        end
          
        %% GSVD Functions
        
        function [U,Sigma,V] = Compute_HDSA_GSVD(this,k,p,q)
           kpp = k + p;
           [Q_RS,WQ_RS] = this.RandSubspace(kpp,q);
           
           disp('Projecting onto sampled subspace')
           B = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
           for k = 1:kpp
               vec = this.Apply_Inv_Hessian_RS(WQ_RS(:,k),this.u_star,this.z_star);
               u_k = -this.Apply_u_u_Hessian_FS(this.Apply_Solution_Operator_Jacobian(vec,this.u_star,this.z_star),this.u_star,this.z_star);
               z_k = -this.Apply_z_Mass(vec);
               B.Set_Vector(k,1,0,u_k,z_k,this.g,this.M_z);
           end
           % B(:,k) = kron(Bu(:,k),this.M_z) + kron(this.g,Bz(:,k));
           
           disp('Applying theta mass matrix inverse to projected vectors')
           TinvB = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
           c = 1+this.beta;
           u = c*this.Linv_g;
           z = this.Apply_N(this.M_z) - this.Einv_M_z;
           a = (1-(this.beta/(1+this.beta)));         
           for k = 1:kpp               
               uk = c*this.Apply_L_Operator_Inv(B.u_k(:,k));
               zk = this.Apply_N(B.z_k(:,k));
               bk = -(this.Einv_M_z'*B.z_k(:,k));
               TinvB.Set_Vector(k,a,bk,uk,zk,u,z);
           end

           disp('Orthogonalizing projected vectors')    
           U = (TinvB.u_k'*B.u_k);
           C1 = (TinvB.a*B.a).*U + U.*(TinvB.z'*B.z);
           U = TinvB.u_k'*B.u*ones(1,kpp);
           C2 = (TinvB.a*ones(kpp,1)*B.b_k').*U + U.*(ones(kpp,1)*TinvB.z'*B.z_k);
           U = ones(kpp,1)*TinvB.u'*B.u_k;
           C3 = (TinvB.b_k*B.a*ones(1,kpp)).*U + U.*(TinvB.z_k'*B.z*ones(1,kpp));
           U = (TinvB.u'*B.u);
           C4 = (TinvB.b_k*B.b_k')*U + U.*(TinvB.z_k'*B.z_k);
           
           C = C1 + C2 + C3 + C4;
           R_B = chol(C);
           
           QB = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
           for k = 1:size(R_B,2)
               e = zeros(size(R_B,2),1);
               e(k) = 1;
               x = linsolve(R_B,e);
               QB.Set_Vector(k,TinvB.a,TinvB.b_k'*x,TinvB.u_k*x,TinvB.z_k*x,TinvB.u,TinvB.z);
           end

           [U,Sigma,Vt] = svd(R_B','econ');
           U = Q_RS*U;
           
           V = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
           for k = 1:kpp
              bk = QB.b_k'*Vt(:,k);
              uk = QB.u_k*Vt(:,k);
              zk = QB.z_k*Vt(:,k);
              V.Set_Vector(k,QB.a,bk,uk,zk,QB.u,QB.z);
           end
           
        end
        
        function [Q,WQ] = RandSubspace(this,kpp,q)
            Y = zeros(this.n,kpp);
            disp('Starting first round of subspace iteration matvecs')
            Omega_u = sqrt(1+this.M_z'*this.M_z)*randn(this.m,kpp);
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
                Y = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
                WY = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
                
                c = 1+this.beta;
                Wu = c*this.Linv_g;
                Wz = this.Apply_N(this.M_z) - this.Einv_M_z;
                a = (1-(this.beta/(1+this.beta)));
                    
                for k = 1:kpp
                    v = this.Apply_Inv_Hessian_RS(WQ(:,k),this.u_star,this.z_star);
                    uk = -this.Apply_u_u_Hessian_FS(this.Apply_Solution_Operator_Jacobian(v,this.u_star,this.z_star),this.u_star,this.z_star);
                    zk = -this.Apply_z_Mass(v);
                    Y.Set_Vector(k,1,0,uk,zk,this.g,this.M_z);
                    
                    Wuk = c*this.Apply_L_Operator_Inv(uk);
                    Wzk = this.Apply_N(zk);
                    bk = -(this.Einv_M_z'*zk);
                    WY.Set_Vector(k,a,bk,Wuk,Wzk,Wu,Wz);
                end
                
                disp('Orthogonalizing row space samples')
                U = (Y.u_k'*WY.u_k);
                C1 = (Y.a*WY.a).*U + U.*(Y.z'*WY.z);
                U = Y.u_k'*WY.u*ones(1,kpp);
                C2 = (Y.a*ones(kpp,1)*WY.b_k').*U + U.*(ones(kpp,1)*Y.z'*WY.z_k);
                U = ones(kpp,1)*Y.u'*WY.u_k;
                C3 = (Y.b_k*WY.a*ones(1,kpp)).*U + U.*(Y.z_k'*WY.z*ones(1,kpp));
                U = (Y.u'*WY.u);
                C4 = (Y.b_k*WY.b_k')*U + U.*(Y.z_k'*WY.z_k);
                
                C = C1 + C2 + C3 + C4;        
                R = chol(C);
                
                Q = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
                WQ = Model_Error_Kronecker_Vector(this.m,this.n,kpp);
                for k = 1:size(R,2)
                    e = zeros(size(R,2),1);
                    e(k) = 1;
                    x = linsolve(R,e);
                    Q.Set_Vector(k,Y.a,Y.b_k'*x,Y.u_k*x,Y.z_k*x,Y.u,Y.z);
                    WQ.Set_Vector(k,WY.a,WY.b_k'*x,WY.u_k*x,WY.z_k*x,WY.u,WY.z);
                end
                
                disp('Starting second subspace iteration round of matvecs')
                Y = zeros(this.n,kpp);
                for k = 1:kpp
                    vec1 = (WQ.a+this.M_z'*WQ.z)*this.Apply_Solution_Operator_Jacobian_Transpose(this.Apply_u_u_Hessian_FS(WQ.u_k(:,k),this.u_star,this.z_star),this.u_star,this.z_star);
                    vec2 = (this.g'*WQ.u_k(:,k))*this.Apply_z_Mass(WQ.z);
                    vec3 = (WQ.b_k(k)+this.M_z'*WQ.z_k(:,k))*this.Apply_Solution_Operator_Jacobian_Transpose(this.Apply_u_u_Hessian_FS(WQ.u,this.u_star,this.z_star),this.u_star,this.z_star);
                    vec4 = (this.g'*WQ.u)*this.Apply_z_Mass(WQ.z_k(:,k));
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
        
        function N_v = Apply_N(this,v)
            tmp1 = this.Apply_z_Mass_Inv(v);
            
            tmp1 = tmp1 + (tmp1'*this.gamma_inv_z_star)*this.z_star;
            G_tmp1 = (1/(1+this.beta))*(tmp1'*this.gamma_inv_z_star)*this.gamma_inv_z_star;
            tmp2 = tmp1./this.z_cov - G_tmp1;
            tmp2 = (1/(1+this.beta))*tmp2;
            
            N_v = this.Apply_z_Mass_Inv(tmp2);
        end
        
        function Mzinv_v = Apply_z_Mass_Inv(this,b)
            [Mzinv_v,flag] = pcg(@(z)this.Apply_z_Mass(z),b,10^-8,length(b));
            if flag~=0
                disp('CG Solver Error')
            end
        end
                
    end
end

