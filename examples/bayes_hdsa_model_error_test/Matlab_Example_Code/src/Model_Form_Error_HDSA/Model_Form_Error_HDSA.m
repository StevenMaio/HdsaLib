classdef Model_Form_Error_HDSA < handle
    
    properties
        m;
        n;
        u_star;
        z_star;
        M_z;
        gamma_inv_z_star;
        Einv_M_z;
        beta;
        g;
        Linv_g;
        verbose;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        [Hinv_v] = Apply_Inv_Hessian_RS(obj,v,u,z);
        
        [g] = Compute_u_Gradient_FS(obj,u,z);
        
        [H_v] = Apply_u_u_Hessian_FS(obj,v,u,z);

        [J_v] = Apply_Solution_Operator_Jacobian(obj,v,u,z);
        
        [J_v] = Apply_Solution_Operator_Jacobian_Transpose(obj,v,u,z);

        [Mz_v] = Apply_z_Mass_Mat(obj,v);
        
        [Mz_v] = Apply_z_Mass_Mat_Inv(obj,v);
        
        [Linv_v] = Apply_L_Mat_Inv(obj,v);
        
        [Ginv_v] = Apply_Gamma_Mat_Inv(obj,v);
        
    end
    
    methods
        function obj = Model_Form_Error_HDSA(u_star,z_star)
            obj.m = length(u_star);
            obj.n = length(z_star);
            obj.u_star = u_star;
            obj.z_star = z_star;
            obj.verbose = true;
        end
        
        function [] = Precompute_Data(obj)
            obj.M_z = obj.Apply_z_Mass_Mat(obj.z_star);
            obj.g = obj.Compute_u_Gradient_FS(obj.u_star,obj.z_star);
            obj.Linv_g = obj.Apply_L_Mat_Inv(obj.g);
            obj.gamma_inv_z_star = obj.Apply_Gamma_Mat_Inv(obj.z_star);
            obj.beta = obj.z_star'*obj.gamma_inv_z_star;
            obj.Einv_M_z = (1/(1+obj.beta))*obj.Apply_z_Mass_Mat_Inv(obj.gamma_inv_z_star);
        end
                
        %% GSVD Functions
        
        function [U,Sigma,V] = Compute_HDSA_GSVD(obj,k,p,q)
            % Comments will mirrow the algorithm in our model form error paper
            
           % Precompute data needed in the algorithm
           obj.Precompute_Data();
           kpp = k + p;
           
           % Starting call to Randomized GSVD for model error sensitivities
           % algorithm from paper here, Lines 1-10 are in RandSubspace
           [Q_RS,WQ_RS] = obj.RandSubspace(kpp,q);
           
           if obj.verbose
               disp('Projecting onto sampled subspace')
           end
           % Compute W = B^T*H^{-1}*M_z*Q in Line 11
           B = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
           vec = obj.Apply_Inv_Hessian_RS(WQ_RS,obj.u_star,obj.z_star);
           Uk = -obj.Apply_u_u_Hessian_FS(obj.Apply_Solution_Operator_Jacobian(vec,obj.u_star,obj.z_star),obj.u_star,obj.z_star);
           Zk = -obj.Apply_z_Mass_Mat(vec);
           B.Set_Vectors(1,zeros(kpp,1),Uk,Zk,obj.g,obj.M_z);
           % B(:,k) = kron(Bu(:,k),obj.M_z) + kron(obj.g,Bz(:,k));
           
           if obj.verbose
               disp('Applying theta mass matrix inverse to projected vectors')
           end
           % Compute M_\theta^{-1}*W for Line 12
           TinvB = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
           c = (1+obj.beta);
           u = c*obj.Linv_g;
           z = 0*obj.z_star; 
           a = (1-(obj.beta/(1+obj.beta)));         
           Uk = c*obj.Apply_L_Mat_Inv(B.u_k);
           Zk = obj.Apply_N(B.z_k);
           Bk = -(obj.Einv_M_z'*B.z_k);
           TinvB.Set_Vectors(a,Bk,Uk,Zk,u,z);

           if obj.verbose
               disp('Orthogonalizing projected vectors')
           end
           % Compute CholQR in Line 12
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
           
           QB = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);           
           X = linsolve(R_B,eye(size(R_B,2)));
           QB.Set_Vectors(TinvB.a,TinvB.b_k'*X,TinvB.u_k*X,TinvB.z_k*X,TinvB.u,TinvB.z);

           % Compute SVD in Line 13
           [U,Sigma,Vt] = svd(R_B','econ');

           % Compute matrix-products in Line 14
           U = Q_RS*U;
           V = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);           
           Bk = QB.b_k'*Vt;
           Uk = QB.u_k*Vt;
           Zk = QB.z_k*Vt;
           V.Set_Vectors(QB.a,Bk,Uk,Zk,QB.u,QB.z);
           
        end
        
        function [Q,WQ] = RandSubspace(obj,kpp,q)
            if obj.verbose
                disp('Starting first round of subspace iteration matvecs')
            end
            % Random matrix generation in Line 2
            Omega_u = sqrt(1+obj.M_z'*obj.M_z)*randn(obj.m,kpp);
            Omega_z = sqrt(obj.g'*obj.g)*randn(obj.n,kpp);
            
            % Product H^{-1}*B*Omega in Line 3
            Bk1 = obj.Apply_Solution_Operator_Jacobian_Transpose(obj.Apply_u_u_Hessian_FS(Omega_u,obj.u_star,obj.z_star),obj.u_star,obj.z_star);
            Bk2 = obj.Apply_z_Mass_Mat(Omega_z);
            Y = obj.Apply_Inv_Hessian_RS(Bk1+Bk2,obj.u_star,obj.z_star);
            
            if obj.verbose
                disp('Orthogonalizing range space samples')
            end
            % CholQR in Line 4
            [Q,~,WQ] = obj.CholQR(Y,@(z)obj.Apply_z_Mass_Mat(z));
            
            % Loop in Line 5
            for j = 1:q
            
                if obj.verbose
                    disp('Starting first subspace iteration round of matvecs')
                end
                
                % Product Y = B^T*H^{-1}*M*Q in Line 6 and M_theta^{-1}*Y
                % needed for Line 7
                Y = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
                WY = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
                
                c = (1+obj.beta);
                Wu = c*obj.Linv_g;
                Wz = 0*obj.z_star;
                a = (1-(obj.beta/(1+obj.beta)));
                V = obj.Apply_Inv_Hessian_RS(WQ,obj.u_star,obj.z_star);
                Uk = -obj.Apply_u_u_Hessian_FS(obj.Apply_Solution_Operator_Jacobian(V,obj.u_star,obj.z_star),obj.u_star,obj.z_star);
                Zk = -obj.Apply_z_Mass_Mat(V); 
                WUk = c*obj.Apply_L_Mat_Inv(Uk);
                WZk = obj.Apply_N(Zk);
                Bk = -(obj.Einv_M_z'*Zk);
                Y.Set_Vectors(1,zeros(kpp,1),Uk,Zk,obj.g,obj.M_z);
                WY.Set_Vectors(a,Bk,WUk,WZk,Wu,Wz);
                
                if obj.verbose
                    disp('Orthogonalizing row space samples')
                end
                % CholQR in Line 7
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
                
                Q = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
                WQ = Model_Form_Error_Kronecker_Vector(obj.m,obj.n,kpp);
                X = linsolve(R,eye(size(R,2)));
                Q.Set_Vectors(Y.a,Y.b_k'*X,Y.u_k*X,Y.z_k*X,Y.u,Y.z);
                WQ.Set_Vectors(WY.a,WY.b_k'*X,WY.u_k*X,WY.z_k*X,WY.u,WY.z);
                
                if obj.verbose
                    disp('Starting second subspace iteration round of matvecs')
                end
                % Product Y = H^{-1}*B*M_theta^{-1}*Q in Line 8
                tmp2 = obj.Apply_z_Mass_Mat(WQ.z);
                tmp3 = obj.Apply_Solution_Operator_Jacobian_Transpose(obj.Apply_u_u_Hessian_FS(WQ.u,obj.u_star,obj.z_star),obj.u_star,obj.z_star);
                vec1 = (WQ.a+obj.M_z'*WQ.z)*obj.Apply_Solution_Operator_Jacobian_Transpose(obj.Apply_u_u_Hessian_FS(WQ.u_k,obj.u_star,obj.z_star),obj.u_star,obj.z_star);
                vec2 = tmp2*(obj.g'*WQ.u_k);
                vec3 = tmp3*(WQ.b_k'+obj.M_z'*WQ.z_k);
                vec4 = (obj.g'*WQ.u)*obj.Apply_z_Mass_Mat(WQ.z_k);
                vec = vec1 + vec2 + vec3 + vec4;
                Y = obj.Apply_Inv_Hessian_RS(vec,obj.u_star,obj.z_star);
                
                if obj.verbose
                    disp('Orthogonalizing range space samples')
                end
                % CholQR in Line 9
                [Q,~,WQ] = obj.CholQR(Y,@(z)obj.Apply_z_Mass_Mat(z));
            
            end
        end
        
        % Cholesky QR algorithm
        function [Q,R,WQ] = CholQR(obj,A,weighting_mat)
            [Z,RA] = qr(A,0);
            X = weighting_mat(Z);
            C = Z'*X;
            RC = chol(C);
            R = RC*RA;
            x = linsolve(RC,eye(size(RC,2)));
            Q = Z*x;
            WQ = X*x;
        end
        
        % Weighting matrix helper function                 
        function N_v = Apply_N(obj,v)
            tmp1 = obj.Apply_z_Mass_Mat_Inv(v);
            tmp2 = (1/(1+obj.beta))*obj.Apply_Gamma_Mat_Inv(tmp1);
            N_v = obj.Apply_z_Mass_Mat_Inv(tmp2);
        end
        
        function X_v = Apply_X(obj,v)
            v = v + obj.z_star*(v'*obj.gamma_inv_z_star)';
            G_tmp = (1/(1+obj.beta))*obj.gamma_inv_z_star*(v'*obj.gamma_inv_z_star)';
            tmp = obj.Apply_Gamma_Mat_Inv(v) - G_tmp;
            X_v = (1/(1+obj.beta))*tmp;
        end
                
        %% Postprocessing
        function [z_pert,delta_I,delta_L_u1,delta_L_z1,delta_L_u2,delta_L_z2] = Evalute_delta(obj,c,U,Sigma,V,zbar)
            % Returns solution perturbation z_pert, i.e. new solution is zbar + z_pert
            % Returns perturbation
            % delta(z,sum c_N*\theta_N) = delta_I + delta_L_u1*delta_L_z1'*(z-zbar) + delta_L_u2*delta_L_z2'*(z-zbar)
            
            z_pert = U*(Sigma*c);
            
            Mz = obj.Apply_z_Mass_Mat(V.z);
            Mz_k = obj.Apply_z_Mass_Mat(V.z_k);
            
            delta_I = (V.a+zbar'*Mz)*V.u_k*c + (c'*V.b_k+zbar'*Mz_k*c)*V.u;
            
            delta_L_u1 = V.u_k*c;
            delta_L_z1 = Mz;
            delta_L_u2 = V.u;
            delta_L_z2 = Mz_k*c;
        end
        
        %% Fit to model form error data
        function [theta_fit] = Fit_delta(obj,Z,Y)
            N = size(Y,2);
            MZ = obj.Apply_z_Mass_Mat(Z);
            XZ = obj.Apply_X(Z);
            F = (1+obj.beta) - Z'*obj.gamma_inv_z_star - obj.gamma_inv_z_star'*Z + Z'*obj.Apply_Gamma_Mat_Inv(Z);
            
            [f,a] = eig(F);
            
            b_k = zeros(N,1);
            u_k = zeros(obj.m,N);
            z_k = zeros(obj.n,N);
            for k = 1:N
                b_k(k) = (sum(f(:,k))-obj.Einv_M_z'*MZ*f(:,k));
                u_k(:,k) = ( (1+obj.beta)/a(k,k) )*Y*f(:,k);
                z_k(:,k) = -sum(f(:,k))*obj.Einv_M_z + obj.Apply_z_Mass_Mat_Inv(XZ*f(:,k));
            end
            
            theta_fit = Model_Form_Error_N_Kronecker_Vector(obj.m,obj.n,N);
            theta_fit.Set_Vector(b_k,u_k,z_k);
        end
        
        function [z_pert] = Update_Solution_with_Fit_delta(obj,theta_fit)
            coeff_1 = theta_fit.b_k + theta_fit.z_k'*obj.M_z;
            tmp1 = obj.Apply_u_u_Hessian_FS(theta_fit.u_k,obj.u_star,obj.z_star);
            tmp1 = obj.Apply_Solution_Operator_Jacobian_Transpose(tmp1,obj.u_star,obj.z_star);
            
            coeff_2 = theta_fit.u_k'*obj.g;
            tmp2 = obj.Apply_z_Mass_Mat(theta_fit.z_k);
            
            z_pert = -obj.Apply_Inv_Hessian_RS(tmp1*coeff_1 + tmp2*coeff_2,obj.u_star,obj.z_star);             
        end
        
    end
end

