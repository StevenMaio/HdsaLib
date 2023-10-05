classdef Process_Perturbations
    
    properties
        theta_sing_vecs_u;
        theta_sing_vecs_uk;
        theta_sing_vecs_z;
        theta_sing_vecs_zk;
        theta_sing_vecs_a;
        theta_sing_vecs_bk;
        M_theta_sing_vecs_z;
        M_theta_sing_vecs_zk;
        z_sing_vecs;
        sing_vals;
        zbar;
    end
    
    methods
        function obj = Process_Perturbations(zbar)
            obj.theta_sing_vecs_u = load('theta_Singular_Vector_u_1.txt');
            obj.theta_sing_vecs_uk = load('theta_Singular_Vector_uk_1.txt');
            obj.theta_sing_vecs_z = load('theta_Singular_Vector_z_1.txt');
            obj.theta_sing_vecs_zk = load('theta_Singular_Vector_zk_1.txt');
            obj.theta_sing_vecs_a = load('theta_Singular_Vector_a_1.txt');
            obj.theta_sing_vecs_bk = load('theta_Singular_Vector_bk_1.txt');
            obj.M_theta_sing_vecs_z = load('theta_Singular_Vector_Mz_1.txt');
            obj.M_theta_sing_vecs_zk = load('theta_Singular_Vector_Mzk_1.txt');
            obj.z_sing_vecs = load('z_Singular_Vector_1.txt');
            sing_vals = load('Singular_Values_1.txt');
            if size(sing_vals,1)==1
               sing_vals = sing_vals'; % Ensure that it is a column vector 
            end
            obj.sing_vals = sing_vals;
            obj.zbar = zbar;
        end
        
        function [z_pert_samps, delta_I_samps,delta_L_u1_samps,delta_L_z1_samps,delta_L_u2_samps,delta_L_z2_samps] = Sample_Nominal_Perturbations(this,S)
            n = length(this.zbar);
            m = length(this.theta_sing_vecs_u);
            z_pert_samps = zeros(n,S);
            delta_I_samps = zeros(m,S);
            delta_L_u1_samps = zeros(m,S);
            delta_L_z1_samps = zeros(n,S);
            delta_L_u2_samps = zeros(m,S);
            delta_L_z2_samps = zeros(n,S);
            K = size(this.theta_sing_vecs_uk,2);
            for k = 1:S
                c = randn(K,1);
                c = c/norm(c);
                [z_pert_samps(:,k),delta_I_samps(:,k),delta_L_u1_samps(:,k),delta_L_z1_samps(:,k),delta_L_u2_samps(:,k),delta_L_z2_samps(:,k)] = this.Evaluate_Perturbation(c);
            end
        end
        
        function [z_pert,delta_I,delta_L_u1,delta_L_z1,delta_L_u2,delta_L_z2] = Evaluate_Perturbation(this,c)
            % Returns solution perturbation z_pert, i.e. new solution is zbar + z_pert
            % Returns perturbation 
            % delta(z,sum c_N*\theta_N) = delta_I + delta_L_u1*delta_L_z1'*(z-zbar) + delta_L_u2*delta_L_z2'*(z-zbar)
            
            z_pert = this.z_sing_vecs*(c.*this.sing_vals);
            
            delta_I = (this.theta_sing_vecs_a+this.zbar'*this.M_theta_sing_vecs_z)*this.theta_sing_vecs_uk*c ...
                    + (c'*this.theta_sing_vecs_bk+this.zbar'*this.M_theta_sing_vecs_zk*c)*this.theta_sing_vecs_u;
               
            delta_L_u1 = this.theta_sing_vecs_uk*c;
            delta_L_z1 = this.M_theta_sing_vecs_z;
            delta_L_u2 = this.theta_sing_vecs_u;
            delta_L_z2 = this.M_theta_sing_vecs_zk*c; 
        end
        
    end
end

