classdef Model_Form_Error_N_Kronecker_Vector < handle
    % Implementation of a vector of the form
    % y = sum_{k=1}^N [ b_k*u_k ; kron(u_k,z_k) ] 
    
    properties
        b_k; % Vector
        u_k; % Matrix
        z_k; % Matrix
        m; % integer
        n; % integer
        N; % integer
    end
    
    methods
        function obj = Model_Form_Error_N_Kronecker_Vector(m,n,N)
            obj.b_k = zeros(N,1);
            obj.u_k = zeros(m,N);
            obj.z_k = zeros(n,N);
            obj.m = m;
            obj.n = n;
            obj.N = N;
        end
        
        function [] = Set_Vector(obj,b_k,u_k,z_k)
            obj.b_k = b_k;
            obj.u_k = u_k;
            obj.z_k = z_k;
        end
        
        function [theta] = Construct_Vector(obj)
            theta = zeros(obj.m*(obj.n+1),1);
            for k = 1:obj.N
                theta(1:obj.m) = theta(1:obj.m) + obj.b_k(k)*obj.u_k(:,k);
                theta((obj.m+1):end) = theta((obj.m+1):end) + kron(obj.u_k(:,k),obj.z_k(:,k));
            end
        end
           
        function [delta] = Evaluate_Delta(obj,z,hdsa)
           delta = zeros(obj.m,1);
           for k = 1:obj.N
               delta = delta + (obj.b_k(k) + obj.z_k(:,k)'*hdsa.Apply_z_Mass_Mat(z))*obj.u_k(:,k);
           end
        end
        
    end
end

