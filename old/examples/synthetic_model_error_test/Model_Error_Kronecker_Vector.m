classdef Model_Error_Kronecker_Vector < handle
    % Implementation of K vectors of the form
    % y_k = [ a*u_k ; kron(u_k,z) ] + [ b_k*u ; kron(u,z_k) ]
    
    properties
        a; % Scalar
        b_k; % Vector
        u_k; % Matrix
        z_k; % Matrix
        u; % Vector
        z; % Vector
        m; % integer
        n; % integer
        K; % integer
    end
    
    methods
        function obj = Model_Error_Kronecker_Vector(m,n,K)
            obj.a = 0;
            obj.b_k = zeros(K,1);
            obj.u_k = zeros(m,K);
            obj.z_k = zeros(n,K);
            obj.u = zeros(m,1);
            obj.z = zeros(n,1);
            obj.m = m;
            obj.n = n;
            obj.K = K;
        end
        
        function [] = Set_Vector(this,k,a,b_k,u_k,z_k,u,z)
            this.a = a;
            this.b_k(k) = b_k;
            this.u_k(:,k) = u_k;
            this.z_k(:,k) = z_k;
            this.u = u;
            this.z = z;
        end
        
        function [Y] = Construct_Vectors(this)
            Y = zeros((this.m+1)*this.n,this.K);
            for k = 1:this.K
               Y(:,k) = [ this.a*this.u_k(:,k) ; kron(this.u_k(:,k),this.z) ] + [ this.b_k(k)*this.u ; kron(this.u,this.z_k(:,k)) ];
            end
        end
        
    end
end

