classdef Model_Form_Error_Kronecker_Vector < handle
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
        function obj = Model_Form_Error_Kronecker_Vector(m,n,K)
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
        
        function [] = Set_Vector(obj,k,a,b_k,u_k,z_k,u,z)
            obj.a = a;
            obj.b_k(k) = b_k;
            obj.u_k(:,k) = u_k;
            obj.z_k(:,k) = z_k;
            obj.u = u;
            obj.z = z;
        end
        
        function [] = Set_Vectors(obj,a,b_k,u_k,z_k,u,z)
            obj.a = a;
            if size(b_k,2)>1
               b_k = b_k'; 
            end
            obj.b_k = b_k;
            obj.u_k = u_k;
            obj.z_k = z_k;
            obj.u = u;
            obj.z = z;
        end
        
        function [val] = Dot_Product(obj,obj_k,v,v_k)
           val = ( obj.a*v.a + obj.z'*v.z ) * ( obj.u_k(:,obj_k)'*v.u_k(:,v_k) );
           val = val + ( obj.a*v.b_k(v_k) + obj.z'*v.z_k(:,v_k) ) * ( obj.u_k(:,obj_k)'*v.u );
           val = val + ( obj.b_k(obj_k)*v.a + obj.z_k(:,obj_k)'*v.z ) * ( obj.u'*v.u_k(:,v_k) );
           val = val + ( obj.b_k(obj_k)*v.b_k(v_k) + obj.z_k(:,obj_k)'*v.z_k(:,v_k) ) * ( obj.u'*v.u );
        end
        
        function [Y] = Construct_Vectors(obj)
            Y = zeros(obj.m*(obj.n+1),obj.K);
            for k = 1:obj.K
               Y(:,k) = [ obj.a*obj.u_k(:,k) ; kron(obj.u_k(:,k),obj.z) ] + [ obj.b_k(k)*obj.u ; kron(obj.u,obj.z_k(:,k)) ];
            end
        end
        
    end
end

