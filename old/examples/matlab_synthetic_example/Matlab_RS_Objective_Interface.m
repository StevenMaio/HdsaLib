classdef Matlab_RS_Objective_Interface
    
    properties
        a;
    end
    
    methods
        function obj = Matlab_RS_Objective_Interface()
            obj.a = (1:10)';
        end
        
        function [val] = value(this,z,theta,update)
            val = sum( (z(1:10)-a.*theta(1:10)).^2 ) + sum(z(11:end).^2);
        end
        
        function [grad] = gradient_z(this,z,theta,update)
            grad = zeros(25,1);
            grad(1:10) = 2*(z(1:10)-this.a.*theta(1:10));
            grad(11:25) = 2*z(11:25);
        end
        
        function [hv] = hessVec_z_z(this,v,z,theta,update)
            hv = 2*v;
        end
        
        function [hv] = hessVec_z_theta(this,v,z,theta,update)
          hv = zeros(25,1);
          hv(1:10) = -2*this.a.*v(1:10);
        end
        
        function [hv] = hessVec_theta_z(this,v,z,theta,update)
          hv = zeros(20,1);
          hv(1:10) = -2*this.a.*v(1:10);
        end
        
    end
end
