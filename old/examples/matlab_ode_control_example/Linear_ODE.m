classdef Linear_ODE < ODE_Optimization
    
    properties
        alpha;
    end
    
    methods
        function obj = Linear_ODE(initial_time,final_time,time_mesh_dofs,control_dim,alpha)
                 obj = obj@ODE_Optimization(initial_time,final_time, time_mesh_dofs, control_dim);
                 obj.alpha = alpha;
        end
         
    end
        
    methods(Access = protected)
        
        function val = Pointwise_Tracking_Objective(this,t,y,z,theta)
            val = .5*(interp1(this.time_mesh,y(:,1),t)-3).^2;
        end
        
        function val = Final_Time_Objective(this,y,z,theta)
           val =  .5*(y(end,2)-10)^2;
        end
        
        function val = Regularization_Objective(this,z,theta)
           val =  .5*this.alpha*norm(z)^2;
        end
        
        function grad = Pointwise_Tracking_Objective_State_Grad(this,t,y,z,theta)
           grad = zeros(size(y,2),1); 
           grad(1) = interp1(this.time_mesh,y(:,1),t) - 3;
        end
        
        function grad = Final_Time_Objective_State_Grad(this,y,z,theta)
            grad = zeros(2,1);
            grad(2) = y(end,2)-10;
        end
   
        function grad = Regularization_Objective_Control_Grad(this,z,theta)
            grad = this.alpha*z;
        end
        
        function y0 = Set_ODE_IC(this)
            y0 = zeros(2,1);
            y0(1) = 2;
        end
        
        function dy = f(this,t,y,z,theta)
            forcing = 0;
            for k = 1:this.control_dim
                forcing = forcing + z(k)*this.Control_Basis_Func(t,k); 
            end
            dy = zeros(2,1);
            dy(1) = theta(1)*y(1) + theta(2)*y(2) + forcing;
            dy(2) = theta(3)*y(1) + theta(4)*y(2);
        end
        
        function Jac = dfdz(this,t,y,z,theta)
            Jac = zeros(2,this.control_dim);
            for k = 1:length(z)
                Jac(1,k) = this.Control_Basis_Func(t,k); 
            end
        end

        function Jac = dfdy(this,t,y,z,theta)
        	Jac = zeros(2,2);
            Jac(1,1) = theta(1);
            Jac(1,2) = theta(2);
            Jac(2,1) = theta(3);
            Jac(2,2) = theta(4);
        end        
        
    end
    
   
end
