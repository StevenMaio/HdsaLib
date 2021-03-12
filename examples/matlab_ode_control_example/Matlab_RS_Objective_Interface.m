classdef Matlab_RS_Objective_Interface
    
    properties
        solver;
    end
    
    methods
        function obj = Matlab_RS_Objective_Interface()
            
            initial_time = 0;
            final_time = 1;
            time_mesh_dofs = 1000;
            control_dim = 5;
            alpha = 10^-4;

            obj.solver = Linear_ODE(initial_time,final_time,time_mesh_dofs,control_dim,alpha);
        end
        
        function [val] = value(this,z,theta,update)
            val = this.solver.value(z,theta,update);      
        end
        
        function [grad] = gradient_z(this,z,theta,update)
            grad = this.solver.gradient_z(z,theta,update);
        end
        
    end
end

