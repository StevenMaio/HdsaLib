classdef Matlab_RS_Objective_Interface
    
    properties
        solver;
    end
    
    methods
        function obj = Matlab_RS_Objective_Interface()
            
            addpath(genpath('/ascldap/users/joshart/Documents/dasco/Matlab_Trajectory_Analysis_Codes/myAutomaticDifferentiation/'))
            
            initial_time = 0;
            final_time = 1;
            time_mesh_dofs = 100;
            control_dofs = 5;
            num_controls = 2;
            alpha = 0;

            obj.solver = Mirzaei_Thesis_ODE(initial_time,final_time,time_mesh_dofs,control_dofs,num_controls,alpha);
        end
        
        function [val] = value(this,z,theta,update)
            val = this.solver.value(z,theta,update);      
        end
        
        function [grad] = gradient_z(this,z,theta,update)
            grad = this.solver.gradient_z(z,theta,update);
        end
        
    end
end

