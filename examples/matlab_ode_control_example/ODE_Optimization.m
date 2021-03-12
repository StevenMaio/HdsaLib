classdef ODE_Optimization
    
    properties
        ti;
        tf;
        time_mesh;
        h;
        m;
        w;
        n;
        control_dim;    
        U;
    end
    
    methods
        function obj = ODE_Optimization(initial_time,final_time, time_mesh_dofs, control_dim)
            obj.ti = initial_time;
            obj.tf = final_time;
            obj.time_mesh = linspace(obj.ti,obj.tf,time_mesh_dofs);
            obj.h = 1/time_mesh_dofs;
            obj.m = time_mesh_dofs;
            obj.w = (obj.tf-obj.ti)*(1/(obj.m-1))*ones(obj.m,1);
            obj.w(1) = .5*obj.w(1); obj.w(end) = .5*obj.w(end);
            obj.n = length(obj.Set_ODE_IC());
            obj.control_dim = control_dim;
            obj.U = zeros(length(obj.time_mesh),obj.control_dim);
            for k = 1:obj.control_dim
               obj.U(:,k) = obj.Basis_Func(k); 
            end
        end
        
        function [val] = value(this,z,theta,update)
            val = this.cost(z,theta);      
        end
        
        function [grad] = gradient_z(this,z,theta,update)
            [~,grad] = this.cost(z,theta);
        end
        
        function [time_mesh,y,z] = Optimize(this,theta)
            options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'OptimalityTolerance',10^-8);
            z0 = zeros(this.control_dim,1);
            z = fminunc(@(z)this.cost(z,theta),z0,options); 
            y = this.ODE_Solver(z,theta); 
            time_mesh = this.time_mesh;
        end
        
        function [] = Finite_Difference_Test(this,theta)
            z = rand(this.control_dim,1);
            grad = this.gradient_z(z,theta,false);
            val = this.value(z,theta,false);
            FD = zeros(this.control_dim,9);
            h_vec = 10.^-(2:10);
            L = length(h_vec);
            for i = 1:L
                delta = h_vec(i);
                for k = 1:this.control_dim
                   ek = zeros(this.control_dim,1);
                   ek(k) = 1;
                   valk = this.value(z+delta*ek,theta,false);
                   FD(k,i) = (valk-val)/delta;
                end 
            end
            
            for k = 1:this.control_dim
               figure,
               hold on
               plot(1:L,FD(k,:),'o')
               plot(1:L,grad(k)*ones(L,1))
               xticklabels(h_vec)
               xlabel('Step Size')
               ylabel('Partial Derivative')
            end
        end
        
    end
        
    methods (Abstract, Access = protected)

        val = Pointwise_Tracking_Objective(this,t,y,z,theta);
        
        val = Final_Time_Objective(this,y,z,theta);
        
        val = Regularization_Objective(this,z,theta);
        
        grad = Pointwise_Tracking_Objective_State_Grad(this,t,y,z,theta);
        
        grad = Final_Time_Objective_State_Grad(this,y,z,theta);

        grad = Regularization_Objective_Control_Grad(this,z,theta);
        
        y0 = Set_ODE_IC(this);
        
        dy = f(this,t,y,z,theta);
        
        Jac = dfdz(this,t,y,z,theta);

        Jac = dfdy(this,t,y,z,theta);     
        
    end
    
    methods (Access = protected)
        
        function val = Objective_Eval(this,y,z,theta)
            val = this.Pointwise_Tracking_Objective(this.time_mesh,y,z,theta)*this.w ...
                  + this.Final_Time_Objective(y,z,theta) ...
                  + this.Regularization_Objective(z,theta);
        end
        
        function [val,grad] = cost(this,z,theta)
            y = this.ODE_Solver(z,theta);
            val = this.Objective_Eval(y,z,theta);
            if nargout > 1
                lambda = this.Adjoint_ODE_Solve(z,y,theta); 
                grad = zeros(this.control_dim,1);
                for k = 1:(this.m-1)
                   jac = this.dfdz(this.time_mesh(k),y,z,theta);
                   grad = grad + (lambda(k+1,:)*jac)';
                end
                grad = -grad + this.Regularization_Objective_Control_Grad(z,theta);
            end
        end
        
        function y = ODE_Solver(this,z,theta)
           y = zeros(this.n,this.m);
           y(:,1) = this.Set_ODE_IC();
           for k = 2:this.m
              y(:,k) = y(:,k-1) + this.h*this.f(this.time_mesh(k-1),y(:,k-1),z,theta); 
           end
           y = y';
        end
        
        function [val] = Control_Basis_Func(this,t,i)
            val = interp1(this.time_mesh,this.U(:,i),t);
        end
        
        function [u] = Basis_Func(this,i)
              u = chebyshevT(i-1,2*(this.time_mesh-this.ti)/(this.tf-this.ti)-1);
        end
        
        function [lambda] = Adjoint_ODE_Solve(this,z,y,theta)
           lambda = zeros(this.n,this.m);
           lambda(:,this.m) = -this.h*this.w(this.m)*this.Pointwise_Tracking_Objective_State_Grad(this.time_mesh(this.m),y,z,theta) ...
                              -this.h*this.Final_Time_Objective_State_Grad(y,z,theta); 
           for k = (this.m-1):-1:1
              lambda(:,k) = lambda(:,k+1) ...
                  + this.h*this.dfdy(this.time_mesh(k),y,z,theta)'*lambda(:,k+1) ...
                  - this.h*this.w(k)*this.Pointwise_Tracking_Objective_State_Grad(this.time_mesh(k),y,z,theta);
           end
           lambda = lambda';
        end
        
    end

end

