classdef ODE_Optimization
    
    properties
        ti;
        tf;
        time_mesh;
        h;
        m;
        w;
        n;
        control_dofs;
        num_controls;
        control_dim;    
        U;
    end
    
    methods
        function obj = ODE_Optimization(initial_time,final_time, time_mesh_dofs, control_dofs,num_controls)
            obj.ti = initial_time;
            obj.tf = final_time;
            obj.time_mesh = linspace(obj.ti,obj.tf,time_mesh_dofs)';
            obj.h = final_time/time_mesh_dofs;
            obj.m = time_mesh_dofs;
            obj.w = (obj.tf-obj.ti)*(1/(obj.m-1))*ones(obj.m,1);
            obj.w(1) = .5*obj.w(1); obj.w(end) = .5*obj.w(end);
            obj.n = length(obj.Set_ODE_IC());
            obj.control_dofs = control_dofs;
            obj.num_controls = num_controls;
            obj.control_dim = num_controls*control_dofs;
            obj.U = zeros(length(obj.time_mesh),obj.control_dofs);
            for k = 1:obj.control_dofs
               obj.U(:,k) = obj.Basis_Func(k); 
            end
        end
        
        function [val] = value(this,z,theta,update)
            val = this.cost(z,theta);      
        end
        
        function [grad] = gradient_z(this,z,theta,update)
            [~,grad] = this.cost(z,theta);
        end
        
        function [hv] = hessVec_z_z(this,v,z,theta,update)
            hv = zeros(length(v),1);
        end
        
        function [time_mesh,y,z] = Optimize(this,theta,iteration_limit)
            options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'OptimalityTolerance',10^-8,'MaxIterations',iteration_limit);
            z0 = rand(this.control_dim,1);
            z = fminunc(@(z)this.cost(z,theta),z0,options); 
            y = this.ODE_Solver(z,theta); 
            time_mesh = this.time_mesh;
        end
        
        function [time_mesh,y] = Evaluate_Uncontrolled_State(this,theta)
            z = zeros(this.control_dim,1);
            y = this.ODE_Solver(z,theta); 
            time_mesh = this.time_mesh;
        end
        
        function [] = Finite_Difference_Gradient_Test(this,theta)
            gradnorm = 10^10;
            while gradnorm > 10^9
                z = rand(this.control_dim,1);
                grad = this.gradient_z(z,theta,false);
                gradnorm = norm(grad);
                if gradnorm > 10^9
                   disp(['Gradient norm =',num2str(gradnorm),', trying another sample']) 
                end
            end
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
        
         function [] = Finite_Difference_Hessian_Test(this,theta)
            gradnorm = 10^10;
            while gradnorm > 10^9
                z = rand(this.control_dim,1);
                v = rand(this.control_dim,1);
                hv = this.hessVec_z_z(v,z,theta,false);
                grad = this.gradient_z(z,theta,false);
                gradnorm = norm(grad);
                if gradnorm > 10^9
                   disp(['Gradient norm =',num2str(gradnorm),', trying another sample']) 
                end
            end
            FD = zeros(this.control_dim,9);
            h_vec = 10.^-(2:10);
            L = length(h_vec);
            scale = norm(z)/norm(v);
            for i = 1:L
                delta = h_vec(i);
                grad_pert = this.gradient_z(z+scale*delta*v,theta,false);
                FD(:,i) = (grad_pert-grad)/(delta*scale);
            end
            
            for k = 1:this.control_dim
               figure,
               hold on
               plot(1:L,FD(k,:),'o')
               plot(1:L,hv(k)*ones(L,1))
               xticklabels(h_vec)
               xlabel('Step Size')
               ylabel('Component of hessVec')
            end
        end
        
    end
        
    methods (Abstract, Access = protected)

        % The ODE solution is a mxn matrix where m is the number of time
        % steps and n is the number of DoFs in the ODE
        % z is the controller and is a vector of length control_dim
        % theta are the parameters and is a vector
        
        % We define a pointwise and vectorized tracking type objective to
        % accomodate the AD codes
        
        % Define pointwise evaluation of tracking type objective
        % t should be a scalar
        % y should be a vector of length n (fixed time snapshot)
        % val should be a scalar
        val = Pointwise_Tracking_Objective(this,t,y,z,theta);
        
        % Define vectorized evaluation of tracking type objective
        % t should be a vector of length m
        % y should be a matrix of size mxn
        % val should be a vector of size mx1
        val = Vectorized_Tracking_Objective(this,t,y,z,theta);
        
        % Define the final time objective
        % y should be a vector of length n (solution at final time)
        % val should be a scalar
        val = Final_Time_Objective(this,y,z,theta);
        
        % Define the regularization function
        val = Regularization_Objective(this,z,theta);
        
        % Set ODE initial conditions
        % y0 should be a vector of length n
        y0 = Set_ODE_IC(this);
        
        % Define ODE system y'=f via the RHS f
        % t should be a scalar
        % y should be a vecctor of length n
        % dy should be a vector of length n
        % dy should be constructed by defining each component separately
        % and using "vertcat" to concatenate them into a vector, this is
        % necessary to enable AD to compute the Jacobian of f
        dy = f(this,t,y,z,theta);  
        
    end
    
    methods (Access = protected)
        
        function grad = Pointwise_Tracking_Objective_State_Grad(this,t,y,z,theta)
            yAD = myAD(y);
            g = this.Pointwise_Tracking_Objective(t,yAD,z,theta);
            grad = getderivs(g)';
        end
        
        function grad = Final_Time_Objective_State_Grad(this,y,z,theta)
            yAD = myAD(y(end,:));
            g = this.Final_Time_Objective(yAD,z,theta);
            grad = getderivs(g)';
        end
        
        function grad = Regularization_Objective_Control_Grad(this,z,theta)
            zAD = myAD(z);
            g = this.Regularization_Objective(zAD,theta);
            grad = getderivs(g)';
        end
        
        function Jac = dfdy(this,t,y,z,theta)
           yAD = myAD(y);
           J = this.f(t,yAD,z,theta);
           Jac = getderivs(J);
        end
        
        function Jac = dfdz(this,t,y,z,theta)
           zAD = myAD(z);
           J = this.f(t,y,zAD,theta);
           Jac = getderivs(J);
        end
        
        function val = Objective_Eval(this,y,z,theta)
            val = this.Vectorized_Tracking_Objective(this.time_mesh,y,z,theta)'*this.w ...
                  + this.Final_Time_Objective(y(end,:),z,theta) ...
                  + this.Regularization_Objective(z,theta);
        end
        
        function [val,grad] = cost(this,z,theta)
            y = this.ODE_Solver(z,theta);
            val = this.Objective_Eval(y,z,theta);
            if nargout > 1
                lambda = this.Adjoint_ODE_Solve(z,y,theta); 
                grad = zeros(this.control_dim,1);
                for k = 1:(this.m-1)
                   yk = interp1(this.time_mesh,y,this.time_mesh(k));
                   jac = this.dfdz(this.time_mesh(k),yk,z,theta);
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
           ym = interp1(this.time_mesh,y,this.time_mesh(this.m));
           lambda(:,this.m) = -this.h*this.w(this.m)*this.Pointwise_Tracking_Objective_State_Grad(this.time_mesh(this.m),ym,z,theta) ...
                              -this.h*this.Final_Time_Objective_State_Grad(y,z,theta);
           for k = (this.m-1):-1:1
               yk = interp1(this.time_mesh,y,this.time_mesh(k));
               lambda(:,k) = lambda(:,k+1) ...
                  + this.h*this.dfdy(this.time_mesh(k),yk,z,theta)'*lambda(:,k+1) ...
                  - this.h*this.w(k)*this.Pointwise_Tracking_Objective_State_Grad(this.time_mesh(k),yk,z,theta);
           end
           lambda = lambda';
        end
        
    end

end

