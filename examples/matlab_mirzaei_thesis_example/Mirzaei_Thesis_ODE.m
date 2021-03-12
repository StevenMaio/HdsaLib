classdef Mirzaei_Thesis_ODE < ODE_Optimization
    
    properties
        alpha;
        final_time;
    end
    
    methods
        function obj = Mirzaei_Thesis_ODE(initial_time,final_time,time_mesh_dofs,control_dofs,num_controls,alpha)
                 obj = obj@ODE_Optimization(initial_time,final_time, time_mesh_dofs, control_dofs,num_controls);
                 obj.alpha = alpha;
                 obj.final_time = final_time;
        end
         
    end
        
    methods(Access = protected)
        
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
        function val = Pointwise_Tracking_Objective(this,t,y,z,theta)
            % This is necessary so that val is an AD type when y is an AD type
            val = .5*50*( ( (y(1) - 8000*exp(-.3*t) )./8000 ).^2 );
        end
        
        % Define vectorized evaluation of tracking type objective
        % t should be a vector of length m
        % y should be a matrix of size mxn
        % val should be a vector of size mx1
        function val = Vectorized_Tracking_Objective(this,t,y,z,theta)
              val = .5*50*( ( (y(:,1) - 8000*exp(-.3*t) )./8000 ).^2 );
        end
        
        % Define the final time objective
        % y should be a vector of length n (solution at final time)
        % val should be a scalar
        function val = Final_Time_Objective(this,y,z,theta)
           val =  -.5*y(4)./5000;
        end
        
        % Define the regularization function
        function val = Regularization_Objective(this,z,theta)
           val =  .5*this.alpha*norm(z)^2;
        end
        
        % Set ODE initial conditions
        % y0 should be a vector of length n
        function y0 = Set_ODE_IC(this)
            y0 = zeros(6,1);
            y0(1) = 8000; y0(2) = 0; y0(3) = 0; y0(4) = 5000; y0(5) = -2*pi/6; y0(6) = 0;
        end
        
        % Define ODE system y'=f via the RHS f
        % t should be a scalar
        % y should be a vecctor of length n
        % dy should be a vector of length n
        % dy should be constructed by defining each component separately
        % and using "vertcat" to concatenate them into a vector, this is
        % necessary to enable AD to compute the Jacobian of f
        function dy = f(this,t,y,z,theta)
            % ODE states indexing (compare with page 16 of OPTIMAL FEEDBACK CONTROL
            % DESIGN FOR A HYPERSONIC REENTRY VEHICLE)
            % y(1) = h (altitude)
            % y(2) = theta (longitude)
            % y(3) = phi (latitude)
            % y(4) = v (velocity)
            % y(5) = gamma (Flight path angle)
            % y(6) = psi (Atmospheric relative heading angle)

            % states
            h = y(1);
            theta_long = y(2);
            phi = y(3);
            v = y(4);
            gamma = y(5);
            psi = y(6);

            % Fine alpha and sigma at time t
            alphat = 0;
            sigmat = 0;
            for k = 1:this.control_dofs
                alphat = alphat + z(k)*this.Control_Basis_Func(t,k);
                sigmat = sigmat + z(k+this.control_dofs)*this.Control_Basis_Func(t,k);
            end

            % hard coded parameters
            re = 6.378*10^6;
            mu = 3.986*10^14;
            rho0 = 1.2;
            H = 7500;
            m = 340;
            Aref = .3;

            % parameters
            Cd = (1.7*alphat^2 + .06)*(1 + theta(1) + theta(2)*alphat + theta(3)*alphat^2);
            Cl = 1.6*alphat*(1 + theta(4) + theta(5)*alphat + theta(6)*alphat^2);
            r = h + re;
            rho = rho0*exp(-h./H)*(1 + theta(7) + theta(8)*(h./H));
            D = (1/2)*rho*v^2*Cd*Aref;
            L = (1/2)*rho*v^2*Cl*Aref;

            dy1 = v*sin(gamma);
            dy2 = (v*cos(gamma)*sin(psi))./(r*cos(phi));
            dy3 = (1./r)*v*cos(gamma)*cos(psi);
            dy4 = -D./m - mu*sin(gamma)./(r^2);
            dy5 = L*cos(sigmat)./(m*v) - mu*cos(gamma)./(v*r^2) + (v./r)*cos(gamma);
            dy6 = L*sin(sigmat)./(m*v*cos(gamma)) + (v./r)*cos(gamma)*sin(psi)*tan(phi);
            
            dy = vertcat(dy1,dy2,dy3,dy4,dy5,dy6);
        end

    end
    
   
end
