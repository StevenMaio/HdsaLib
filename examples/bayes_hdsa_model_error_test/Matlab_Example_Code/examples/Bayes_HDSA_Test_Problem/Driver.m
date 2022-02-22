%%
clear
close all
clc
addpath('../../src/Optimization/')
addpath('../../src/Model_Form_Error_HDSA/')

obj = Example();
obj.verbose = false;
z0 = rand(51,1)+1;

% Optimization
[u_star,z_star] = obj.Optimize(z0);

% HDSA
hdsa = Example_HDSA(u_star,z_star,obj,obj.M,obj.L,obj.Gamma);
hdsa.Precompute_Data();

num_prior_samps = 400;
num_post_samps = 400;
bayes_hdsa = Example_Bayes_HDSA(hdsa,num_prior_samps,num_post_samps,obj.M,obj.L,obj.Gamma);

[z_prior_samps,delta_prior_samps] = bayes_hdsa.Compute_Prior_Samples();
t = linspace(0,1,hdsa.m)';

figure,
hold on
for k = 1:num_prior_samps
    plot(t,delta_prior_samps(:,k),'color',[.9,.9,.9],'LineWidth',2)
end
title('Prior discrepancy samples','Interpreter','latex')
set(gca,'fontsize', 18)

figure,
hold on
plot(t,z_star,'color','blue','LineWidth',2)
for k = 1:num_prior_samps
    plot(t,z_prior_samps(:,k),'color',[.9,.9,.9],'LineWidth',2)
end
plot(t,z_star,'color','blue','LineWidth',2)
legend('$\overline{z}$','Interpreter','latex')
title('Prior $z$ samples','Interpreter','latex')
set(gca,'fontsize', 18)

%%
obj_hifi = Example_HiFi(obj);
N = 3;
Z = zeros(hdsa.n,N);
Z(:,1) = z_star;
for k = 2:N
   pert = bayes_hdsa.Apply_Sqrt_Gamma_Mat_Inv(randn(length(z_star),1));
   pert = .2*norm(z_star)*pert/norm(pert);
   Z(:,k) = z_star + pert; 
end
Y = obj_hifi.State_Solve(Z) - obj.State_Solve(Z);

alpha = 1e-2;
bayes_hdsa.Compute_Posterior_Data(Z,Y,alpha);

[delta_map_Z,delta_samps_Z,delta_map_Z0] = bayes_hdsa.Compute_Discrepancy_Posterior_Samples();

figure,
hold on
plot(t,delta_map_Z(:,1),'color','red','LineWidth',2)
plot(t,Y(:,1),'color','black','LineWidth',2) 
for k = 1:num_post_samps
   plot(t,delta_samps_Z(:,k,1),'Color',[.9,.9,.9],'LineWidth',2) 
end
plot(t,delta_map_Z(:,1),'color','red','LineWidth',2)
plot(t,Y(:,1),'color','black','LineWidth',2) 
legend({'$\overline{\delta}(z_1)$','$y_1$'},'Interpreter','latex')
title('$\delta(z_1,\theta)$ posterior distribution','Interpreter','latex')
set(gca,'fontsize', 18)

%%
[z_pert,Z_samps] = bayes_hdsa.Compute_Optimal_Solution_Update_Samples();

[u_star_hifi,z_star_hifi] = obj_hifi.Optimize(z0);
norm(z_star - z_star_hifi)/norm(z_star_hifi)
norm((z_star+z_pert) - z_star_hifi)/norm(z_star_hifi)
t = linspace(0,1,hdsa.m)';
figure,
hold on
plot(t,z_star,'color','blue','LineWidth',2)
plot(t,z_star+z_pert,'color','red','LineWidth',2)
plot(t,z_star_hifi,'color','black','LineWidth',2)
for k = 1:num_post_samps
   plot(t,z_star+Z_samps(:,k),'color',[.9,.9,.9],'LineWidth',2) 
end
plot(t,z_star,'color','blue','LineWidth',2)
plot(t,z_star+z_pert,'color','red','LineWidth',2)
plot(t,z_star_hifi,'color','black','LineWidth',2)
xlim([-.05,1])
legend({'$\overline{z}$','$\overline{z}+F_\theta(0)\overline{\theta}$','$z^\star$'},'Interpreter','latex')
title('Optimal solution posterior distribution','Interpreter','latex')
set(gca,'fontsize', 18)

%% Code verification test
verification_test = false;
if verification_test
    
    m = 51;
    n = 51;
    p = m*(n+1);
    Mtheta = zeros(p,p);
    Mtheta(1:m,1:m) = obj.L;
    Mtheta(1:m,(m+1):end) = kron(obj.L,hdsa.M_z');
    Mtheta((m+1):end,1:m) = kron(obj.L,hdsa.M_z);
    Mtheta((m+1):end,(m+1):end) = kron(obj.L,obj.M*(obj.Gamma + z_star*z_star')*obj.M);
    A = zeros(m*N,p);
    for ell = 1:N
        A_ell = [eye(m) , kron(eye(m),Z(:,ell)'*obj.M)];
        A(((ell-1)*m+1):(m*ell),:) = A_ell;
    end
    
    Sigma_inv = Mtheta + (1/alpha)*A'*A;
    b = Y(:);
    
    theta_bar_true = (1/alpha)*linsolve(Sigma_inv,A'*b);
    delta_bar_true_Z = zeros(m,N);
    for ell = 1:N
        delta_bar_true(:,ell) = [eye(m) , kron(eye(m),Z(:,ell)'*obj.M)]*theta_bar_true;
    end
    
    norm(delta_bar_true-delta_map_Z,'fro')/norm(delta_bar_true,'fro')
    
    Sigma_inv_chol = chol(Sigma_inv);
    delta_samps_true = zeros(m,num_post_samps,N);
    theta_samps = linsolve(Sigma_inv_chol,randn(p,num_post_samps)) + theta_bar_true;
    for ell = 1:N
        delta_samps_true(:,:,ell) = [eye(m) , kron(eye(m),Z(:,ell)'*obj.M)]*theta_samps;
    end
    
    figure,
    hold on
    plot(t,delta_bar_true(:,1),'color','red')
    plot(t,Y(:,1),'color','black')
    for k = 1:num_post_samps
        plot(t,delta_samps_true(:,k,1),'Color',[.9,.9,.9])
    end
    plot(t,delta_bar_true(:,1),'color','red')
    plot(t,Y(:,1),'color','black')
    legend({'$\delta(z_1,\theta_{MAP})$','Discrepancy Data'},'Interpreter','latex')
    title('Evaluation of posterior $\delta(z_1,\theta)$ at data point $z_1$','Interpreter','latex')
    
    tmp1 = [eye(m) , kron(eye(m),z_star'*obj.M)]*theta_bar_true;
    tmp2 = hdsa.Apply_u_u_Hessian_FS(tmp1,u_star,z_star);
    tmp3 = hdsa.Apply_Solution_Operator_Jacobian_Transpose(tmp2,u_star,z_star);
    tmp4 = [zeros(n,n), kron(hdsa.g',obj.M)]*theta_bar_true;
    z_pert_true = -hdsa.Apply_Inv_Hessian_RS(tmp3+tmp4,u_star,z_star);
    
    tmp1 = [eye(m) , kron(eye(m),z_star'*obj.M)]*theta_samps;
    tmp2 = hdsa.Apply_u_u_Hessian_FS(tmp1,u_star,z_star);
    tmp3 = hdsa.Apply_Solution_Operator_Jacobian_Transpose(tmp2,u_star,z_star);
    tmp4 = [zeros(n,n), kron(hdsa.g',obj.M)]*theta_samps;
    Z_samps_true = -hdsa.Apply_Inv_Hessian_RS(tmp3+tmp4,u_star,z_star);
    figure,
    hold on
    plot(t,z_star,'color','red')
    plot(t,z_star+z_pert_true,'color','blue')
    plot(t,z_star_hifi,'color','black')
    for k = 1:num_post_samps
        plot(t,z_star+Z_samps_true(:,k),'color',[.9,.9,.9])
    end
    plot(t,z_star,'color','red')
    plot(t,z_star+z_pert_true,'color','blue')
    plot(t,z_star_hifi,'color','black')
    legend({'Nominal solution','MAP Point Prediction','HiFi solution'})
    title('Evaluation of optimal solution distribution')
    
end