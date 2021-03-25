clear
clc
close all
rng(1234)

n_mesh = 50;
num_data = 20;
source_nodes = 2:49;
beta = 10^-7;
noise = 0.05;
Pe = 5;
lambda = 2.e-1; 
alpha = .001;
smoothing_coeff = .02;

hdsa = Model_Error_HDSA_React_Diff(num_data,n_mesh,beta,noise,Pe,lambda,source_nodes,smoothing_coeff);

[u_star,z_star] = hdsa.Solve_Inv_Prob();
[u_hifi,z_hifi] = hdsa.Solve_HiFi_Inv_Prob();
t = linspace(0,1,n_mesh);
z_true = exp(-10*(t-.5).^2);
z_true = z_true(source_nodes)';

z_cov = 1*ones(length(source_nodes),1);
hdsa.HDSA_Setup(u_star,z_star,alpha,z_cov);

S = 6;
p = 10;
q = 1;
[U,Sigma,V] = hdsa.Compute_HDSA_GSVD(S,p,q);

% axsize = 200;
% figure('Position', [100 680 3*axsize 2*axsize]);
% figure('Position', [750 680 3*axsize 2*axsize]);
% figure('Position', [1400 680 3*axsize 2*axsize]);
% figure(1)
% plot(t,u_star,'color','red','LineWidth',3)
% xlabel('x')
% title('State')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(1)],'epsc')
% 
% figure(2)
% hold on
% plot(t(source_nodes),z_star,'color','red','LineWidth',3)
% plot(t(source_nodes),z_true,'color','black','LineWidth',3)
% title('Source Estimate')
% xlabel('x')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(2)],'epsc')
% 
% figure(3)
% hold on
% plot(t,0*u_star,'color','red','LineWidth',3)
% title('Model Error')
% xlabel('x')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(3)],'epsc')
% 
% pause()
% 
% axsize = 200;
% figure('Position', [100 680 3*axsize 2*axsize]);
% figure('Position', [750 680 3*axsize 2*axsize]);
% figure('Position', [1400 680 3*axsize 2*axsize]);
% figure(4)
% hold on
% plot(t,u_star,'color','red','LineWidth',3)
% plot(t(find(diag(hdsa.W_misfit)==1)),hdsa.d(find(diag(hdsa.W_misfit)==1)),'o','color','black','MarkerSize',8)
% xlabel('x')
% title('State')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(4)],'epsc')
% 
% figure(5)
% hold on
% plot(t(source_nodes),z_star,'color','red','LineWidth',3)
% title('Source Estimate')
% xlabel('x')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(5)],'epsc')
% 
% figure(6)
% hold on
% plot(t,0*u_star,'color','red','LineWidth',3)
% title('Model Error')
% xlabel('x')
% set(gca,'FontSize',18)
% saveas(gcf,['fig',num2str(6)],'epsc')
% 
% pause()
% 
% for j = 1:4
%     h = 1;
%     theta = h*V(:,j);
%     z_pert = z_star+h*Sigma(j,j)*U(:,j);
%     delta = hdsa.Construct_delta(theta);
%     d_z_pert = delta*z_pert;
%     d_z_star = delta*z_star;
%     z_pert_a = zeros(n_mesh,1);
%     z_pert_a(source_nodes) = z_pert;
%     z_pert_a(end) = 1;
%     u_pert = hdsa.A_hat_inv*z_pert_a + delta*z_pert;
%     axsize = 200;
%     figure('Position', [100 680 3*axsize 2*axsize]);
%     figure('Position', [750 680 3*axsize 2*axsize]);
%     figure('Position', [1400 680 3*axsize 2*axsize]);
%     fig_num = 7 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t,u_star,'color','red','LineWidth',3)
%     plot(t,u_pert,'color','blue','LineWidth',3)
%     plot(t(find(diag(hdsa.W_misfit)==1)),hdsa.d(find(diag(hdsa.W_misfit)==1)),'o','color','black','MarkerSize',8)
%     xlabel('x')
%     title('State')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%     
%     fig_num = 8 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t(source_nodes),z_star,'color','red','LineWidth',3)
%     plot(t(source_nodes),z_pert,'color','blue','LineWidth',3)
%     title('Source Estimate')
%     xlabel('x')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%     
%     fig_num = 9 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t,d_z_star,'color','red','LineWidth',3)
%     plot(t,d_z_pert,'color','blue','LineWidth',3)
%     title('Model Error')
%     xlabel('x')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%     
%     pause()
%     
%     figure('Position', [100 680 3*axsize 2*axsize]);
%     figure('Position', [750 680 3*axsize 2*axsize]);
%     figure('Position', [1400 680 3*axsize 2*axsize]);
%     [u_hat,z_hat] = hdsa.Solve_Perturbed_Inverse_Problem(theta);
%     
%     fig_num = 10 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t,u_star,'color','red','LineWidth',3)
%     plot(t,u_pert,'color','blue','LineWidth',3)
%     plot(t,u_hat,'color','black','LineWidth',3)
%     plot(t(find(diag(hdsa.W_misfit)==1)),hdsa.d(find(diag(hdsa.W_misfit)==1)),'o','color','black','MarkerSize',8)
%     xlabel('x')
%     title('State')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%         
%     fig_num = 11 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t(source_nodes),z_star,'color','red','LineWidth',3)
%     plot(t(source_nodes),z_pert,'color','blue','LineWidth',3)
%     plot(t(source_nodes),z_hat,'color','black','LineWidth',3)
%     title('Source Estimate')
%     xlabel('x')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%         
%     fig_num = 12 + 6*(j-1);
%     figure(fig_num)
%     hold on
%     plot(t,d_z_star,'color','red','LineWidth',3)
%     plot(t,d_z_pert,'color','blue','LineWidth',3)
%     title('Model Error')
%     xlabel('x')
%     set(gca,'FontSize',18)
%     saveas(gcf,['fig',num2str(fig_num)],'epsc')
%     
%     pause()
%     
% end
%  
