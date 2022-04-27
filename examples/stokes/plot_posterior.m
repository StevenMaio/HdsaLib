clear
close all
clc
rng(21324)

write_to_file = true;

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

N = 3*length(nodes);

load('opt_z_mean.txt')
load('opt_z_samples.txt')
load('delta_mean_at_z1.txt')
load('delta_samples_at_z1.txt')
load('Y.txt')
load('control_read.txt')
load('objective_fun_vals.txt')

M_mat = load('M.txt');
M_mat = diag(M_mat);

state = load('state_read.txt');

data_obj = importdata('NS_Opt_state.txt', ' ', 2);  %% we need to skip the first two lines
ns_opt_state = data_obj.data;

data_obj = importdata('NS_Opt_control.txt', ' ', 2);  %% we need to skip the first two lines
ns_opt_control = data_obj.data;

data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
hifi_state = data_obj.data;
data_obj = importdata('updated_hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
updated_hifi_state = data_obj.data;

delta_at_z1_std = std(delta_samples_at_z1,[],2);
opt_z_std = std(opt_z_samples,[],2);

if write_to_file
    cd figures/
end

m = min([delta_mean_at_z1(1:3:N);Y(1:3:N)]);
M = max([delta_mean_at_z1(1:3:N);Y(1:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('$x$-velocity mean discrepancy $\delta_{v_x}(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'x_velocity_mean_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('$x$-velocity discrepancy data $v_x(z_1)-\tilde{v}_x(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'x_velocity_data_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(1:3:N));
shading interp;
view(0,90)
title('X-Velocity Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([delta_mean_at_z1(2:3:N);Y(2:3:N)]);
M = max([delta_mean_at_z1(2:3:N);Y(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('$y$-velocity mean discrepancy $\delta_{v_y}(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'y_velocity_mean_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('$y$-velocity discrepancy data $v_y(z_1)-\tilde{v}_y(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'y_velocity_data_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(2:3:N));
shading interp;
view(0,90)
title('Y-Velocity Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([delta_mean_at_z1(3:3:N);Y(3:3:N)]);
M = max([delta_mean_at_z1(3:3:N);Y(3:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Pressure mean discrepancy $\delta_p(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'pressure_mean_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Pressure discrepancy data $p(z_1)-\tilde{p}(z_1)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'pressure_data_discrepancy','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(3:3:N));
shading interp;
view(0,90)
title('Pressure Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([opt_z_mean(1:3:N);control_read(1:3:N); ns_opt_control(1:3:N)]);
M = max([opt_z_mean(1:3:N);control_read(1:3:N);ns_opt_control(1:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_mean(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated controller $z_x$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'updated_x_controller','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), control_read(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal controller $z_x$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'nominal_x_controller','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), ns_opt_control(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High-fidelity optimal controller $z_x$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'hifi_x_controller','epsc')
end

m = min([opt_z_mean(2:3:N);control_read(2:3:N); ns_opt_control(2:3:N)]);
M = max([opt_z_mean(2:3:N);control_read(2:3:N);ns_opt_control(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_mean(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated controller $z_y$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'updated_y_controller','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), control_read(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal controller $z_y$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'nominal_y_controller','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), ns_opt_control(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High-fidelity optimal controller $z_y$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'hifi_y_controller','epsc')
end

vec = opt_z_mean;
A = sqrt(M_mat)*(opt_z_samples-vec);
[V,Sigma,~] = svd(A,'econ');
V = inv(sqrt(M_mat))*V;
z_samps = V'*M_mat*(opt_z_samples-vec);
nom_z = V'*M_mat*(control_read-vec);
updated_z = V'*M_mat*(opt_z_mean-vec);
ns_z = V'*M_mat*(ns_opt_control-vec);
%[V,Sigma,~] = svd(opt_z_samples-vec,'econ');
%z_samps = V'*(opt_z_samples-vec);
%nom_z = V'*(control_read-vec);
%updated_z = V'*(opt_z_mean-vec);
%ns_z = V'*(ns_opt_control-vec);
%[~,J] = sort(abs(ns_z),'descend');
J = 1:4;

k = J(1);
figure,
hold on
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
bin_M = max(plt.BinCounts+1);
clf
hold on
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
legend({'Nominal', 'Updated','High-fidelity'},'Location','northeast')
title('First mode coefficient')
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_1_histogram','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(1:3:N,k));
shading interp;
view(0,90)
title('First mode $z_x$ controller','Interpreter','latex')
caxis([min(V(:,k)),max(V(:,k))])
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_1_zx','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(2:3:N,k));
shading interp;
view(0,90)
title('First mode $z_y$ controller','Interpreter','latex')
colorbar(); axis equal
caxis([min(V(:,k)),max(V(:,k))])
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_1_zy','epsc')
end

k = J(2);
figure,
hold on
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
bin_M = max(plt.BinCounts+1);
clf
hold on
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
legend({'Nominal', 'Updated','High-fidelity'},'Location','northeast')
title('Second mode coefficient')
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_2_histogram','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(1:3:N,k));
shading interp;
view(0,90)
title('Second mode $z_x$ controller','Interpreter','latex')
caxis([min(V(:,k)),max(V(:,k))])
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_2_zx','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(2:3:N,k));
shading interp;
view(0,90)
title('Second mode $z_y$ controller','Interpreter','latex')
colorbar(); axis equal
caxis([min(V(:,k)),max(V(:,k))])
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_2_zy','epsc')
end

k = J(3);
figure,
hold on
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
bin_M = max(plt.BinCounts+1);
clf
hold on
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
legend({'Nominal', 'Updated','High-fidelity'},'Location','northeast')
title('Third mode coefficient')
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_3_histogram','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(1:3:N,k));
shading interp;
view(0,90)
title('Third mode $z_x$ controller','Interpreter','latex')
caxis([min(V(:,k)),max(V(:,k))])
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_3_zx','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(2:3:N,k));
shading interp;
view(0,90)
title('Third mode $z_y$ controller','Interpreter','latex')
colorbar(); axis equal
caxis([min(V(:,k)),max(V(:,k))])
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_3_zy','epsc')
end

k = J(4);
figure,
hold on
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
bin_M = max(plt.BinCounts+1);
clf
hold on
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
plt = histogram(z_samps(k,:),'FaceColor',[.8,.8,.8]);
plot([nom_z(k),nom_z(k)],[0,bin_M],'LineWidth',3,'color','blue')
plot([updated_z(k),updated_z(k)],[0,bin_M],'LineWidth',3,'color','black')
plot([ns_z(k),ns_z(k)],[0,bin_M],'--','LineWidth',3,'color','red')
legend({'Nominal', 'Updated','High-fidelity'},'Location','northeast')
title('Fourth mode coefficient')
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_4_histogram','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(1:3:N,k));
shading interp;
view(0,90)
title('Fourth mode $z_x$ controller','Interpreter','latex')
caxis([min(V(:,k)),max(V(:,k))])
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_4_zx','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), V(2:3:N,k));
shading interp;
view(0,90)
title('Fourth mode $z_y$ controller','Interpreter','latex')
colorbar(); axis equal
caxis([min(V(:,k)),max(V(:,k))])
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'posterior_mode_4_zy','epsc')
end

m = min([state(1:3:N); hifi_state(1:3:N);updated_hifi_state(1:3:N);ns_opt_state(1:3:N)]);
M = max([state(1:3:N); hifi_state(1:3:N);updated_hifi_state(1:3:N);ns_opt_state(1:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), state(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal low-fidelity state $\tilde{v}_x(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_lofi_state_v_x','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal high-fidelity state $v_x(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_hifi_state_v_x','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), updated_hifi_state(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated high-fidelity state $v_x(\overline{z}+F_\theta(0)\overline{\theta})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'updated_hifi_state_v_x','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2),ns_opt_state(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Optimal high-fidelity state $v_x(z^\star)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'optimal_hifi_state_v_x','epsc')
end

m = min([state(2:3:N); hifi_state(2:3:N);updated_hifi_state(2:3:N);ns_opt_state(2:3:N)]);
M = max([state(2:3:N); hifi_state(2:3:N);updated_hifi_state(2:3:N);ns_opt_state(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal low-fidelity state $\tilde{v}_y(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_lofi_state_v_y','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal high-fidelity state $v_y(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_hifi_state_v_y','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), updated_hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated high-fidelity state $v_y(\overline{z}+F_\theta(0)\overline{\theta})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'updated_hifi_state_v_y','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2),ns_opt_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Optimal high-fidelity state $v_y(z^\star)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'optimal_hifi_state_v_y','epsc')
end

m = min([state(3:3:N); hifi_state(3:3:N);updated_hifi_state(3:3:N);ns_opt_state(3:3:N)]);
M = max([state(3:3:N); hifi_state(3:3:N);updated_hifi_state(3:3:N);ns_opt_state(3:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), state(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal low-fidelity state $\tilde{p}(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_lofi_state_p','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal high-fidelity state $p(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'nominal_hifi_state_p','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), updated_hifi_state(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated high-fidelity state $p(\overline{z}+F_\theta(0)\overline{\theta})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'updated_hifi_state_p','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2),ns_opt_state(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Optimal high-fidelity state $p(z^\star)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gca,'optimal_hifi_state_p','epsc')
end

if write_to_file 
    cd ../
end
