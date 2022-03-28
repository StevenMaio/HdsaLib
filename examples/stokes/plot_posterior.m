clear
close all
clc

write_to_file = false;

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

m = min([opt_z_mean(1:3:N);control_read(1:3:N)]);
M = max([opt_z_mean(1:3:N);control_read(1:3:N)]);

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
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_std(1:3:N));
shading interp;
view(0,90)
title('Updated controller $z_x$ standard deviation','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'x_velocity_control_std','epsc')
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

m = min([opt_z_mean(2:3:N);control_read(2:3:N)]);
M = max([opt_z_mean(2:3:N);control_read(2:3:N)]);

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
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_std(2:3:N));
shading interp;
view(0,90)
title('Updated controller $z_y$ standard deviation','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'y_velocity_control_std','epsc')
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

m = min([hifi_state(2:3:N);updated_hifi_state(2:3:N)]);
M = max([hifi_state(2:3:N);updated_hifi_state(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High-fidelity state $v_y(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'nominal_hifi_state','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), updated_hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High-fidelity state $v_y(\overline{z}+F_\theta(0)\overline{\theta})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'updated_hifi_state','epsc')
end

figure,


if write_to_file 
    cd ../
end
