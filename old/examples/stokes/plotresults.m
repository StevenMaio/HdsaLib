clear
close all
clc

write_to_file = true;

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('uncontrolled_state.txt', ' ', 2);  %% we need to skip the first two lines
uncon_state = data_obj.data;
N = 3*length(nodes);
% Extract x-velocity
uncon_Ux = uncon_state(1:3:N);
% Extract y-velocity
uncon_Uy = uncon_state(2:3:N);
% Extract pressure
uncon_P  = uncon_state(3:3:N);

state = load('state_read.txt');
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
hifi_state = data_obj.data;
% Extract x-velocity
hifi_Ux = hifi_state(1:3:N);
% Extract y-velocity
hifi_Uy = hifi_state(2:3:N);
% Extract pressure
hifi_P  = hifi_state(3:3:N);

control = load('control_read.txt');
% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
ZP  = control(3:3:N);

if write_to_file 
    cd figures
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_Ux);
shading interp;
view(0,90)
title('Uncontrolled $x$-velocity $\tilde{v}_x(0)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'uncontrolled_X_velocity','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_Uy);
shading interp;
view(0,90)
title('Uncontrolled $y$-velocity $\tilde{v}_y(0)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'uncontrolled_Y_velocity','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_P);
shading interp;
view(0,90)
title('Uncontrolled pressure $\tilde{p}(0)$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'uncontrolled_pressure','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
title('Controlled $x$-velocity $\tilde{v}_x(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'controlled_X_velocity','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
title('Controlled $y$-velocity $\tilde{v}_y(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'controlled_Y_velocity','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
title('Controlled pressure $\tilde{p}(\overline{z})$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)
if write_to_file
    saveas(gcf,'controlled_pressure','epsc')
end

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_Ux);
shading interp;
view(0,90)
title('High-fidelity $x$-velocity $v_x$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_Uy);
shading interp;
view(0,90)
title('High-fidelity $y$-velocity $v_y$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_P);
shading interp;
view(0,90)
title('High-fidelity pressure $p$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zx);
shading interp;
view(0,90)
title('$x$-velocity controller $z_x$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zy);
shading interp;
view(0,90)
title('$y$-velocity controller $z_y$','Interpreter','latex')
colorbar(); axis equal
axis tight
set(gca,'FontSize',18)

if write_to_file
    cd ../
end