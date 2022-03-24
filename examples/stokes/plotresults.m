clear
close all
clc

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

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
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

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
ZP  = control(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_Ux);
shading interp;
view(0,90)
title('Uncontrolled X-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_Uy);
shading interp;
view(0,90)
title('Uncontrolled Y-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncon_P);
shading interp;
view(0,90)
title('Uncontrolled Pressure')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
title('X-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
title('Y-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), P);
shading interp;
view(0,90)
title('Pressure')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_Ux);
shading interp;
view(0,90)
title('HiFi X-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_Uy);
shading interp;
view(0,90)
title('HiFi Y-Velocity')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_P);
shading interp;
view(0,90)
title('HiFi Pressure')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zx);
shading interp;
view(0,90)
title('X-Velocity Control')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zy);
shading interp;
view(0,90)
title('Y-Velocity Control')
colorbar(); axis equal
axis tight