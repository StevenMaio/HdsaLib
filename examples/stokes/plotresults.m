%%
clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('uncontrolled_state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_uncontrolled_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Uncontrolled x-velocity')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Uncontrolled y-velocity')

% figure,
% trisurf(adj, nodes(:,1), nodes(:,2), P);
% shading interp;
% view(0,90)
% axis equal
% axis tight
% colorbar()
% title('Uncontrolled pressure')
% 
% figure,
% quiver(nodes(:,1), nodes(:,2), Ux, Uy);
% axis equal
% axis tight
% title('Uncontrolled velocity field')

%%
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
data_obj = importdata('map_state.txt', ' ', 9);  %% we need to skip the first 9 lines
map_state = data_obj.data;
map_state = map_state(1:2:end)+1;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Controlled x-velocity')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Controlled y-velocity')

% figure,
% trisurf(adj, nodes(:,1), nodes(:,2), P);
% shading interp;
% view(0,90)
% axis equal
% axis tight
% colorbar()
% title('Controlled pressure')
% 
% figure,
% quiver(nodes(:,1), nodes(:,2), Ux, Uy);
% axis equal
% axis tight
% title('Controlled velocity field')

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
data_obj = importdata('map_control.txt', ' ', 9);  %% we need to skip the first 9 lines
map_control = data_obj.data;
map_control = map_control(1:2:end)+1;
[tmp, control_permute] = sort(map_control);
control = control(control_permute);  %% we need to permute the control according to parallel maps

% Extract x-velocity
Zx = control(1:3:N);
% Extract y-velocity
Zy = control(2:3:N);
% Extract pressure
P  = control(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zx);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('X Controller')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Zy);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Y Controller')

%%
data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
[tmp, state_permute] = sort(map_state);
state = state(state_permute);  %% we need to permute the state according to parallel maps

N = 3*length(nodes);
% Extract x-velocity
Ux = state(1:3:N);
% Extract y-velocity
Uy = state(2:3:N);
% Extract pressure
P  = state(3:3:N);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Ux);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Controlled high-fidelity x-velocity')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Uy);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Controlled high-fidelity y-velocity')