clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
colorbar();
view(0,90)
axis square
title('State')

data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
hifi_state = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state);
shading interp;
colorbar();
view(0,90)
axis square
title('High Fidelity State')

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
colorbar();
view(0,90)
axis square
title('Control')

