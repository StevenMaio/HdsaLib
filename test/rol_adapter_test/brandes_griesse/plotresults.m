clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
colorbar;
view(0,90)
axis square
title('State')

data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;
figure(2)
trisurf(adj, nodes(:,1), nodes(:,2), control);
shading interp;
colorbar;
view(0,90)
axis square
title('Control')
