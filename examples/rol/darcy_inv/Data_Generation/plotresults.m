clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('true_u.txt', ' ', 2);  %% we need to skip the first two lines
true_u = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_u);
shading interp;
view(0,90)
colorbar
axis square
title('State')

data_obj = importdata('true_z.txt', ' ', 2);  %% we need to skip the first two lines
true_z = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_z);
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability')
