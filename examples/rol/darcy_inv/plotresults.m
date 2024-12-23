clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('optimal_u.txt', ' ', 2);  %% we need to skip the first two lines
estimate_u = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_u);
shading interp;
view(0,90)
colorbar
axis square
title('State')

data_obj = importdata('optimal_z.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z);
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability Estimate')

data_obj = importdata('initial_iterate.txt', ' ', 2);  %% we need to skip the first two lines
initial_z = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), initial_z);
shading interp;
view(0,90)
colorbar
axis square
title('Initial Log Permeability')

data_obj = importdata('Data_Generation/true_z.txt', ' ', 2);  %% we need to skip the first two lines
true_z = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_z);
shading interp;
view(0,90)
colorbar
axis square
title('True Log Permeability')
