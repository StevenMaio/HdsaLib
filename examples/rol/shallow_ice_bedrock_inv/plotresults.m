clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('Data_Generation/true_bed.txt', ' ', 2);  %% we need to skip the first two lines
beta = data_obj.data(1:3:end);
figure,
trisurf(adj, nodes(:,1), nodes(:,2), beta);
shading interp;
view(0,90)
colorbar
axis square
title('True Parameter')

data_obj = importdata('optimal_z.txt', ' ', 2);  %% we need to skip the first two lines
beta = data_obj.data(1:3:end);
figure,
trisurf(adj, nodes(:,1), nodes(:,2), beta);
shading interp;
view(0,90)
colorbar
axis square
title('Estimated Parameter')