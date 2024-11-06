clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

true_z = 5*(1+nodes(:,2)).*sin(pi*((nodes(:,1)-.2*nodes(:,2)))).^2 - 1;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_z);
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability')

writematrix(true_z,'Log_Permability.txt','Delimiter',' ')