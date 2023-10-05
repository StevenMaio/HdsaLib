clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

load('delta_mean_at_z1.txt')
load('delta_samples_at_z1.txt')
load('Y.txt')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(:,1));
shading interp;
colorbar();
view(0,90)
axis square
title('Discrepancy Data')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1);
shading interp;
colorbar();
view(0,90)
axis square
title('Discrepancy Fit')

S = 5;
for k = 1:S
    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), delta_samples_at_z1(:,k));
    shading interp;
    colorbar();
    view(0,90)
    axis square
    title('Discrepancy Sample')
end