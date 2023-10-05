clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
load('z_prior_samples.txt')
load('delta_prior_samples.txt')

S = 5;
for k = 1:S
    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), z_prior_samples(:,k));
    shading interp;
    colorbar();
    view(0,90)
    axis square
    title('Prior Control Sample')

    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), delta_prior_samples(:,k));
    shading interp;
    colorbar();
    view(0,90)
    axis square
    title('Prior Discrepancy Sample')
end