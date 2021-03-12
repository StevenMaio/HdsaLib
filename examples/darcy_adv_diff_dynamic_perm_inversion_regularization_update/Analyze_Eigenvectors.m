clear
close all
clc

eval = load('LIS_evals_1.txt');
evecs = load('LIS_evecs_1.txt');

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

n = size(evecs,2);
for k = 1:n
    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), evecs(1:2:end,k));
    shading interp;
    view(2)
    colorbar
    axis square;
    title(['Eigenvector number ',num2str(k),' with eigenvalue \lambda=',num2str(eval(k))])
end