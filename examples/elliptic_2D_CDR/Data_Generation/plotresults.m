clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

state = load('true_state.txt');
source = load('true_source.txt');

axsize = 200;
figure('Position', [100 680 3*axsize 2*axsize]);
figure('Position', [1000 680 3*axsize 2*axsize]);

figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
colorbar
axis square
title('True State')

figure(2)
trisurf(adj, nodes(:,1), nodes(:,2), source);
shading interp;
view(0,90)
colorbar
axis square
title('True Source')
