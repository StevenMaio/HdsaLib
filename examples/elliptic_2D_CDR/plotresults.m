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

data_obj = importdata('source.txt', ' ', 2);  %% we need to skip the first two lines
source = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), source);
shading interp;
colorbar();
view(0,90)
axis square
title('Source')

noisy_state = load('data.txt');
figure,
trisurf(adj, nodes(:,1), nodes(:,2), noisy_state);
shading interp;
colorbar();
view(0,90)
axis square
title('Noisy State')

true_state = load('Data_Generation/true_state.txt');
true_source = load('Data_Generation/true_source.txt');
true_adj = load('Data_Generation/cell_to_node_quad.txt') + 1;
true_nodes = load('Data_Generation/nodes.txt');

figure,
trisurf(true_adj, true_nodes(:,1), true_nodes(:,2), true_state);
shading interp;
colorbar();
view(0,90)
axis square
title('True State')

figure,
trisurf(true_adj, true_nodes(:,1), true_nodes(:,2), true_source);
shading interp;
colorbar();
view(0,90)
axis square
title('True Source')