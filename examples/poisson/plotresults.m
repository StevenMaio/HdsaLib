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

true_state = nodes(:,1).*(1-nodes(:,1)).*nodes(:,2).*(1-nodes(:,2));
true_source = 2*(nodes(:,1)-nodes(:,1).^2 + nodes(:,2) - nodes(:,2).^2);

disp('Relative State Error')
disp(norm(true_state-state)/norm(true_state))
disp('Relative Source Error')
disp(norm(true_source-source)/norm(true_source))

figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_state);
shading interp;
colorbar();
view(0,90)
axis square
title('True State')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_source);
shading interp;
colorbar();
view(0,90)
axis square
title('True Source')