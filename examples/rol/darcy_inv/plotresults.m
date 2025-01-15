clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

adj_fine = load('Data_Generation/cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes_fine = load('Data_Generation/nodes.txt');  %% load node coordinates

data_obj = importdata('optimal_u.txt', ' ', 2);  %% we need to skip the first two lines
estimate_u = data_obj.data;
data_obj = importdata('optimal_z.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z = data_obj.data;
data_obj = importdata('optimal_u_perturbed.txt', ' ', 2);  %% we need to skip the first two lines
estimate_u_perturbed = data_obj.data;
data_obj = importdata('optimal_z_perturbed.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z_perturbed = data_obj.data;
data_obj = importdata('initial_iterate.txt', ' ', 2);  %% we need to skip the first two lines
initial_z = data_obj.data;
data_obj = importdata('Data_Generation/true_z.txt', ' ', 2);  %% we need to skip the first two lines
true_z = data_obj.data;

z_min = min([true_z;estimate_z;estimate_z_perturbed;initial_z]);
z_max = max([true_z;estimate_z;estimate_z_perturbed;initial_z]);
u_min = min([estimate_u;estimate_u_perturbed]);
u_max = max([estimate_u;estimate_u_perturbed]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_u);
clim([u_min,u_max])
shading interp;
view(0,90)
colorbar
axis square
title('State')


figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z);
clim([z_min,z_max])
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability Estimate')


figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_u_perturbed);
clim([u_min,u_max])
shading interp;
view(0,90)
colorbar
axis square
title('Reoptimized State')


figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z_perturbed);
clim([z_min,z_max])
shading interp;
view(0,90)
colorbar
axis square
title('Reoptimized Log Permeability Estimate')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), initial_z);
clim([z_min,z_max])
shading interp;
view(0,90)
colorbar
axis square
title('Initial Log Permeability')

figure,
trisurf(adj_fine, nodes_fine(:,1), nodes_fine(:,2), true_z);
clim([z_min,z_max])
shading interp;
view(0,90)
colorbar
axis square
title('True Log Permeability')
