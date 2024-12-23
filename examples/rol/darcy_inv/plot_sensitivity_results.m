clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('optimal_z.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z = data_obj.data;
data_obj = importdata('optimal_z_perturbed.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z_pert = data_obj.data;
data_obj = importdata('z_star_fe.txt', ' ', 2);  %% we need to skip the first two lines
estimate_z_fe = data_obj.data;

c_min = min([estimate_z;estimate_z_pert;estimate_z_fe]);
c_max = max([estimate_z;estimate_z_pert;estimate_z_fe]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z);
shading interp;
clim([c_min,c_max])
view(0,90)
colorbar
axis square
title('Nominal Optimal Solution')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z_pert);
shading interp;
clim([c_min,c_max])
view(0,90)
colorbar
axis square
title('Perturbed Optimal Solution')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), estimate_z_fe);
shading interp;
clim([c_min,c_max])
view(0,90)
colorbar
axis square
title('FE Perturbed Optimal Solution')