clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
N = 3*length(nodes);

load('opt_z_mean.txt')
load('opt_z_samples.txt')
load('delta_mean_at_z1.txt')
load('delta_samples_at_z1.txt')
load('Y.txt')
load('control_read.txt')

data_obj = importdata('hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
hifi_state = data_obj.data;
data_obj = importdata('updated_hifi_state.txt', ' ', 2);  %% we need to skip the first two lines
updated_hifi_state = data_obj.data;

delta_at_z1_std = std(delta_samples_at_z1,[],2);
opt_z_std = std(opt_z_samples,[],2);

m = min([delta_mean_at_z1(1:3:N);Y(1:3:N)]);
M = max([delta_mean_at_z1(1:3:N);Y(1:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('X-Velocity Mean Discrepancy')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('X-Velocity Discrepancy Data')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(1:3:N));
shading interp;
view(0,90)
title('X-Velocity Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([delta_mean_at_z1(2:3:N);Y(2:3:N)]);
M = max([delta_mean_at_z1(2:3:N);Y(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Y-Velocity Mean Discrepancy')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Y-Velocity Discrepancy Data')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(2:3:N));
shading interp;
view(0,90)
title('Y-Velocity Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([delta_mean_at_z1(3:3:N);Y(3:3:N)]);
M = max([delta_mean_at_z1(3:3:N);Y(3:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_mean_at_z1(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Pressure Mean Discrepancy')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), Y(3:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Pressure Discrepancy Data')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_at_z1_std(3:3:N));
shading interp;
view(0,90)
title('Pressure Discrepancy Standard Deviation')
colorbar(); axis equal
axis tight

m = min([opt_z_mean(1:3:N);control_read(1:3:N)]);
M = max([opt_z_mean(1:3:N);control_read(1:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_mean(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated X-Velocity Control')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_std(1:3:N));
shading interp;
view(0,90)
title('X-Velocity Controller Standard Deviation')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), control_read(1:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal X-Velocity Control')
colorbar(); axis equal
axis tight

m = min([opt_z_mean(2:3:N);control_read(2:3:N)]);
M = max([opt_z_mean(2:3:N);control_read(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_mean(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Updated Y-Velocity Control')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_std(2:3:N));
shading interp;
view(0,90)
title('Y-Velocity Controller Standard Deviation')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), control_read(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('Nominal Y-Velocity Control')
colorbar(); axis equal
axis tight

m = min([hifi_state(2:3:N);updated_hifi_state(2:3:N)]);
M = max([hifi_state(2:3:N);updated_hifi_state(2:3:N)]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High Fidelity Y-Velocity at Nominal Control')
colorbar(); axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), updated_hifi_state(2:3:N));
shading interp;
view(0,90)
caxis([m,M])
title('High Fidelity Y-Velocity at Updated Control')
colorbar(); axis equal
axis tight
