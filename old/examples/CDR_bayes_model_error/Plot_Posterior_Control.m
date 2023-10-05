clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

load('opt_z_mean.txt')
load('opt_z_samples.txt')
load('control_read.txt')

lower = min([min(control_read),min(opt_z_mean),min(min(opt_z_samples))]);
upper = max([max(control_read),max(opt_z_mean),max(max(opt_z_samples))]);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), control_read);
shading interp;
colorbar();
caxis([lower,upper])
view(0,90)
axis square
title('Nominal Control')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), opt_z_mean);
shading interp;
colorbar();
caxis([lower,upper])
view(0,90)
axis square
title('Updated Control Mean')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), std(opt_z_samples,[],2));
shading interp;
colorbar();
view(0,90)
axis square
title('Updated Control Standard Deviation')

S = 5;
for k = 1:S
    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), opt_z_samples(:,k));
    shading interp;
    caxis([lower,upper])
    colorbar();
    view(0,90)
    axis square
    title('Updated Control Sample')
end