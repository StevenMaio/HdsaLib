clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

%true_z = 5*(1+nodes(:,2)).*sin(pi*((nodes(:,1)-.1*nodes(:,2)))).^2 - 1;
x = nodes(:,1);
y = nodes(:,2);
true_z =  9*exp(-10*(x-.5).^2).*exp(-20*(y-.5).^2).*exp(-20*(x-.5).*(y-.5));
true_z =  true_z + 3*exp(-10*(x-.7).^2).*exp(-20*(y-.7).^2).*exp(-20*(x-.7).*(y-.7));
true_z =  true_z + 4*exp(-10*(x-.3).^2).*exp(-20*(y-.3).^2).*exp(-20*(x-.3).*(y-.3));
true_z = true_z - 4;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_z);
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability')

writematrix(true_z,'Log_Permability.txt','Delimiter',' ')