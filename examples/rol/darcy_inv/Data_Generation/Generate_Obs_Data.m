clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('true_u.txt', ' ', 2);  %% we need to skip the first two lines
true_u = data_obj.data;

obs_locations = [];
xobs = .1:.2:.9;
yobs = .1:.2:.9;

noise_level = 0.02;

for x = xobs
   for y = yobs
    pt = intersect(find(abs(nodes(:,1)-x)<1.e-8),find(abs(nodes(:,2)-y)<1.e-8));
    obs_locations = [obs_locations,pt];
   end
end

obs_data = true_u(obs_locations).*(1 + noise_level*randn(length(obs_locations),1));

figure,
trisurf(adj, nodes(:,1), nodes(:,2), true_u);
shading interp;
view(0,90)
colorbar
axis square
title('State')
hold on
for k = 1:length(obs_locations)
    scatter3(nodes(obs_locations(k),1),nodes(obs_locations(k),2),10,36,'black','filled')
end

m_coarse = 51;
x_coarse = linspace(0,1,m_coarse )';
y_coarse = linspace(0,1,m_coarse )';
nodes_coarse = zeros(m_coarse^2,2);
for i = 1:m_coarse
    I = (1+(i-1)*m_coarse):(i*m_coarse);
    nodes_coarse(I,1) = x_coarse;
    nodes_coarse(I,2) = y_coarse(i);
end

obs_locations = [];
for x = xobs
   for y = yobs
    pt = intersect(find(abs(nodes_coarse(:,1)-x)<1.e-8),find(abs(nodes_coarse(:,2)-y)<1.e-8));
    obs_locations = [obs_locations,pt];
   end
end


writematrix(obs_data,'obs_data.txt','Delimiter',' ')
writematrix(obs_locations,'obs_locations.txt','Delimiter',' ')