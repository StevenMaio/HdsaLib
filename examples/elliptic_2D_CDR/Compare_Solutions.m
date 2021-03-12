clear
close all
clc

% With reaction
cd with_reaction/

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state_with = data_obj.data;
data_obj = importdata('source.txt', ' ', 2);  %% we need to skip the first two lines
source_with = data_obj.data;

% without reaction
cd ../without_reaction
data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state_without = data_obj.data;
data_obj = importdata('source.txt', ' ', 2);  %% we need to skip the first two lines
source_without = data_obj.data;

noisy_state = load('data.txt');

cd ..

true_state = load('Data_Generation/true_state.txt');
true_source = load('Data_Generation/true_source.txt');
true_adj = load('Data_Generation/cell_to_node_quad.txt') + 1;
true_nodes = load('Data_Generation/nodes.txt');

axsize = 200;
figure('Position', [10 680 3*axsize 2*axsize]);
figure('Position', [650 680 3*axsize 2*axsize]);
figure('Position', [1300 680 3*axsize 2*axsize]);
figure('Position', [10 30 3*axsize 2*axsize]);
figure('Position', [650 30 3*axsize 2*axsize]);
figure('Position', [1300 30 3*axsize 2*axsize]);

state_min = min([state_with;state_without;noisy_state]);
state_max = max([state_with;state_without;noisy_state]);

source_min = min([source_with;source_without;true_source]);
source_max = max([source_with;source_without;true_source]);

figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), state_with);
shading interp;
colorbar();
caxis([state_min,state_max])
view(0,90)
axis square
title('State with reaction')

figure(2)
trisurf(adj, nodes(:,1), nodes(:,2), state_without);
shading interp;
colorbar();
caxis([state_min,state_max])
view(0,90)
axis square
title('State without reaction')

figure(3)
trisurf(adj, nodes(:,1), nodes(:,2), noisy_state);
shading interp;
colorbar();
caxis([state_min,state_max])
view(0,90)
axis square
title('State data')

figure(4)
trisurf(adj, nodes(:,1), nodes(:,2), source_with);
shading interp;
colorbar();
caxis([source_min,source_max])
view(0,90)
axis square
title('Source with reaction')

figure(5)
trisurf(adj, nodes(:,1), nodes(:,2), source_without);
shading interp;
colorbar();
caxis([source_min,source_max])
view(0,90)
axis square
title('Source without reaction')

figure(6)
trisurf(true_adj, true_nodes(:,1), true_nodes(:,2), true_source);
shading interp;
colorbar();
caxis([source_min,source_max])
view(0,90)
axis square
title('True Source')

