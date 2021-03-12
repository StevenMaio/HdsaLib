clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

load clean_true_state.txt

axsize = 200;
figure('Position', [50 680 3*axsize 2*axsize]);
figure('Position', [700 680 3*axsize 2*axsize]);
figure('Position', [1350 680 3*axsize 2*axsize]);

data_obj = importdata('true_beta.txt', ' ', 2);  %% we need to skip the first two lines
beta = data_obj.data;
figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), beta(1:2:end));
shading interp;
view(0,90)
colorbar
axis square
title('Sliding Coefficient')
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
 
figure(2)
trisurf(adj, nodes(:,1), nodes(:,2), clean_true_state(1:2:end));
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
colorbar
axis square
title('x velocity')
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

figure(3)
trisurf(adj, nodes(:,1), nodes(:,2), clean_true_state(2:2:end));
shading interp;
view(2);
axis square;
xlabel('x');
ylabel('y');
colorbar
axis square
title('y velocity')
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;


