clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

active_sensors = load('Active_Sensors.txt');

nt = 31;
T = 10;
Tmesh   = linspace(0,T,nt+1);

load('data.txt');

axsize = 200;
figure('Position', [10 680 3*axsize 2*axsize]);
figure('Position', [650 680 3*axsize 2*axsize]);
figure('Position', [1300 680 3*axsize 2*axsize]);
figure('Position', [10 30 3*axsize 2*axsize]);
figure('Position', [650 30 3*axsize 2*axsize]);
figure('Position', [1300 30 3*axsize 2*axsize]);

data_obj = importdata('beta.txt', ' ', 2);  %% we need to skip the first two lines
beta = data_obj.data;

figure(4)
trisurf(adj, nodes(:,1), nodes(:,2), beta(1:3:end));
shading interp;
view(0,90)
colorbar
axis square
title('Sliding Coefficient')

adj_f = load('Data_Generation/cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes_f = load('Data_Generation/nodes.txt');  %% load node coordinates
data_obj = importdata('Data_Generation/true_beta.txt', ' ', 2);  %% we need to skip the first two lines
beta_f = data_obj.data;
figure(1)
trisurf(adj_f, nodes_f(:,1), nodes_f(:,2), beta_f(1:3:end));
shading interp;
view(0,90)
colorbar
axis square
title('True Sliding Coefficient')

for i=1:nt

  figure(2)
  trisurf(adj, nodes(:,1), nodes(:,2), data(i,2:3:end));
  shading interp;
  view(2);
  axis square;
  title(['Observed x-velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
  figure(3)
  trisurf(adj, nodes(:,1), nodes(:,2), data(i,3:3:end));
  shading interp;
  view(2);
  axis square;
  title(['Observed y-velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
  data_state = importdata(['state_',int2str(i-1),'.txt'], ' ', 2);  %% we need to skip the first two lines
  xvel = data_state.data(2:3:end);
  yvel = data_state.data(3:3:end);
  figure(5)
  trisurf(adj, nodes(:,1), nodes(:,2), xvel);
  shading interp;
  view(2);
  axis square;
  title(['x-velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  figure(6)
  trisurf(adj, nodes(:,1), nodes(:,2), yvel);
  shading interp;
  view(2);
  axis square;
  title(['y-velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow

end

