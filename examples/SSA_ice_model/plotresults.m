clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

active_sensors = load('Active_Sensors.txt');

n_p = 49; % nsy*(nsx-2)
n_c = 81; % nsy*nsx
n = n_p + n_c;
I_p = 1:(n_p);
I_c = (I_p(end)+1):(n_p+n_c);

nt = 24;
T = .16;
Tmesh   = linspace(0,T,nt+1);

load('data.txt');
pressure = data(2,1:2:end);

axsize = 200;
figure('Position', [10 680 3*axsize 2*axsize]);
figure('Position', [650 680 3*axsize 2*axsize]);
figure('Position', [1300 680 3*axsize 2*axsize]);
figure('Position', [10 30 3*axsize 2*axsize]);
figure('Position', [650 30 3*axsize 2*axsize]);
figure('Position', [1300 30 3*axsize 2*axsize]);

figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), pressure);
shading interp;
view(2);
axis square;
title('Observed Pressure','fontsize',16);
xlabel('x');
ylabel('y');
colorbar
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

data_state = importdata('state_1.txt', ' ', 2);  %% we need to skip the first two lines
pressure = data_state.data(1:2:end);
figure(4)
trisurf(adj, nodes(:,1), nodes(:,2), pressure);
shading interp;
view(2);
axis square;
title('Pressure','fontsize',16);
xlabel('x');
ylabel('y');
colorbar
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;

data_obj1 = importdata('perm.txt', ' ', 2);  %% we need to skip the first two lines
perm1 = data_obj1.data;
perm1 = perm1(1:2:end);
data_obj2 = importdata('Data_Generation/true_perm.txt', ' ', 2);  %% we need to skip the first two lines
perm2 = data_obj2.data;
perm2 = perm2(1:2:end);

figure(3)
trisurf(adj, nodes(:,1), nodes(:,2), perm1);
shading interp;
view(0,90)
colorbar
axis square
title('Log Permeability Field')

adj_f = load('Data_Generation/cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes_f = load('Data_Generation/nodes.txt');  %% load node coordinates
figure(6)
trisurf(adj_f, nodes_f(:,1), nodes_f(:,2), perm2);
shading interp;
view(0,90)
colorbar
axis square
title('True Log Permeability Field')

for i=1:nt-1

  
  contaminant = data(i,2:2:end);
  figure(2)
  trisurf(adj, nodes(:,1), nodes(:,2), contaminant);
  shading interp;
  view(2);
  axis square;
  title(['Observed contaminant at time t = ',num2str(Tmesh(i+1))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
  data_state = importdata(['state_',int2str(i),'.txt'], ' ', 2);  %% we need to skip the first two lines
  contaminant = data_state.data(2:2:end);
  figure(5)
  trisurf(adj, nodes(:,1), nodes(:,2), contaminant);
  shading interp;
  view(2);
  axis square;
  title(['Contaminant at time t = ',num2str(Tmesh(i+1))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow

end

