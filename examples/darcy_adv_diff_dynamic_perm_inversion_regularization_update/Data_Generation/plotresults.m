clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

state = load('clean_true_state.txt');

nt = 151;
T = .5;
Tmesh   = linspace(0,T,nt);

axsize = 200;
figure('Position', [100 680 3*axsize 2*axsize]);
figure('Position', [1000 680 3*axsize 2*axsize]);
figure('Position', [100 30 3*axsize 2*axsize]);

pressure = state(2,1:2:end);
figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), pressure);
shading interp;
view(2);
axis square;
title('True pressure','fontsize',16);
xlabel('x');
ylabel('y');
colorbar
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
drawnow

data_obj = importdata('true_perm.txt', ' ', 2);  %% we need to skip the first two lines
perm = data_obj.data;
perm = perm(1:2:end);
figure(3)
trisurf(adj, nodes(:,1), nodes(:,2), perm);
shading interp;
view(2)
colorbar
axis square
title('True Log Permeability Field')
  
for i=1:nt

  contaminant = state(i,2:2:end);

  figure(2)
  trisurf(adj, nodes(:,1), nodes(:,2), contaminant);
  shading interp;
  view(2);
  axis square;
  title(['True contaminant at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
end


