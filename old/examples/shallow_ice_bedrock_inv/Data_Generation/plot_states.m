% script to plot s, vx, vy over time
clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

load clean_true_state.txt
%load coarse_data.txt
%clean_true_state = coarse_data;
nt = 21; 
%nt = 41;
%T = 0.5;
T = 1.0;
%T = 10.0;
Tmesh   = linspace(0,T,nt);

% axsize = 200;
% figure('Position', [100 680 3*axsize 2*axsize]);
% figure('Position', [1000 680 3*axsize 2*axsize]);
% figure('Position', [100 50 3*axsize 2*axsize]);
% figure('Position', [1000 50 3*axsize 2*axsize]);
% 
% data_obj = importdata('true_bed.txt', ' ', 2);  %% we need to skip the first two lines
% beta = data_obj.data(1:3:end);
% figure(1)
% trisurf(adj, nodes(:,1), nodes(:,2), beta);
% shading interp;
% view(0,90)
% colorbar
% axis square
% title('True Parameter')


times = [floor(nt/4) floor(nt/2) floor(nt*3/4) nt];
for i = times

  state = clean_true_state(i,1:3:end);

  %figure(2)
  fig1 = figure;
  trisurf(adj, nodes(:,1), nodes(:,2), state);
  shading interp;
  view(2);
  axis square;
  title(['True state at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  fname = ['/ascldap/users/wmreese/Documents/State_Snapshots_GIS/51x51_mesh_T_1_nt_21/s_',num2str(i)];
  print(fig1,fname,'-depsc')
  %drawnow
  
  xvel = clean_true_state(i,2:3:end); 
  %figure(3)
  fig2 = figure;
  trisurf(adj, nodes(:,1), nodes(:,2), xvel);
  shading interp;
  view(2);
  axis square;
  title(['True x velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  fname = ['/ascldap/users/wmreese/Documents/State_Snapshots_GIS/51x51_mesh_T_1_nt_21/vx_',num2str(i)];
  print(fig2,fname,'-depsc')
  %drawnow
  
  yvel = clean_true_state(i,3:3:end); 
  %figure(4)
  fig3 = figure;
  trisurf(adj, nodes(:,1), nodes(:,2), yvel);
  shading interp;
  view(2);
  axis square;
  title(['True y velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  fname = ['/ascldap/users/wmreese/Documents/State_Snapshots_GIS/51x51_mesh_T_1_nt_21/vy_',num2str(i)];
  print(fig3,fname,'-depsc')
  %drawnow
  
 % pause(.5);
end

