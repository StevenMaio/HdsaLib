clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

load clean_true_state.txt
%load coarse_data.txt
%clean_true_state = coarse_data;
%nt = 21; 
%nt = 41;
nt = 121;
%T = 0.5;
%T = 0.75;
%T = 1.0;
%T = 2.0;
%T = 5.0;
T = 10.0;
Tmesh   = linspace(0,T,nt);

axsize = 200;
figure('Position', [100 680 3*axsize 2*axsize]);
figure('Position', [1000 680 3*axsize 2*axsize]);
figure('Position', [100 50 3*axsize 2*axsize]);
figure('Position', [1000 50 3*axsize 2*axsize]);

data_obj = importdata('true_bed.txt', ' ', 2);  %% we need to skip the first two lines
beta = data_obj.data(1:3:end);
figure(1)
trisurf(adj, nodes(:,1), nodes(:,2), beta);
shading interp;
view(0,90)
colorbar
axis square
title('True Parameter')
  
for i=1:nt

  state = clean_true_state(i,1:3:end);

  figure(2)
  H = state(:)-beta;
  trisurf(adj, nodes(:,1), nodes(:,2), H);
  shading interp;
  view(2);
  axis square;
  title(['True thickness at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
  xvel = clean_true_state(i,2:3:end); 
  figure(3)
  trisurf(adj, nodes(:,1), nodes(:,2), xvel);
  shading interp;
  view(2);
  axis square;
  title(['True x velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
  
  yvel = clean_true_state(i,3:3:end); 
  figure(4)
  trisurf(adj, nodes(:,1), nodes(:,2), yvel);
  shading interp;
  view(2);
  axis square;
  title(['True y velocity at time t = ',num2str(Tmesh(i))],'fontsize',16);
  xlabel('x');
  ylabel('y');
  colorbar
  set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
  drawnow
%   
% H = state(:)-beta;
% binH = 0*H;
% binH(H == 0) = 1 ;
% figure(1)
% hold on
%   trisurf(adj, nodes(:,1), nodes(:,2), binH);
%   k = find(nodes(:,1)< 23);
%   for j = 1:length(k)
%       scatter3(nodes(k(j),1),nodes(k(j),2),1.1*max(binH))
%   end
%   shading interp;
%   view(2);
%   axis square;
%   title(['Thickness at time t = ',num2str(Tmesh(i))],'fontsize',16);
%   xlabel('x');
%   ylabel('y');
%   colorbar
%   set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
%   drawnow
%   
%   figure(2)
%   trisurf(adj, nodes(:,1), nodes(:,2), min(0,H));
%   shading interp;
%   view(2);
%   axis square;
%   title(['Thickness at time t = ',num2str(Tmesh(i))],'fontsize',16);
%   xlabel('x');
%   ylabel('y');
%   colorbar
%   set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');% tightfig;
%   drawnow
  pause(.5);
  %pause;
  % This is a loop to see if the thickness goes below zero
   %H = state(:) - beta;
   
   
%   if (min(H) < 0)
%       pause 
%       min(H)
%       figure(5);
%       trisurf(adj,nodes(:,1),nodes(:,2),min(0,H));
%       view(0,90)
%       shading interp;
%       colorbar
%       fig = figure(5);
%   end
end

