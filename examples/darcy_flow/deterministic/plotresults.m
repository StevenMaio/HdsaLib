clear
clc
close all

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

target = load('Target_Data.txt');

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
view(0,90)
colorbar
axis square
title('State')
hold on 
for k = 1:size(target,1)
   scatter3(nodes(target(k,1)+1,1),nodes(target(k,1)+1,2),1.1*max(state),40,'filled','MarkerFaceColor','black')
end

data_obj = importdata('uncontrolled_state.txt', ' ', 2);  %% we need to skip the first two lines
uncontrolled_state = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), uncontrolled_state);
shading interp;
view(0,90)
colorbar
axis square
title('Uncontrolled State')
hold on 
for k = 1:size(target,1)
   scatter3(nodes(target(k,1)+1,1),nodes(target(k,1)+1,2),1.1*max(uncontrolled_state),40,'filled','MarkerFaceColor','black')
end

D = [uncontrolled_state(target(:,1))'; state(target(:,1))';target(:,2)']';
figure,
bar(D)
hold on
xlim([0,13])

load('control.txt')
xc = [.2, .2, .5, .5, .8, .8];
yc = [.5, .8, .5, .8, .5, .8];
c = 0;
for i = 1:length(xc)
   c = c + control(i)*1000*exp(-1000*( (nodes(:,1)-xc(i)).^2 + (nodes(:,2)-yc(i)).^2 ) ); 
end
figure,
trisurf(adj, nodes(:,1), nodes(:,2), c);
shading interp;
view(0,90)
colorbar
axis square
title('Control')
