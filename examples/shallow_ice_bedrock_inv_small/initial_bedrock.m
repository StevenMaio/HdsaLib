% read in the bedrock topography
% smooth it and spit out matlab_initial_iterate.txt

% load initial condition and nodes 
b_initial = load('Bedrock_Topography.txt');
nodes = load('nodes.txt');  %% load node coordinates
adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

% this is where we will smooth it
m = length(b_initial);
alpha = 1.e-4; % smoothing level; alpha = 5 gives I, less than 1e-6 - gives mat of 1's (too much smoothing)
x = nodes(:,1); y= nodes(:,2);
D = zeros(m,m);
for i = 1:m
    for j = 1:m
        D(i,j) = (x(i)-x(j))^2 + (y(i)-y(j))^2 ;
    end
end

L = exp(-alpha.*D);

% normalize so that we compute averages
L_norm = L./repmat(sum(L,2),1,size(L,2));
smooth_b_initial = L_norm*b_initial; 
m1 = min(b_initial);
m2 = min(smooth_b_initial);
m = min(m1,m2);
M1 = max(b_initial);
M2 = max(smooth_b_initial);
M = max(M1,M2);



% plot to see smoothing 
fig1 =figure(1);
trisurf(adj, nodes(:,1), nodes(:,2),b_initial)
shading interp;
colorbar
view(0,90);
caxis([m,M])
title('Before Smoothing')
fname = '/ascldap/users/wmreese/Documents/Experiment_Data/ip_small/before_smooth';
print(fig1,fname,'-depsc')

fig2 = figure(2);
trisurf(adj, nodes(:,1), nodes(:,2),smooth_b_initial)
shading interp;
colorbar
view(0,90);
caxis([m,M])
title('After Smoothing')
fname = '/ascldap/users/wmreese/Documents/Experiment_Data/ip_small/after_smooth';
print(fig2,fname,'-depsc')


fileID = fopen('matlab_initial_iterate.txt', 'w');
fprintf(fileID,'%12.8f\n', smooth_b_initial);
fclose(fileID);
