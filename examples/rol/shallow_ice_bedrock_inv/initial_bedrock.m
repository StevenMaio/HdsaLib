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

% plot to see smoothing 
fig1 =figure(1);
trisurf(adj, nodes(:,1), nodes(:,2),b_initial)
shading interp;
colorbar
view(0,90);
title('Before Smoothing')
fname = '/ascldap/users/wmreese/Documents/Experiment_Data/ip_smoothed_T_2_nt_25_alpha_1e-4_greg_1e-3_freg_1e-7/before_smooth';
%print(fig1,fname,'-depsc')

smooth_b_initial = L_norm*b_initial;
fig2 = figure(2);
trisurf(adj, nodes(:,1), nodes(:,2),smooth_b_initial)
shading interp;
colorbar
view(0,90);
title('After Smoothing')
fname = '/ascldap/users/wmreese/Documents/Experiment_Data/ip_smoothed_T_2_nt_25_alpha_1e-4_greg_1e-3_freg_1e-7/after_smooth';
%print(fig2,fname,'-depsc')


fileID = fopen('matlab_initial_iterate.txt', 'w');
fprintf(fileID,'%12.8f\n', smooth_b_initial);
fclose(fileID);
