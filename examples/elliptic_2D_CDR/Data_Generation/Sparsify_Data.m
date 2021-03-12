clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
state = load('true_state.txt');

state_sparse = reshape(state,sqrt(length(state)),sqrt(length(state)));
state_sparse = state_sparse(1:2:end,1:2:end);
state_sparse = state_sparse(:);
state_sparse = state_sparse.*(1 + 0.02*randn(length(state_sparse),1));

cd ..
writematrix(state_sparse,'data.txt','Delimiter','tab');