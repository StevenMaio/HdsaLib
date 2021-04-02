clear
close all
clc

mu = load('Prior_Mean.txt');
M = load('Mass_Mat.txt');
A = load('Elliptic_Operator.txt');

% alpha = .1;
% gamma = 1;
% K = (A-alpha*M)/gamma;
% 
% alpha = .5;
% gamma = 1;
% A = gamma*K + alpha*M;

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates

M_sqrt = sqrtm(M);

S = 5;
n = size(A,1);
X = zeros(n,S);

for k = 1:S
   X(:,k) = mu + linsolve(A,M_sqrt*randn(n,1)); 
end

for k = 1:S
    figure,
    trisurf(adj, nodes(:,1), nodes(:,2), X(:,k));
    shading interp;
    colorbar();
    view(0,90)
    axis square
    title('Prior Sample')
end
