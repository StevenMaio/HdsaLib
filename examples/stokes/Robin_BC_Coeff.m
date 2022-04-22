%% Control prior
clear
close all
clc

% -gamma*Lap + alpha*I
gamma = 10^-4;
alpha = 1;
m = 2000;

kappa = sqrt(alpha/gamma);

[X,Y] = meshgrid(linspace(0,1,m)');
W = 0*X + 1;
W(1,:) = .5*W(1,:);
W(end,:) = .5*W(end,:);
W(:,1) = .5*W(:,1);
W(:,end) = .5*W(:,end);
W = W/sum(sum(W));

% x=0 boundary
x_coord = 0;
y_nodes = linspace(0,1,500)';
y_nodes = y_nodes(2:end-1);
n = [-1 ; 0];

L = length(y_nodes);
beta_tilde = zeros(L,1);

for k = 1:L
    y = [x_coord ; y_nodes(k)];
    
    r = sqrt( (X-y(1)).^2 + (Y-y(2)).^2 );
    bess0 = besselk(0,kappa*r);
    bess1 = besselk(1,kappa*r);
    z = (y(1)-X)*n(1) + (y(2)-Y)*n(2);
    
    beta_tilde(k) = (1/2)*kappa*sum(sum((bess0.^2+bess1.^2).*z.*W))/sum(sum(r.*bess1.*bess0.*W));
end

beta = max(0,beta_tilde);
y_nodes = [0; y_nodes; 1];
beta = [beta(1);beta];
beta = [beta;beta(end)];

figure,
plot(y_nodes,beta)

writematrix(y_nodes,'prior_control_robin_coeff_input.txt')
writematrix(beta,'prior_control_robin_coeff.txt')

%% Discrepancy prior
clear
clc

% -gamma*Lap + alpha*I
gamma = 10^-2;
alpha = .05;
m = 2000;

kappa = sqrt(alpha/gamma);

[X,Y] = meshgrid(linspace(0,1,m)');
W = 0*X + 1;
W(1,:) = .5*W(1,:);
W(end,:) = .5*W(end,:);
W(:,1) = .5*W(:,1);
W(:,end) = .5*W(:,end);
W = W/sum(sum(W));

% x=0 boundary
x_coord = 0;
y_nodes = linspace(0,1,500)';
y_nodes = y_nodes(2:end-1);
n = [-1 ; 0];

L = length(y_nodes);
beta_tilde = zeros(L,1);

for k = 1:L
    y = [x_coord ; y_nodes(k)];
    
    r = sqrt( (X-y(1)).^2 + (Y-y(2)).^2 );
    bess0 = besselk(0,kappa*r);
    bess1 = besselk(1,kappa*r);
    z = (y(1)-X)*n(1) + (y(2)-Y)*n(2);
    
    beta_tilde(k) = (1/2)*kappa*sum(sum((bess0.^2+bess1.^2).*z.*W))/sum(sum(r.*bess1.*bess0.*W));
end

beta = max(0,beta_tilde);
y_nodes = [0; y_nodes; 1];
beta = [beta(1);beta];
beta = [beta;beta(end)];

figure,
plot(y_nodes,beta)

writematrix(y_nodes,'prior_discrepancy_robin_coeff_input.txt')
writematrix(beta,'prior_discrepancy_robin_coeff.txt')
