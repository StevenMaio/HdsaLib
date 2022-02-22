clear
close all
clc

m = 1;
n = 1;

grad_J_u = randn(1,m);
Mz = randn(n,1);

Sigma = [ (Mz'*Mz)*eye(m) , grad_J_u'*Mz' ; Mz*grad_J_u , (grad_J_u*grad_J_u')*eye(n) ];

N = 10000;

% Omega_u = zeros(m,N);
% Omega_z = zeros(n,N);
% for k = 1:N
%    Omega_z(:,k) = norm(grad_J_u)*randn(n,1);
%    Omega_u(:,k) = ((Mz'*Omega_z(:,k))/(grad_J_u*grad_J_u'))*grad_J_u' + norm(Mz)*(eye(m) - (1/(grad_J_u*grad_J_u'))*(grad_J_u'*grad_J_u))*randn(m,1);
% end

Omega_z = norm(grad_J_u)*randn(n,N);
Omega_u = grad_J_u'*((Mz'*Omega_z)/(grad_J_u*grad_J_u')) + norm(Mz)*(eye(m) - (1/(grad_J_u*grad_J_u'))*(grad_J_u'*grad_J_u))*randn(m,N);

Omega_u_ind = norm(Mz)*randn(m,N);
Omega_z_ind = norm(grad_J_u)*randn(n,N);

R = mvnrnd(zeros(m+n,1),Sigma,N);

Omega_u_R = R(:,1:m)';
Omega_z_R = R(:,(m+1):end)';

figure,
scatter(Omega_u,Omega_z,'o')

figure,
scatter(Omega_u_R,Omega_z_R)