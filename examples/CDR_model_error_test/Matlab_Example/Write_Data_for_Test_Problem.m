clear
clc
close all
rng(1234)

n_mesh = 50;
num_data = 20;
source_nodes = 2:49;
beta = 10^-7;
noise = 0.05;
Pe = 5;
lambda = 2.e-1; 
alpha = .001;
smoothing_coeff = .02;

hdsa = Model_Error_HDSA_React_Diff(num_data,n_mesh,beta,noise,Pe,lambda,source_nodes,smoothing_coeff);

[u_star,z_star] = hdsa.Solve_Inv_Prob();
[u_hifi,z_hifi] = hdsa.Solve_HiFi_Inv_Prob();
t = linspace(0,1,n_mesh);
z_true = exp(-10*(t-.5).^2);
z_true = z_true(source_nodes)';

z_cov = 1*ones(length(source_nodes),1);
hdsa.HDSA_Setup(u_star,z_star,alpha,z_cov);

w = ones(n_mesh,1); w(1) = .5; w(end) = .5; w = w/sum(w);
Mu = diag(w);
w = ones(length(source_nodes),1); w(1) = .5; w(end) = .5; w = w/sum(w);
Mz = diag(w);

% List of data to export
% A_hat, B_source_nodes, W_misfit, d, R, R_z, Mu, Mz
A = hdsa.A_hat;
B = hdsa.B_source_nodes;
D = hdsa.W_misfit;
data = hdsa.d;
L = hdsa.L;
Linv = inv(L);
R = hdsa.R_z;

alpha = alpha;
beta = beta;
Gamma = z_cov;

% We are solving the discrete problem
% min_{u,z} (u-d)^T*D*(u-d) + beta*z^T*R^T*R*z
% s.t. A*u = B*z
% with alpha, C, W, Gamma defined for HDSA

A = A';
A = A(:);
writematrix(A,'A.txt')
B = B';
B = B(:);
writematrix(B,'B_source.txt')
D = D';
D = D(:);
writematrix(D,'D.txt')
R = R';
R = R(:);
writematrix(R,'R.txt')
Mz = Mz';
Mz = Mz(:);
writematrix(Mz,'Mz.txt')
writematrix(u_star,'opt_state.txt')
writematrix(z_star,'opt_solution.txt')
writematrix(data,'data.txt')
L = L';
L = L(:);
writematrix(L,'L.txt')
Linv = Linv';
Linv = Linv(:);
writematrix(Linv,'Linv.txt')
writematrix(Gamma,'Gamma.txt')
writematrix(beta,'beta.txt')
writematrix(alpha,'alpha.txt')

 