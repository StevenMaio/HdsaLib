clear
clc
close all
rng(1234)

n_mesh = 50;
num_data = 20;
source_nodes = 2:49;
beta_reg = 10^-7;
noise = 0.05;
Pe = 5;
lambda = 2.e-1; 
smoothing_coeff = .02;

hdsa = Model_Error_HDSA_React_Diff(num_data,n_mesh,beta_reg,noise,Pe,lambda,source_nodes,smoothing_coeff);

[u_star,z_star] = hdsa.Solve_Inv_Prob();
[u_hifi,z_hifi] = hdsa.Solve_HiFi_Inv_Prob();
t = linspace(0,1,n_mesh);
z_true = exp(-10*(t-.5).^2);
z_true = z_true(source_nodes)';

z_cov = 1*ones(length(source_nodes),1);
hdsa.HDSA_Setup(u_star,z_star,z_cov);

S = 6;
p = 10;
q = 1;
[U,Sigma,V] = hdsa.Compute_HDSA_GSVD(S,p,q);

load Singular_Values_1.txt
max(abs(diag(Sigma(1:S,1:S))-Singular_Values_1'))

load z_Singular_Vector_1.txt
for k = 1:S
   min(norm(U(:,k)-z_Singular_Vector_1(:,k))/norm(U(:,k)),norm(U(:,k)+z_Singular_Vector_1(:,k))/norm(U(:,k)))
end

load theta_Singular_Vector_a_1.txt
abs(V.a-theta_Singular_Vector_a_1)/norm(V.a)

load theta_Singular_Vector_bk_1.txt
load theta_Singular_Vector_uk_1.txt
load theta_Singular_Vector_zk_1.txt
for k = 1:S
    min(abs(theta_Singular_Vector_bk_1(k)-V.b_k(k)),abs(theta_Singular_Vector_bk_1(k)+V.b_k(k)))
    min(norm(theta_Singular_Vector_uk_1(:,k)-V.u_k(:,k))/norm(V.u_k(:,k)),norm(theta_Singular_Vector_uk_1(:,k)+V.u_k(:,k))/norm(V.u_k(:,k)))
    min(norm(theta_Singular_Vector_zk_1(:,k)-V.z_k(:,k))/norm(V.z_k(:,k)),norm(theta_Singular_Vector_zk_1(:,k)+V.z_k(:,k))/norm(V.z_k(:,k)))
end
load theta_Singular_Vector_u_1.txt
load theta_Singular_Vector_z_1.txt
min(norm(theta_Singular_Vector_u_1-V.u)/norm(V.u),norm(theta_Singular_Vector_u_1+V.u)/norm(V.u))
min(norm(theta_Singular_Vector_z_1-V.z)/norm(V.z),norm(theta_Singular_Vector_z_1+V.z)/norm(V.z))

