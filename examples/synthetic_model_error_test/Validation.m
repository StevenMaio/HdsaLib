%%
clear
close all
clc

solver = load('HDSA_Results.mat','solver').solver;
ustar = load('HDSA_Results.mat','ustar').ustar;
zstar = load('HDSA_Results.mat','zstar').zstar;

n = solver.n_mesh;
d = solver.d;
I = solver.I;
M = solver.M;
L = solver.L;
D = solver.D;
reg_beta = solver.reg_beta;
Gamma = solver.Gamma;

Sz = diag(3*zstar.^2); % Problem Specific
g = zeros(n,1);
g(I) = -(d-ustar(I));
v = zeros(n,1);
v(I) = 1;
H = diag(v);
X = kron(g',M);
X = [zeros(n,n),X];
A = kron(eye(n),zstar'*M);
A = [eye(n),A];
B = Sz'*H*A + X;
v = zeros(n,1);
v(I) = 9*zstar(I).^4 - 6*zstar(I).*(d-zstar(I).^3); % Problem Specific
Hr = diag(v) + reg_beta*D;

B = -B; % Because of code conventions

D = zeros(n,n*(n+1));
for i = 1:n
   ei = zeros(n,1);
   ei(i) = 1;
   D(i,:) = B'*linsolve(Hr,ei); 
end

E = M*(zstar*zstar'+Gamma)*M;
Mtheta = [L,kron(L,zstar'*M);kron(L,M*zstar),kron(L,E)]; 

U = load('HDSA_Results.mat','U').U;
V = load('HDSA_Results.mat','V').V;
Sigma = load('HDSA_Results.mat','Sigma').Sigma;
k = load('HDSA_Results.mat','k').k;

U = U(:,1:k);
V = V(:,1:k);
Sigma = Sigma(1:k,1:k);

norm(D*V-U*Sigma,'fro')/norm(D*V,'fro')
norm(V'*Mtheta*V-eye(size(V,2)))
norm(U'*M*U-eye(size(U,2)))

Ltheta = chol(Mtheta);
Ltheta = Ltheta'; % Mtheta = Ltheta*Ltheta'
Lz = chol(M);
Lz = Lz'; 
Test = linsolve(Ltheta,D'*Lz);
Test = Test';
[Utest,Sigmatest,Vtest] = svd(Test);
Utest = linsolve(Lz',Utest);
Vtest = linsolve(Ltheta',Vtest);

norm(Sigma-Sigmatest(1:size(V,2),1:size(V,2)))
e_vec = zeros(size(V,2),1);
for k = 1:size(V,2)
   e_vec(k) = min(norm(V(:,k)-Vtest(:,k)),norm(V(:,k)+Vtest(:,k)))/norm(V(:,k)); 
end
e_vec'
e_vec = zeros(size(V,2),1);
for k = 1:size(V,2)
   e_vec(k) = min(norm(U(:,k)-Utest(:,k)),norm(U(:,k)+Utest(:,k)))/norm(U(:,k)); 
end
e_vec'

%%
load Singular_Values_1.txt
norm(diag(Sigmatest(1:2,1:2))'-Singular_Values_1)

load z_Singular_Vector_1.txt
U'*solver.M*z_Singular_Vector_1
S = size(U,2);
for k = 1:S
    min(norm(U(:,k)-z_Singular_Vector_1(:,k))/norm(U(:,k)),norm(U(:,k)+z_Singular_Vector_1(:,k))/norm(U(:,k)))
end

load theta_Singular_Vector_zk_1.txt
load theta_Singular_Vector_z_1.txt
load theta_Singular_Vector_uk_1.txt
load theta_Singular_Vector_u_1.txt
load theta_Singular_Vector_bk_1.txt
load theta_Singular_Vector_a_1.txt

theta_Singular_Vector_1 = zeros(size(V));
for k = 1:size(V,2)
   theta_Singular_Vector_1(:,k) = [ theta_Singular_Vector_a_1*theta_Singular_Vector_uk_1(:,k) ; kron(theta_Singular_Vector_uk_1(:,k),theta_Singular_Vector_z_1) ] ...
       + [ theta_Singular_Vector_bk_1(k)*theta_Singular_Vector_u_1 ; kron(theta_Singular_Vector_u_1,theta_Singular_Vector_zk_1(:,k)) ];
end

V'*Mtheta*theta_Singular_Vector_1
for k = 1:S
    min(norm(V(:,k)-theta_Singular_Vector_1(:,k))/norm(V(:,k)),norm(V(:,k)+theta_Singular_Vector_1(:,k))/norm(V(:,k)))
end

load theta_Singular_Vector_Mz_1.txt
load theta_Singular_Vector_Mzk_1.txt

norm(theta_Singular_Vector_Mz_1-M*theta_Singular_Vector_z_1)
norm(theta_Singular_Vector_Mzk_1-M*theta_Singular_Vector_zk_1)

% val_oper = struct();
% val_oper.B = B;
% val_oper.Mz = M;
% val_oper.Mtheta = Mtheta;
% val_oper.Hr = Hr;
% val_oper.Sz = Sz;
% val_oper.H = H;
% save('val_oper.mat','val_oper')
