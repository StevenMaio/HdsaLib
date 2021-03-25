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
beta = solver.beta;
z_cov = solver.z_cov;
alpha = solver.alpha;

Sz = diag(3*zstar.^2); % Problem Specific
g = zeros(n,1);
g(I) = -(d-ustar(I));
v = zeros(n,1);
v(I) = 1;
H = diag(v);
X = kron(g',M);
A = kron(eye(n),zstar'*M);
B = Sz'*H*A + X;
v = zeros(n,1);
v(I) = 9*zstar(I).^4 - 6*zstar(I).*(d-zstar(I).^3); % Problem Specific
Hr = diag(v) + beta*D;

B = -B; % Because of code conventions

D = zeros(n,n^2);
for i = 1:n
   ei = zeros(n,1);
   ei(i) = 1;
   D(i,:) = B'*linsolve(Hr,ei); 
end

E = M*(zstar*zstar'+diag(z_cov))*M;
coeff = alpha^2*ustar'*L*ustar;
Mtheta = (1/coeff)*kron(L,E);

U = load('HDSA_Results.mat','U').U;
V = load('HDSA_Results.mat','V').V;
Sigma = load('HDSA_Results.mat','Sigma').Sigma;

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

load Singular_Values_1.txt
norm(diag(Sigmatest(1:2,1:2))'-Singular_Values_1)

load z_Singular_Vector_1.txt
U'*solver.M*z_Singular_Vector_1
norm(U-z_Singular_Vector_1,'fro')/norm(U,'fro')

e_vec = zeros(size(V,2),1);
for k = 1:size(V,2)
   e_vec(k) = min(norm(Utest(:,k)-z_Singular_Vector_1(:,k)),norm(Utest(:,k)+z_Singular_Vector_1(:,k)))/norm(Utest(:,k)); 
end
e_vec'

load theta_Singular_Vector_u1_1.txt
load theta_Singular_Vector_u2_1.txt
load theta_Singular_Vector_z1_1.txt
load theta_Singular_Vector_z2_1.txt

theta_vecs = kron(theta_Singular_Vector_u1_1,theta_Singular_Vector_z1_1) + kron(theta_Singular_Vector_u2_1,theta_Singular_Vector_z2_1);
e_vec = zeros(size(V,2),1);
for k = 1:size(V,2)
   e_vec(k) = min(norm(Vtest(:,k)-theta_vecs(:,k)),norm(Vtest(:,k)+theta_vecs(:,k)))/norm(Vtest(:,k)); 
end
e_vec'
