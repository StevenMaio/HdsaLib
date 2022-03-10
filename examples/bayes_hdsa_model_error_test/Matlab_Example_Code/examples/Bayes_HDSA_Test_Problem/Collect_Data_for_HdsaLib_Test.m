rng(12342)
Driver
close all
clc

m = length(z_star);

W = eye(m);
beta = obj.reg_beta;
D = obj.Gamma_bc_inv/sqrt(40);
data = obj.d;
u_star = obj.State_Solve(z_star);
Gamma = hdsa.Gamma;
L = hdsa.L;
Gamma_inv = linsolve(hdsa.Gamma,eye(m));
Gamma_inv_Sqrt = linsolve(sqrtm(hdsa.Gamma),eye(m));
Linv = linsolve(hdsa.L,eye(m));
Linv_Sqrt = linsolve(sqrtm(L),eye(m));
Mz = hdsa.M;

cd Data_Export/

A = L';
A = A(:);
writematrix(A,'L.txt')

A = Linv';
A = A(:);
writematrix(A,'Linv.txt')

A = Linv_Sqrt';
A = A(:);
writematrix(A,'Linv_Sqrt.txt')

A = W';
A = A(:);
writematrix(A,'W.txt')

A = D;
A = A(:);
writematrix(A,'D.txt')

A = Mz;
A = A(:);
writematrix(A,'Mz.txt')

A = Gamma;
A = A(:);
writematrix(A,'Gamma.txt')

A = Gamma_inv;
A = A(:);
writematrix(A,'Gamma_inv.txt')

A = Gamma_inv_Sqrt;
A = A(:);
writematrix(A,'Gamma_inv_Sqrt.txt')

A = Z';
A = A(:);
writematrix(A,'Z.txt')

A = Y';
A = A(:);
writematrix(A,'Y.txt')

writematrix(u_star,'opt_state.txt')
writematrix(z_star,'opt_solution.txt')
writematrix(data,'data.txt')
writematrix(beta,'beta.txt')

cd ../