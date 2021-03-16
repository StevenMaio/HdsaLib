clear
close all
clc

load HDSA_Results.mat

L = solver.L;
Linv = linsolve(L,eye(solver.n_mesh));

W = zeros(solver.n_mesh,solver.n_mesh);
for k = 1:length(solver.I)
   W(solver.I(k),solver.I(k)) = 1; 
end

data = zeros(solver.n_mesh,1);
data(solver.I) = solver.d;

A = L';
A = A(:);
writematrix(A,'L.txt')

A = Linv';
A = A(:);
writematrix(A,'Linv.txt')

A = W';
A = A(:);
writematrix(A,'W.txt')

A = solver.D;
A = A(:);
writematrix(A,'D.txt')

A = solver.M;
A = A(:);
writematrix(A,'Mz.txt')

writematrix(solver.u_star,'opt_state.txt')
writematrix(solver.z_star,'opt_solution.txt')
writematrix(data,'data.txt')
writematrix(solver.beta,'beta.txt')
writematrix(solver.z_cov,'Gamma.txt')