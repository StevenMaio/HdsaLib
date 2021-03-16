
clear
close all
clc

solver = Model_Error_HDSA_Syn_Test();
[ustar,zstar] = solver.Solve_Inv_Prob();
x = solver.x;

alpha = .08;
z_cov = 1*ones(length(x),1);
solver.HDSA_Setup(ustar,zstar,alpha,z_cov);

k = 2;
p = 14;
q = 1;
[U,Sigma,V] = solver.Compute_HDSA_GSVD(k,p,q);

axsize = 200;
figure('Position', [100 280 3*axsize 2*axsize]);
figure('Position', [750 280 3*axsize 2*axsize]);
figure('Position', [1400 280 3*axsize 2*axsize]);
figure('Position', [100 900 3*axsize 2*axsize]);
figure('Position', [750 900 3*axsize 2*axsize]);
figure('Position', [1400 900 3*axsize 2*axsize]);

figure(1)
hold on
plot(x,zstar)
plot(x,solver.ztrue)
xlabel('x')
ylabel('z')
legend({'Estimated Parameter','True Parameter'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

figure(4)
hold on
plot(x,ustar,x,solver.utrue)
scatter(x(solver.I),solver.d,'o')
xlabel('x')
ylabel('u')
legend({'Estimated State','True State'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

figure(2)
hold on
plot(x,zstar,x,zstar+Sigma(1,1)*U(:,1))
title('First Singular Vector')
xlabel('x')
ylabel('z')
legend({'Nominal','Perturbed'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

figure(3)
hold on
plot(x,zstar,x,zstar+Sigma(2,2)*U(:,2))
title('Second Singular Vector')
xlabel('x')
ylabel('z')
legend({'Nominal','Perturbed'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

Theta = reshape(V(:,1),length(x),length(x))';
d1 = Theta*solver.Apply_z_Mass(zstar);
figure(5)
hold on
plot(x,ustar,x,ustar+d1)
title('First Singular Vector Perturbation')
xlabel('x')
ylabel('u')
legend({'Nominal','Perturbed'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

Theta = reshape(V(:,2),length(x),length(x))';
d2 = Theta*solver.Apply_z_Mass(zstar);
figure(6)
hold on
plot(x,ustar,x,ustar+d2)
title('Second Singular Vector Perturbation')
xlabel('x')
ylabel('u')
legend({'Nominal','Perturbed'})
set(gca, 'FontSize', 16); set(gcf, 'Color', 'White');

save('HDSA_Results.mat');