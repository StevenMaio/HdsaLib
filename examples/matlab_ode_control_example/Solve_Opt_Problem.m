clear
clc
close all
addpath ../

tic;

theta = [1,2,2,1]';

finite_diff_check = false;
optimize = true;

initial_time = 0;
final_time = 1;
time_mesh_dofs = 1000;
control_dim = 5;
alpha = 10^-4;

solver = Linear_ODE(initial_time,final_time,time_mesh_dofs,control_dim,alpha);

if finite_diff_check
    solver.Finite_Difference_Test(theta);
end

if optimize
    [time_mesh,y,z] = solver.Optimize(theta);

    z_fun = zeros(length(time_mesh),1);
    for k = 1:length(z)
       z_fun = z_fun + z(k)*solver.U(:,k);
    end

    axsize = 200;
    figure('Position', [10 780 3*axsize 2*axsize]);
    figure('Position', [650 780 3*axsize 2*axsize]);
    figure('Position', [1300 780 3*axsize 2*axsize]);
    figure(1)
    plot(time_mesh,y(:,1),time_mesh,3*ones(length(time_mesh),1))
    title('y_1')
    figure(2)
    plot(time_mesh,y(:,2),time_mesh,10*time_mesh)
    title('y_2')
    figure(3)
    plot(time_mesh,z_fun)
    title('z')

end

toc