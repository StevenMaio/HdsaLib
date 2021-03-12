clear
clc
close all

addpath(genpath('/ascldap/users/joshart/Documents/dasco/Matlab_Trajectory_Analysis_Codes/myAutomaticDifferentiation/'))

tic;

theta = zeros(8,1);

finite_diff_check = false;
solve_uncontrolled = false;
optimize = true;

initial_time = 0;
final_time = 1;
time_mesh_dofs = 100;
control_dofs = 5;
num_controls = 2;
alpha = 0;

solver = Mirzaei_Thesis_ODE(initial_time,final_time,time_mesh_dofs,control_dofs,num_controls,alpha);

if solve_uncontrolled
    [time_mesh,y_uncontrolled] = solver.Evaluate_Uncontrolled_State(theta);
    axsize = 200;
    figure('Position', [10 780 3*axsize 2*axsize]);
    figure('Position', [650 780 3*axsize 2*axsize]);
    figure('Position', [1300 780 3*axsize 2*axsize]);
    figure('Position', [10 130 3*axsize 2*axsize]);
    figure('Position', [650 130 3*axsize 2*axsize]);
    figure('Position', [1300 130 3*axsize 2*axsize]);
    state_names = cell(6,1);
    state_names{1} = 'h';
    state_names{2} = '\theta';
    state_names{3} = '\phi';
    state_names{4} = 'v';
    state_names{5} = '\gamma';
    state_names{6} = '\psi';
    for k = 1:6
        figure(k)
        plot(time_mesh,y_uncontrolled(:,k))
        title(state_names{k})
    end
end

if finite_diff_check
    solver.Finite_Difference_Gradient_Test(theta);
    %solver.Finite_Difference_Hessian_Test(theta);
end

if optimize
    [time_mesh,y,z] = solver.Optimize(theta,100);

    alpha = zeros(length(time_mesh),1);
    for k = 1:control_dofs
       alpha = alpha + z(k)*solver.U(:,k);
    end
    sigma = zeros(length(time_mesh),1);
    for k = 1:control_dofs
       sigma = sigma + z(k+control_dofs)*solver.U(:,k);
    end
    
    [~,y_uncontrolled] = solver.Evaluate_Uncontrolled_State(theta);
    
    axsize = 200;
    figure('Position', [10 780 3*axsize 2*axsize]);
    figure('Position', [650 780 3*axsize 2*axsize]);
    figure('Position', [1300 780 3*axsize 2*axsize]);
    figure('Position', [1950 780 3*axsize 2*axsize]);
    figure('Position', [10 130 3*axsize 2*axsize]);
    figure('Position', [650 130 3*axsize 2*axsize]);
    figure('Position', [1300 130 3*axsize 2*axsize]);
    figure('Position', [1950 130 3*axsize 2*axsize]);
    state_names = cell(6,1);
    state_names{1} = 'h';
    state_names{2} = '\theta';
    state_names{3} = '\phi';
    state_names{4} = 'v';
    state_names{5} = '\gamma';
    state_names{6} = '\psi';
    for k = 1:6
        figure(k)
        plot(time_mesh,y(:,k),time_mesh,y_uncontrolled(:,k))
        title(state_names{k})
        legend('Controlled State','Uncontrolled State')
    end
    figure(7)
    plot(time_mesh,alpha)
    title('\alpha')
    figure(8)
    plot(time_mesh,sigma)
    title('\sigma')
end

toc

save('Optimal_Solution.mat')
save('control.txt','z','-ASCII')

