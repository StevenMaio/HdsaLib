clear
close all
clc

reopt_cost = Compute_Cost_of_Reoptimization('perturbed_optimization.txt');

steps = 1:8;
fe_cost = 0*steps';
me_cost = 0*steps';
fe_grad = 0*steps';
me_grad = 0*steps';
pre_fe_cost = 0*steps';
pre_me_cost = 0*steps';
pre_fe_grad = 0*steps';
pre_me_grad = 0*steps';

count = 1;
for k = steps
    fe_results = Read_Results(['Forward_Euler_Cost_Report_',num2str(k),'.txt']);
    fe_cost(count) = 2*(fe_results.num_B_vector_products + fe_results.num_H_vector_products + fe_results.num_gradient_evaluations);
    fe_grad(count) = fe_results.solution_gradient_norm;
    me_results = Read_Results(['Modified_Euler_Cost_Report_',num2str(k),'.txt']);
    me_cost(count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
    me_grad(count) = me_results.solution_gradient_norm;
    fe_results = Read_Results(['Preconditioned_Forward_Euler_Cost_Report_',num2str(k),'.txt']);
    pre_fe_cost(count) = 2*(fe_results.num_B_vector_products + fe_results.num_H_vector_products + fe_results.num_gradient_evaluations);
    pre_fe_grad(count) = fe_results.solution_gradient_norm;
    me_results = Read_Results(['Preconditioned_Modified_Euler_Cost_Report_',num2str(k),'.txt']);
    pre_me_cost(count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
    pre_me_grad(count) = me_results.solution_gradient_norm;
    count = count + 1;
end

colors = lines(2);

figure,
semilogx(10.^(-steps),pre_fe_cost,'LineWidth',3,'Color',colors(1,:))
hold on
semilogx(10.^(-steps),pre_me_cost,'LineWidth',3,'Color',colors(2,:))
semilogx(10.^(-steps),fe_cost,'--','LineWidth',3,'Color',colors(1,:))
semilogx(10.^(-steps),me_cost,'--','LineWidth',3,'Color',colors(2,:))
xlabel('CG Tolerance')
ylabel('Number of PDE Solves')
xlim(10.^[-8,-2])
legend({'Preconditioned Forward Euler','Preconditioned Modified Euler','Forward Euler','Modified Euler'},'Position',[0.5759    0.1869    0.3125    0.1905],'FontSize',20)
set(gca, 'fontsize', 20);
%saveas(gcf,'vary_tolerance_cost_comparison','epsc')