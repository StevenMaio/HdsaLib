clear
close all
clc

reopt_cost = Compute_Cost_of_Reoptimization('perturbed_optimization.txt');

steps = 1:8;
fe_cost = 0*steps';
me_cost = 0*steps';
fe_grad = 0*steps';
me_grad = 0*steps';

count = 1;
for k = steps
    fe_results = Read_Results(['Forward_Euler_Cost_Report_',num2str(k),'.txt']);
    fe_cost(count) = 2*(fe_results.num_B_vector_products + fe_results.num_H_vector_products + fe_results.num_gradient_evaluations);
    fe_grad(count) = fe_results.solution_gradient_norm;
    me_results = Read_Results(['Modified_Euler_Cost_Report_',num2str(k),'.txt']);
    me_cost(count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
    me_grad(count) = me_results.solution_gradient_norm;
    count = count + 1;
end


figure,
semilogx(10.^(-steps),fe_cost,'LineWidth',3)
hold on
semilogx(10.^(-steps),me_cost,'LineWidth',3)
semilogx(10.^(-steps),reopt_cost*ones(length(steps),1),'--','LineWidth',3)
xlabel('CG Tolerance')
ylabel('Number of PDE Solves')
legend({'Forward Euler','Modified Euler','Reoptimization'},'Position',[0.3705    0.5464    0.3125    0.1905],'FontSize',20)
set(gca, 'fontsize', 20);
saveas(gcf,'vary_tolerance_cost_comparison','epsc')