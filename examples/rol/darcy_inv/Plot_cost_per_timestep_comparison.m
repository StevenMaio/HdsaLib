clear
close all
clc

reopt_cost = Compute_Cost_of_Reoptimization('perturbed_optimization.txt');

steps = 2:9;
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

max([fe_grad;me_grad;pre_fe_grad;pre_me_grad])

colors = lines(2);

figure,
hold on
plot(steps,pre_fe_cost,'LineWidth',3,'Color',colors(1,:))
plot(steps,pre_me_cost,'LineWidth',3,'Color',colors(2,:))
plot(steps,fe_cost,'--','LineWidth',3,'Color',colors(1,:))
plot(steps,me_cost,'--','LineWidth',3,'Color',colors(2,:))
xlim([2,9])
%ylim([200,600])
xlabel('Number of Time Steps')
ylabel('Number of PDE Solves')
legend({'Preconditioned Forward Euler','Preconditioned Modified Euler','Forward Euler','Modified Euler'},'Position',[0.5759    0.1869    0.3125    0.1905],'FontSize',20)
xticks(steps)
set(gca, 'fontsize', 20);
%saveas(gcf,'vary_time_steps_1_cost_comparison','epsc')