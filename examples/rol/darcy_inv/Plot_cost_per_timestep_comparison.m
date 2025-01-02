clear
close all
clc

reopt_cost = Compute_Cost_of_Reoptimization('perturbed_optimization.txt');

steps = 2:9;
fe_cost = 0*steps';
me_cost = 0*steps';
fe_grad = 0*steps';
me_grad = 0*steps';
cd vary_time_steps/

count = 1;
for k = steps
    load(['FE_with_',num2str(k),'_time_steps.mat']);
    fe_cost(count) = 2*(fe_results.num_B_vector_products + fe_results.num_H_vector_products + fe_results.num_gradient_evaluations);
    fe_grad(count) = fe_results.solution_gradient_norm;
    load(['ME_with_',num2str(k),'_time_steps.mat']);
    me_cost(count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
    me_grad(count) = me_results.solution_gradient_norm;
    count = count + 1;
end

cd ../

figure,
hold on
plot(steps,fe_cost,'LineWidth',3)
plot(steps,me_cost,'LineWidth',3)
plot(steps,reopt_cost*ones(length(steps),1),'--','LineWidth',3)
xlim([2,9])
%ylim([200,650])
xlabel('Number of Time Steps')
ylabel('Number of PDE Solves')
legend({'Forward Euler','Modified Euler','Reoptimization'},'Position',[0.5759    0.1869    0.3125    0.1905],'FontSize',20)
xticks(steps)
set(gca, 'fontsize', 20);