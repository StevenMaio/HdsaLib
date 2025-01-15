clear
close all
clc

mesh = [50 , 100 , 200 , 400];

me_cost = 0*mesh';
me_storage = 0*mesh';
me_grad = 0*mesh';

count = 1;
for k = mesh
    cd(['vary_mesh_',num2str(k)]);
    me_results = Read_Results('Modified_Euler_Cost_Report.txt');
    me_cost(count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
    me_grad(count) = me_results.solution_gradient_norm;
    me_storage(count) = me_results.num_vectors_stored;
    count = count + 1;
    cd('../')
end

figure,

yyaxis left; 
plot(mesh,me_cost, '-o', 'DisplayName', 'Dataset 1','LineWidth',3);
ylabel('Number of PDE Solves'); 
ylim([300,500])
yyaxis right; 
plot(mesh,me_storage, '-s', 'DisplayName', 'Dataset 2', 'Color', 'r','LineWidth',3); 
ylabel('Number of Vectors Stored'); 
ylim([0,250])
xlim([49,401])
xlabel('Mesh Resolution')
set(gca, 'fontsize', 20);
saveas(gcf,'vary_mesh_comparison','epsc')