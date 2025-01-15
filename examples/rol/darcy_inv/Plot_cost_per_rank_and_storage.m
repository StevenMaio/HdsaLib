clear
close all
clc

storage = 0:5:20;
rank = 0:5:20;

me_cost = zeros(length(storage),length(rank));
me_storage = zeros(length(storage),length(rank));
me_grad = zeros(length(storage),length(rank));

storage_count = 1;
for k = storage
    rank_count = 1;
    for j = rank
        me_results = Read_Results(['Modified_Euler_Cost_Report_',num2str(j),'_',num2str(k),'.txt']);
        me_cost(storage_count,rank_count) = 2*(me_results.num_B_vector_products + me_results.num_H_vector_products + me_results.num_gradient_evaluations);
        me_grad(storage_count,rank_count) = me_results.solution_gradient_norm;
        me_storage(storage_count,rank_count) = me_results.num_vectors_stored;
        rank_count = rank_count + 1;
    end
    storage_count = storage_count + 1;
end

figure,
hold on
for k = 1:length(rank)
    plot(storage,me_cost(:,k),'LineWidth',3)
end
xlabel('Maximum Update Rank')
ylabel('Number of PDE Solves')
legend({'Initial Rank 0','Initial Rank 5','Initial Rank 10','Initial Rank 15','Initial Rank 20'},'FontSize',20)
set(gca, 'fontsize', 20);
saveas(gcf,'vary_rank_cost_comparison','epsc')

figure,
hold on
for k = 1:length(rank)
    plot(storage,me_storage(:,k),'LineWidth',3)
end
xlabel('Maximum Update Rank')
ylabel('Number of Vectors Stored')
legend({'Initial Rank 0','Initial Rank 5','Initial Rank 10','Initial Rank 15','Initial Rank 20'},'FontSize',20,'Position',[0.5607    0.2893    0.2982    0.3726])
set(gca, 'fontsize', 20);
saveas(gcf,'vary_rank_storage_comparison','epsc')