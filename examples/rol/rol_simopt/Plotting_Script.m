clear
close all
clc

prior_delta = cell(100,1);
for k = 1:100
    prior_delta{k} = zeros(51,3);
    prior_delta{k}(:,1) = load(['prior_discrepancy_sample_',num2str(k),'/Vector_1.txt']);
    prior_delta{k}(:,2) = load(['prior_discrepancy_sample_',num2str(k),'/Vector_2.txt']);
    prior_delta{k}(:,3) = load(['prior_discrepancy_sample_',num2str(k),'/Vector_3.txt']);
end

prior_delta_z_opt = zeros(51,100);
for j = 1:100
    prior_delta_z_opt(:,j) = load(['prior_discrepancy_evaluated_at_z_opt/Vector_',num2str(j),'.txt']);
end

post_delta_mean = zeros(51,3);
for k = 1:3
    post_delta_mean(:,k) = load(['posterior_discrepancy_mean_',num2str(k),'.txt']);
end

post_delta_samples = cell(3,1);
for k = 1:3
    post_delta_samples{k} = zeros(51,100);
    for j = 1:100
        post_delta_samples{k}(:,j) = load(['posterior_discrepancy_samples_',num2str(k),'/Vector_',num2str(j),'.txt']);
    end
end

post_z_mean = load('posterior_update_mean.txt')';
post_z_samples = zeros(51,100);
for k = 1:100
    post_z_samples(:,k) = load(['posterior_update_samples/Vector_',num2str(k),'.txt']);
end

x = linspace(0,1,51)';

figure,
plot(x,prior_delta_z_opt,'LineWidth',3)
title('Prior discrepancy at optimal z')

for k = 1:5
    figure,
    plot(x,prior_delta{k},'LineWidth',3)
    title('Prior discrepancy at sampled z')
end

for k = 1:3
    figure,
    hold on
    plot(x,post_delta_samples{k},'LineWidth',3,'Color',.9*ones(3,1))
    plot(x,post_delta_mean(:,k),'LineWidth',3,'Color','red')
    title('Posterior discrepancy')
end

figure,
hold on
plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3)
plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
plot(x,post_z_samples,'LineWidth',3,'Color',.9*ones(3,1))
plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3)
plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
title('Optimal solution posterior')

