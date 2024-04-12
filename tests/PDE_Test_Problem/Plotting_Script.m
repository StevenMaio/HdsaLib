clear
close all
clc

delta = cell(100,1);
for k = 1:100
    delta{k} = zeros(51,3);
    delta{k}(:,1) = load(['prior_discrepancy_evaluated_at_z_',num2str(k),'/Vector_1.txt']);
    delta{k}(:,2) = load(['prior_discrepancy_evaluated_at_z_',num2str(k),'/Vector_2.txt']);
    delta{k}(:,3) = load(['prior_discrepancy_evaluated_at_z_',num2str(k),'/Vector_3.txt']);
end

delta_z_opt = zeros(51,100);
for j = 1:100
   delta_z_opt(:,j) = load(['prior_discrepancy_evaluated_at_z_opt/Vector_',num2str(j),'.txt']);
end

x = linspace(0,1,51)';

figure,
plot(x,delta_z_opt,'LineWidth',3)

for k = 1:10
   figure,
   plot(x,delta{k},'LineWidth',3)
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

for k = 1:3
    figure
    hold on
    plot(x,post_delta_samples{k},'LineWidth',3,'Color',.9*ones(3,1))
    plot(x,post_delta_mean(:,k),'LineWidth',3,'Color','red')
end

post_z_mean = load('posterior_update_mean.txt');
post_z_samples = zeros(51,100);
for k = 1:100
    post_z_samples(:,k) = load(['posterior_update_samples/Vector_',num2str(k),'.txt']);
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

load reference_solution.mat
figure,
hold on
plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3)
plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
plot(x,z_update_mean,'--','LineWidth',3,'Color','red')
plot(x,z_update_samples,'LineWidth',3,'Color',.9*ones(3,1))
plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3)
plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
plot(x,z_update_mean,'--','LineWidth',3,'Color','red')
title('Reference Solution')
