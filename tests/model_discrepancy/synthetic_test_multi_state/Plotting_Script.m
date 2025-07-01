clear
close all
clc

surpress_figures = false; %true;

diff = [];

prior_delta = cell(100,1);
for k = 1:100
    prior_delta{k} = zeros(102,3);
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_1.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,1)= data_obj.data;
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_2.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,2)= data_obj.data;
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_3.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,3)= data_obj.data;
end
prior_delta_sabl = load('Sabl_Output.mat','prior_delta').prior_delta;

local_diff = 0;
for k = 1:100
    local_diff = max(local_diff,norm(prior_delta{k} - prior_delta_sabl{k}));
end
diff = [diff;local_diff];

prior_delta_z_opt = zeros(102,100);
for j = 1:100
    data_obj = importdata(['prior_discrepancy_evaluated_at_z_opt/Vector_',num2str(j),'.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta_z_opt(:,j)= data_obj.data;
end
prior_delta_z_opt_sabl = load('Sabl_Output.mat','prior_delta_z_opt').prior_delta_z_opt;

local_diff = norm(prior_delta_z_opt - prior_delta_z_opt_sabl);
diff = [diff;local_diff];

post_delta_mean = zeros(102,3);
for k = 1:3
    data_obj = importdata(['posterior_discrepancy_mean_',num2str(k),'.txt'], ' ', 2);  %% we need to skip the first two lines
    post_delta_mean(:,k)= data_obj.data;
end
post_delta_mean_sabl = load('Sabl_Output.mat','post_delta_mean').post_delta_mean;

local_diff = norm(post_delta_mean - post_delta_mean_sabl);
diff = [diff;local_diff];

post_delta_samples = cell(3,1);
for k = 1:3
    post_delta_samples{k} = zeros(102,100);
    for j = 1:100
        data_obj = importdata(['posterior_discrepancy_samples_',num2str(k),'/Vector_',num2str(j),'.txt'], ' ', 2);  %% we need to skip the first two lines
        post_delta_samples{k}(:,j) = data_obj.data;
    end
end
post_delta_samples_sabl = load('Sabl_Output.mat','post_delta_samples').post_delta_samples;

local_diff = 0;
for k = 1:3
    local_diff = max(local_diff,norm(post_delta_samples{k} - post_delta_samples_sabl{k}));
end
diff = [diff;local_diff];

data_obj = importdata('posterior_update_mean.txt', ' ', 2);  %% we need to skip the first two lines
post_z_mean = data_obj.data;
post_z_samples = zeros(51,100);
for k = 1:100
    data_obj = importdata(['posterior_update_samples/Vector_',num2str(k),'.txt'], ' ', 2);  %% we need to skip the first two lines
    post_z_samples(:,k) = data_obj.data;
end
post_z_mean_sabl = load('Sabl_Output.mat','post_z_mean').post_z_mean;
post_z_samples_sabl = load('Sabl_Output.mat','post_z_samples').post_z_samples;

local_diff = norm(post_z_mean - post_z_mean_sabl);
diff = [diff;local_diff];

local_diff = norm(post_z_samples - post_z_samples_sabl);
diff = [diff;local_diff];


if ~surpress_figures

    x = linspace(0,1,51)';

    figure,
    plot(x,prior_delta_z_opt(1:51,:),'LineWidth',3)
    title('HdsaLib')
    figure,
    plot(x,prior_delta_z_opt_sabl(1:51,:),'LineWidth',3)
    title('Sabl')
    pause()
    close all

    for k = 1:5
        figure,
        plot(x,prior_delta{k}(1:51,:),'LineWidth',3)
        title('HdsaLib')
        figure,
        plot(x,prior_delta_sabl{k}(1:51,:),'LineWidth',3)
        title('Sabl')
        pause()
        close all
    end

    for k = 1:3
        figure
        hold on
        plot(x,post_delta_samples{k}(1:51,:),'LineWidth',3,'Color',.9*ones(3,1))
        plot(x,post_delta_mean(1:51,k),'LineWidth',3,'Color','red')
        title('HdsaLib')
        figure
        hold on
        plot(x,post_delta_samples_sabl{k}(1:51,:),'LineWidth',3,'Color',.9*ones(3,1))
        plot(x,post_delta_mean_sabl(1:51,k),'LineWidth',3,'Color','red')
        title('Sabl')
        pause()
        close all
    end

    figure,
    hold on
    plot(x, (1 + x) / (1.1^(2 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
    plot(x,post_z_samples,'LineWidth',3,'Color',.9*ones(3,1))
    plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3)
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
    title('HdsaLib')
    figure,
    hold on
    plot(x, (1 + x) / (1.1^(2 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean_sabl,'--','LineWidth',3,'Color','red')
    plot(x,post_z_samples_sabl,'LineWidth',3,'Color',.9*ones(3,1))
    plot(x, (1 + x) / (1.1^(2 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean_sabl,'--','LineWidth',3,'Color','red')
    title('Sabl')

end
