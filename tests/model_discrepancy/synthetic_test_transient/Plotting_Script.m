clear
close all
clc

surpress_figures = false; %true;

c_low = 0.95;
c_high = 0.93;
n_t = 10;

diff = [];

prior_delta = cell(100,1);
prior_delta_sabl = cell(100,1);
prior_delta_sabl_load = load('Sabl_Output.mat','prior_delta').prior_delta;
for k = 1:100
    prior_delta{k} = zeros(51,3);
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_1_time_8.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,1)= data_obj.data;
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_2_time_8.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,2)= data_obj.data;
    data_obj = importdata(['prior_discrepancy_sample_',num2str(k),'/Vector_3_time_8.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta{k}(:,3)= data_obj.data;
    I = (51*7+1):(51*8);
    prior_delta_sabl{k} = prior_delta_sabl_load{k}(I,:);
end

local_diff = 0;
for k = 1:100
    local_diff = max(local_diff,norm(prior_delta{k} - prior_delta_sabl{k}));
end
diff = [diff;local_diff];

prior_delta_z_opt = zeros(51,100);
for j = 1:100
    data_obj = importdata(['prior_discrepancy_evaluated_at_z_opt/Vector_',num2str(j),'_time_6.txt'], ' ', 2);  %% we need to skip the first two lines
    prior_delta_z_opt(:,j)= data_obj.data;
end
prior_delta_z_opt_sabl = load('Sabl_Output.mat','prior_delta_z_opt').prior_delta_z_opt;
I = (51*5+1):(51*6);
prior_delta_z_opt_sabl = prior_delta_z_opt_sabl(I,:);

local_diff = norm(prior_delta_z_opt - prior_delta_z_opt_sabl);
diff = [diff;local_diff];

post_delta_mean = zeros(51,3);
for k = 1:3
    data_obj = importdata(['posterior_discrepancy_mean_',num2str(k),'_time_5.txt'], ' ', 2);  %% we need to skip the first two lines
    post_delta_mean(:,k)= data_obj.data;
end
post_delta_mean_sabl = load('Sabl_Output.mat','post_delta_mean').post_delta_mean;
I = (51*4+1):(51*5);
post_delta_mean_sabl = post_delta_mean_sabl(I,:);

local_diff = norm(post_delta_mean - post_delta_mean_sabl);
diff = [diff;local_diff];

post_delta_samples = cell(3,1);
post_delta_samples_sabl = cell(3,1);
post_delta_samples_sabl_load = load('Sabl_Output.mat','post_delta_samples').post_delta_samples;
for k = 1:3
    post_delta_samples{k} = zeros(51,100);
    for j = 1:100
        data_obj = importdata(['posterior_discrepancy_samples_',num2str(k),'/Vector_',num2str(j),'_time_5.txt'], ' ', 2);  %% we need to skip the first two lines
        post_delta_samples{k}(:,j) = data_obj.data;
    end
    I = (51*4+1):(51*5);
    post_delta_samples_sabl{k} = post_delta_samples_sabl_load{k}(I,:);
end

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
    plot(x,prior_delta_z_opt,'LineWidth',3)
    title('HdsaLib')
    figure,
    plot(x,prior_delta_z_opt_sabl,'LineWidth',3)
    title('Sabl')
    pause()
    close all

    for k = 1:5
        figure,
        plot(x,prior_delta{k},'LineWidth',3)
        title('HdsaLib')
        figure,
        plot(x,prior_delta_sabl{k},'LineWidth',3)
        title('Sabl')
        pause()
        close all
    end

    for k = 1:3
        figure
        hold on
        plot(x,post_delta_samples{k},'LineWidth',3,'Color',.9*ones(3,1))
        plot(x,post_delta_mean(:,k),'LineWidth',3,'Color','red')
        title('HdsaLib')
        figure
        hold on
        plot(x,post_delta_samples_sabl{k},'LineWidth',3,'Color',.9*ones(3,1))
        plot(x,post_delta_mean_sabl(:,k),'LineWidth',3,'Color','red')
        title('Sabl')
        pause()
        close all
    end

    figure,
    hold on
    plot(x, (c_low/c_high)^((n_t-1)/3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
    plot(x,post_z_samples,'LineWidth',3,'Color',.9*ones(3,1))
    plot(x, (c_low/c_high)^((n_t-1)/3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean,'--','LineWidth',3,'Color','red')
    title('HdsaLib')
    figure,
    hold on
    plot(x, (c_low/c_high)^((n_t-1)/3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean_sabl,'--','LineWidth',3,'Color','red')
    plot(x,post_z_samples_sabl,'LineWidth',3,'Color',.9*ones(3,1))
    plot(x, (c_low/c_high)^((n_t-1)/3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3)
    plot(x,post_z_mean_sabl,'--','LineWidth',3,'Color','red')
    title('Sabl')

end
