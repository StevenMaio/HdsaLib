clear
close all
clc

surpress_figures = false; %true;

n_t = 20;

diff = [];

post_z_mean = zeros(2,n_t);
for j = 1:n_t
    post_z_mean(:,j) = load(['posterior_update_mean_time_',num2str(j),'.txt']);
end
post_z_mean = post_z_mean(:);

post_z_samples = zeros(2*n_t,100);
for k = 1:100
    post_z_samples_tmp = zeros(2,n_t);
    for j = 1:n_t
        post_z_samples_tmp(:,j) = load(['posterior_update_samples/Vector_',num2str(k),'_time_',num2str(j),'.txt']);
    end
    post_z_samples(:,k) = post_z_samples_tmp(:);
end
post_z_mean_sabl = load('Sabl_Output.mat','post_z_mean').post_z_mean;
post_z_samples_sabl = load('Sabl_Output.mat','post_z_samples').post_z_samples;

local_diff = norm(post_z_mean - post_z_mean_sabl)/norm(post_z_samples_sabl);
diff = [diff;local_diff];

local_diff = norm(post_z_samples - post_z_samples_sabl,'fro')/norm(post_z_samples_sabl,'fro');
diff = [diff;local_diff];


if ~surpress_figures

    figure,
    plot(1:40,post_z_mean_sabl,1:40,post_z_mean(:),'--','LineWidth',3)

    figure,
    plot(1:40,post_z_samples,'LineWidth',3,'Color',[.9,.9,.9])
    title('HdsaLib')

    figure,
    plot(1:40,post_z_samples_sabl,'LineWidth',3,'Color',[.9,.9,.9])
    title('Sabl')

end
