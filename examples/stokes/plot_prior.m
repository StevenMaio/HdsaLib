%%
clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
N = 3*length(nodes);

%%
data_obj = importdata('control.txt', ' ', 2);  %% we need to skip the first two lines
control = data_obj.data;

% Extract x-velocity
z_prior_mean_X = control(1:3:N);
% Extract y-velocity
z_prior_mean_Y = control(2:3:N);

z_prior_samples = load('z_prior_samples.txt');
z_prior_X = z_prior_samples(1:3:N,:);
z_prior_Y = z_prior_samples(2:3:N,:);

z_prior_X_std = std(z_prior_X,[],2);
z_prior_Y_std = std(z_prior_Y,[],2);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), z_prior_mean_X);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('X-velocity control prior mean')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), z_prior_mean_Y);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Y-velocity control prior mean')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), z_prior_X_std);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('X-velocity control prior standard deviation')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), z_prior_Y_std);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Y-velocity control prior standard deviation')

%%
plot_prior_control_samples = 5;
if plot_prior_control_samples
    for k = 1:plot_prior_control_samples
        figure,
        trisurf(adj, nodes(:,1), nodes(:,2), z_prior_X(:,k));
        shading interp;
        view(0,90)
        axis equal
        axis tight
        colorbar()
        title('X-velocity control prior sample')

        figure,
        trisurf(adj, nodes(:,1), nodes(:,2), z_prior_Y(:,k));
        shading interp;
        view(0,90)
        axis equal
        axis tight
        colorbar()
        title('Y-velocity control prior sample')
    end
end

%%
discrepancy_prior = load('delta_prior_samples.txt');
discrepancy_prior_Ux = discrepancy_prior(1:3:N,:);
discrepancy_prior_Uy = discrepancy_prior(2:3:N,:);
discrepancy_prior_P = discrepancy_prior(3:3:N,:);
discrepancy_prior_Ux_std = std(discrepancy_prior_Ux,[],2);
discrepancy_prior_Uy_std = std(discrepancy_prior_Uy,[],2);
discrepancy_prior_P_std = std(discrepancy_prior_P,[],2);

figure,
trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_Ux_std);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('X-velocity discrepancy prior standard deviation')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_Uy_std);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Y-velocity discrepancy prior standard deviation')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_P_std);
shading interp;
view(0,90)
axis equal
axis tight
colorbar()
title('Pressure discrepancy prior standard deviation')

%%
plot_prior_discrepany_samples = 5;
if plot_prior_discrepany_samples
    for k = 1:plot_prior_discrepany_samples

        figure,
        trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_Ux(:,k));
        shading interp;
        view(0,90)
        axis equal
        axis tight
        colorbar()
        title('X-velocity discrepancy prior sample')

        figure,
        trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_Uy(:,k));
        shading interp;
        view(0,90)
        axis equal
        axis tight
        colorbar()
        title('Y-velocity discrepancy prior sample')

        figure,
        trisurf(adj, nodes(:,1), nodes(:,2), discrepancy_prior_P(:,k));
        shading interp;
        view(0,90)
        axis equal
        axis tight
        colorbar()
        title('Pressure discrepancy prior sample')

    end
end