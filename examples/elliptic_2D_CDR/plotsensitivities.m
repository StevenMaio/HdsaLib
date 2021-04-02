clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing

nodes = load('nodes.txt');  %% load node coordinates

z_sing_vecs = load('z_Singular_Vector_1.txt');
sing_vals = load('Singular_Values_1.txt');

theta_sing_vecs_u1 = load('theta_Singular_Vector_u1_1.txt');
theta_sing_vecs_u2 = load('theta_Singular_Vector_u2_1.txt');
theta_sing_vecs_z1 = load('theta_Singular_Vector_z1_1.txt');
theta_sing_vecs_z2 = load('theta_Singular_Vector_z2_1.txt');
M_theta_sing_vecs_z1 = load('Mz_theta_Singular_Vector_z1_1.txt');
M_theta_sing_vecs_z2 = load('Mz_theta_Singular_Vector_z2_1.txt');

data_obj = importdata('state.txt', ' ', 2);  %% we need to skip the first two lines
state = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), state);
shading interp;
colorbar();
view(0,90)
axis square
title('State')

data_obj = importdata('source.txt', ' ', 2);  %% we need to skip the first two lines
source = data_obj.data;
figure,
trisurf(adj, nodes(:,1), nodes(:,2), source);
shading interp;
colorbar();
view(0,90)
axis square
title('Source')

noisy_state = load('data.txt');
figure,
trisurf(adj, nodes(:,1), nodes(:,2), noisy_state);
shading interp;
colorbar();
view(0,90)
axis square
title('Noisy State')

z_pert = sing_vals(1)*z_sing_vecs(:,1);
z_eval = source + 0*z_pert;
delta1_zstar = theta_sing_vecs_u1(:,1)*(z_eval'*M_theta_sing_vecs_z1) + theta_sing_vecs_u2*(z_eval'*M_theta_sing_vecs_z2(:,1));
figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta1_zstar);
shading interp;
colorbar();
view(0,90)
axis square
title('1st Model Perturbation')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), z_pert);
shading interp;
colorbar();
view(0,90)
axis square
title('1st Source Perturbation')

z_pert = sing_vals(2)*z_sing_vecs(:,2);
z_eval = source + z_pert;
delta2_zstar = theta_sing_vecs_u1(:,2)*(z_eval'*M_theta_sing_vecs_z1) + theta_sing_vecs_u2*(z_eval'*M_theta_sing_vecs_z2(:,2));
figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta2_zstar);
shading interp;
colorbar();
view(0,90)
axis square
title('2nd Model Perturbation')

figure,
trisurf(adj, nodes(:,1), nodes(:,2),  z_pert);
shading interp;
colorbar();
view(0,90)
axis square
title('2nd Source Perturbation')
