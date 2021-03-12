clear
close all
clc

plotresults

z_sing = load('z_Singular_Vector_1.txt');
theta_sing = load('theta_Singular_Vector_1.txt');
sing_vals = load('Singular_Values_1.txt');

figure,
trisurf(adj, nodes(:,1), nodes(:,2), sing_vals(1)*z_sing(:,1));
shading interp;
colorbar();
view(0,90)
axis square
title('First source singular vector')

figure,
trisurf(adj, nodes(:,1), nodes(:,2), sing_vals(2)*z_sing(:,2));
shading interp;
colorbar();
view(0,90)
axis square
title('Second source singular vector')