clear
close all
clc

adj = load('cell_to_node_quad.txt') + 1;  %% load node adjacency table, increment by 1 for 1-based indexing
nodes = load('nodes.txt');  %% load node coordinates
zbar = load('control_read.txt');
obj = Process_Perturbations(zbar);
K = length(obj.sing_vals);

c = zeros(K,1); c(1) = 1;
[z_pert,delta_I,delta_L_u1,delta_L_z1,delta_L_u2,delta_L_z2] = obj.Evaluate_Perturbation(c);

N = 4*size(nodes,1);
xvel = 1:4:N;
yvel = 2:4:N;
pres = 3:4:N;
thr = 4:4:N;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_I(xvel));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for X-velocity model error at nominal')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_I(yvel));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for Y-velocity model error at nominal')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_I(pres));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for pressure model error at nominal')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), delta_I(thr));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for thermal model error at nominal')
axis equal
axis tight

d = delta_I + delta_L_u1*delta_L_z1'*z_pert + delta_L_u2*delta_L_z2'*z_pert;

figure,
trisurf(adj, nodes(:,1), nodes(:,2), d(xvel));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for X-velocity model error at perturbed')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), d(yvel));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for Y-velocity model error at perturbed')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), d(pres));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for pressure model error at perturbed')
axis equal
axis tight

figure,
trisurf(adj, nodes(:,1), nodes(:,2), d(thr));
shading interp;
view(0,90)
colorbar()
title('Leading singular vector for thermal model error at perturbed')
axis equal
axis tight

zp = z_pert(thr);
figure,
idx = find(nodes(:,1) == 0.0);
plot(nodes(idx,2), zp(idx));
title('Leading singular vector for left control')

figure,
idx = find(nodes(:,1) == 1.0);
plot(nodes(idx,2), zp(idx));
title('Leading singular vector for right control')
