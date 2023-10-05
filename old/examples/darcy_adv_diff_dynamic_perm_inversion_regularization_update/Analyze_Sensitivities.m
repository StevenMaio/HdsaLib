clear
close all
clc
FS = 18;

evals = load('LIS_evals_1.txt');
evecs = load('LIS_evecs_1.txt');
W = load('eigenvector_weightmat_products_1.txt');
B = load('B.txt');

m = size(evecs,2);
n = size(B,2);
S = zeros(n,m);

for k = 1:m
   for j = 1:n
      u = evecs(:,1:k)'*B(:,j)./(1+evals(1:k)');
      S(j,k) = sqrt(u'*W(1:k,1:k)*u); 
   end
end

I = find(evals>1);
k = I(end);

L = 20;
b = 1;
e = 9*16;
source_indices = S(b:e,k);
b = e + 1;
e = b + 1 - 1;
diff_index = S(b:e,k);
b = e + 1;
e = b + (L+1) - 1;
left_bdry_indices = S(b:e,k);
b = e + 1;
e = b + (L+1) - 1;
right_bdry_indices = S(b:e,k);
uncertain_basis_grid = linspace(0,1,L+1);

jet = colormap('jet');
jet(end-5:end,:) = []; jet(1:5,:) = [];
close all

xl = {0.2, 0.2, 0.2, 0.2, 0.4, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.6, 0.8, 0.8, 0.8, 0.8};
yl = {0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8, 0.2, 0.4, 0.6, 0.8};
figure,
hold on
for k = 1:16
    count = 1;
    for j = 1:3
        y_pt = yl{k} + (j-2)*.2/4;
       for i = 1:3
            x_pt = xl{k} + (i-2)*.2/4;
            scatter(x_pt,y_pt,[],source_indices((k-1)*9+count),'filled');
            count = count + 1;
       end
    end
end
title('Tracer Injection Sensitivities','fontsize',FS)
colormap(jet)
colorbar
xticks([0,.25,.5,.75,1])
yticks([0,.25,.5,.75,1])
xlim([0,1])
ylim([0,1])
caxis([0,max(source_indices)])
set(gca, 'FontSize', 16);
 
figure,
plot(uncertain_basis_grid,left_bdry_indices,'ob','MarkerFaceColor','b','Linewidth',2);
xlabel('y'); ylabel('Sensitivity');
title('Left Boundary Sensitivities');
set(gca, 'FontSize', 16);

figure,
plot(uncertain_basis_grid,right_bdry_indices,'ob','MarkerFaceColor','b','Linewidth',2);
xlabel('y'); ylabel('Sensitivity');
title('Right Boundary Sensitivities');
set(gca, 'FontSize', 16);

figure,
hold on
for k = 1:n
   plot(S(k,:)) 
end
I = find(evals<1);
xline(I(1),'-',['\lambda = ',num2str(evals(I(1)))]);
I = find(evals<.5);
xline(I(1),'-',['\lambda = ',num2str(evals(I(1)))]);
I = find(evals<.1);
xline(I(1),'-',['\lambda = ',num2str(evals(I(1)))]);
xlabel('Number of eigenvectors'); ylabel('Sensitivity Index')
title('Sensitivity Indices for Different Likelihood Informed Subspace Dimensions')
set(gca, 'FontSize', 16);
