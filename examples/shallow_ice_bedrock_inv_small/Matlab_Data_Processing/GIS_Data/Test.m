clear
clc
close all

c = load('Log_Basal_Friction.txt');

m = 81;
[X,Y] = meshgrid(linspace(0,1,m));
Z = zeros(m,m);

for i = 1:m
   for j = 1:m
       Z(i,j) = FE_Eval(X(i,j),Y(i,j),c,sqrt(length(c))-1);
   end
end

figure,
surf(X,Y,Z)
view(2)
shading interp
colorbar
title('FE Localized Log Basal Friction')

open Log_Basal_Friction_Figure.fig

v = linspace(0,1,round(1.5*m));
v = v(2:(end-1));
m2 = length(v);
[X2,Y2] = meshgrid(v);
Z_x_diff = zeros(m2,m2);
for i = 1:m2
   for j = 1:m2
       Z_x_diff(i,j) = FE_Eval_x_diff(X2(i,j),Y2(i,j),c,sqrt(length(c))-1);
   end
end
figure,
surf(X2,Y2,Z_x_diff)
view(2)
shading interp
colorbar
title('FE Localized Log Basal Friction x Derivative')
h = 10^-3;
Z_x_pert = zeros(m,m);
for i = 1:m
   for j = 1:m
       Z_x_pert(i,j) = FE_Eval(X(i,j)+h,Y(i,j),c,sqrt(length(c))-1);
   end
end
Z_x_pert(:,end) = Z(:,end);
figure,
surf(X,Y,(Z_x_pert-Z)/h)
view(2)
shading interp
colorbar
title('FE Localized Log Basal Friction x FD')


Z_y_diff = zeros(m2,m2);
for i = 1:m2
   for j = 1:m2
       Z_y_diff(i,j) = FE_Eval_y_diff(X2(i,j),Y2(i,j),c,sqrt(length(c))-1);
   end
end
figure,
surf(X2,Y2,Z_y_diff)
view(2)
shading interp
colorbar
title('FE Localized Log Basal Friction y Derivative')
h = 10^-3;
Z_y_pert = zeros(m,m);
for i = 1:m
   for j = 1:m
       Z_y_pert(i,j) = FE_Eval(X(i,j),Y(i,j)+h,c,sqrt(length(c))-1);
   end
end
Z_y_pert(end,:) = Z(end,:);
figure,
surf(X,Y,(Z_y_pert-Z)/h)
view(2)
shading interp
colorbar
title('FE Localized Log Basal Friction y FD')

