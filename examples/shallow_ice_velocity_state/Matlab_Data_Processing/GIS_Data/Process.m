clear
close all
clc

load gis_2d.mat

basal_sliding = 1./(nvar01*1000); % In m/(Pa*yr)
bed_topo = nvar02*1000; % In m
thickness = nvar04*1000; % In m
forcing = nvar12; % In m/yr
x = x0; % In km
y = y0; % In km

surface_height = bed_topo + thickness; % Same as nvar11 (surface_height in nnames)

xmin = 225;
xmax = 325;
ymin = -2100;
ymax = -2000;
xb = 50;
yb = 50;
I = intersect(intersect(find(x>xmin-xb),find(x<xmax+xb)),intersect(find(y>ymin-yb),find(y<ymax+yb)));

m = 31;
X = linspace(xmin,xmax,m);
Y = linspace(ymin,ymax,m);

figure,
trisurf(blk01',x,y,log(basal_sliding))
hold on
M = 1.1*max(log(basal_sliding));
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Log Basal Sliding')

x_data = x(I);
y_data = y(I);
z_data = log(basal_sliding(I));
model = fit([x_data,y_data],z_data,'linearinterp');
[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Log Basal Sliding')

Z_grid = Z_grid';
fileID = fopen('Log_Basal_Sliding.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

figure,
trisurf(blk01',x,y,bed_topo)
hold on
M = 1.1*max(bed_topo);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Bedrock Topography')

x_data = x(I);
y_data = y(I);
z_data = bed_topo(I);
model = fit([x_data,y_data],z_data,'linearinterp');
[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Bedrock Topography')

Z_grid = Z_grid';
fileID = fopen('Bedrock_Topography.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

figure,
trisurf(blk01',x,y,thickness)
hold on
M = 1.1*max(thickness);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Ice Thickness')

x_data = x(I);
y_data = y(I);
z_data = thickness(I);
model = fit([x_data,y_data],z_data,'linearinterp');
[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Thickness')

% Z_grid = Z_grid';
% fileID = fopen('Thickness.txt', 'w');
% fprintf(fileID,'%12.8f\n', Z_grid(:));
% fclose(fileID);

figure,
trisurf(blk01',x,y,surface_height)
hold on
M = 1.1*max(surface_height);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Surface Height')

x_data = x(I);
y_data = y(I);
z_data = surface_height(I);
model = fit([x_data,y_data],z_data,'linearinterp');
[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Surface Height')

Z_grid = Z_grid';
fileID = fopen('Surface_Height.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

figure,
trisurf(blk01',x,y,forcing)
hold on
M = 1.1*max(forcing);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Accumulation/Ablation Forcing')

x_data = x(I);
y_data = y(I);
z_data = forcing(I);
model = fit([x_data,y_data],z_data,'linearinterp');
[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Accumulation/Ablation Forcing')

Z_grid = Z_grid';
fileID = fopen('Forcing.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);
