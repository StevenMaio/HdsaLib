clear
close all
clc

write_fine_mesh = false;

%load gis_2d.mat
load gis_2d_basal.mat

% basal_sliding = 1./(nvar01*1000); % In m/(Pa*yr)
% bed_topo = nvar02*1000; % In m
% thickness = nvar04*1000; % In m
% forcing = nvar12; % In m/yr
% xvel = nvar06; % In m/yr
% yvel = nvar07; % In m/yr
% x = x0; % In km
% y = y0; % In km
basal_sliding = 1./(nvar03*1000); % In m/(Pa*yr)
bed_topo = nvar04*1000; % In m
thickness = nvar05*1000; % In m
forcing = nvar01; % In m/yr
xvel = nvar06; % In m/yr
yvel = nvar07; % In m/yr
x = x0; % In km
y = y0; % In km

surface_height = bed_topo + thickness; % Same as nvar11 (surface_height in nnames)

xmin = -300;
xmax = 250;
ymin = -2340;
ymax = -1890;
xb = 50;
yb = 50;
I = intersect(intersect(find(x>xmin-xb),find(x<xmax+xb)),intersect(find(y>ymin-yb),find(y<ymax+yb)));

m = 31; % mesh size for extrapolation
n = 21; % mesh size for interpolation after smoothing
bandwidth = 1; % number of nearest neighbors used for averaging 


X = linspace(xmin,xmax,m);
Y = linspace(ymin,ymax,m);

[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));
[X_grid2,Y_grid2] = meshgrid(linspace(xmin,xmax,n),linspace(ymin,ymax,n));

figure,
trisurf(blk01',x,y,surface_height)
hold on
M = 1.1*max(xvel);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
xlabel('x (km)')
ylabel('y (km)')
title('Surface Height (m)')
ax = gca;
set(ax,'xticklabel',[])
set(ax,'yticklabel',[])

x_data = x(I);
y_data = y(I);
z_data = surface_height(I);
model = fit([x_data,y_data],z_data,'linearinterp');
K = find(x_data<-200);
mdl = fitlm([x_data(K),y_data(K)],z_data(K));
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid_nans = Z_grid;
Z_grid_linear = predict(mdl,[X_grid(:),Y_grid(:)]);
J = isnan(Z_grid);
Z_grid(J) = Z_grid_linear(J);
Z_grid = reshape(Z_grid,m,m);

figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Surface Height (m)')
xlabel('x (km)')
ylabel('y (km)')
ax = gca;
set(ax,'xticklabel',[])
set(ax,'yticklabel',[])
axis([xmin xmax, ymin ymax])


Z_grid = Smoothing_Func(Z_grid,bandwidth);

if write_fine_mesh
    local_surface_height = Z_grid;
end

% we interpolate back onto a nxn mesh after the smoothing
Z_grid = interp2(X_grid,Y_grid,Z_grid,X_grid2,Y_grid2);

figure,
surf(X_grid2,Y_grid2,Z_grid)
view(2)
shading interp
colorbar
title('Localized Surface Height with smoothing (m)')
xlabel('x (km)')
ylabel('y (km)')
axis([xmin xmax, ymin ymax])

if ~write_fine_mesh
    local_surface_height = Z_grid;
end


figure,
trisurf(blk01',x,y,thickness)
hold on
M = 1.1*max(xvel);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Thickness')

x_data = x(I);
y_data = y(I);
z_data = thickness(I);
model = fit([x_data,y_data],z_data,'linearinterp');
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid(isnan(Z_grid))=0;
Z_grid = reshape(Z_grid,m,m);

figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Thickness without smoothing')
axis([xmin xmax, ymin ymax])

Z_grid = Smoothing_Func(Z_grid,bandwidth);

if write_fine_mesh
    local_thickness = Z_grid;
end

% interpolate back to nxn mesh
Z_grid = interp2(X_grid,Y_grid,Z_grid,X_grid2,Y_grid2);

figure,
surf(X_grid2,Y_grid2,Z_grid)
view(2)
shading interp
colorbar
title('Localized Thickness with smoothing')
axis([xmin xmax, ymin ymax])

if ~write_fine_mesh
    local_thickness = Z_grid;
end

figure,
trisurf(blk01',x,y,bed_topo)
hold on
M = 1.1*max(xvel);
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


Z_grid = local_surface_height - local_thickness;

figure,
if ~write_fine_mesh
    surf(X_grid2,Y_grid2,Z_grid)
else
    surf(X_grid,Y_grid,Z_grid)
end
view(2)
shading interp
colorbar
title('Bedrock Toopography (m)')
xlabel('x (km)')
ylabel('y (km)')
axis([xmin xmax, ymin ymax])
local_bedrock = Z_grid;

figure,
trisurf(blk01',x,y,log(basal_sliding))
hold on
M = 1.1*max(xvel);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Log Basal Sliding (m)')


x_data = x;
y_data = y;
z_data = log(basal_sliding);
model = fit([x_data,y_data],z_data,'linearinterp');
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);

% smoothing the sliding coefficient
Z_grid = Smoothing_Func(Z_grid,bandwidth);

if write_fine_mesh
    local_basal_sliding = Z_grid;
end

% we interpolate back onto a nxn mesh after the smoothing
Z_grid = interp2(X_grid,Y_grid,Z_grid,X_grid2,Y_grid2);

figure,
surf(X_grid2,Y_grid2,Z_grid)
view(2)
shading interp
colorbar
title('Log Basal Sliding (m)')
xlabel('x (km)')
ylabel('y (km)')
ax = gca;
set(ax,'xticklabel',[])
set(ax,'yticklabel',[])
axis([xmin xmax, ymin ymax])


if ~write_fine_mesh
    local_basal_sliding = Z_grid;
end

figure,
trisurf(blk01',x,y,forcing)
hold on
M = 1.1*max(xvel);
for i = 1:m
    scatter3(X(i),ymax,M,'o','black') 
    scatter3(X(i),ymin,M,'o','black') 
    scatter3(xmin,Y(i),M,'o','black') 
    scatter3(xmax,Y(i),M,'o','black') 
end
view(2)
shading interp
colorbar
title('Forcing')


x_data = x(I);
y_data = y(I);
z_data = forcing(I);
model = fit([x_data,y_data],z_data,'linearinterp');
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
I = isnan(Z_grid);
Z_grid(I) = 0;
Z_grid = reshape(Z_grid,m,m);


% Creating 0 forcing for testing purposes
%Z_grid = 0*Z_grid;
% 

% smoothing the forcing
Z_grid = Smoothing_Func(Z_grid,bandwidth);

if write_fine_mesh
    local_forcing = Z_grid;
end

% interpolate on a nxn mesh after smoothing
Z_grid = interp2(X_grid,Y_grid,Z_grid,X_grid2,Y_grid2);

figure,
surf(X_grid2,Y_grid2,Z_grid)
view(2)
shading interp
colorbar
title('Accumulation/Ablation (m/yr)')
xlabel('x (km)')
ylabel('y (km)')
ax = gca;
set(ax,'xticklabel',[])
set(ax,'yticklabel',[])
axis([xmin xmax, ymin ymax])

if ~write_fine_mesh
    local_forcing = Z_grid;
end

Z_grid = local_basal_sliding';
fileID = fopen('Log_Basal_Sliding.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

Z_grid = local_bedrock';
fileID = fopen('Bedrock_Topography.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

Z_grid = local_surface_height';
fileID = fopen('Surface_Height.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);

Z_grid = local_forcing';
fileID = fopen('Forcing.txt', 'w');
fprintf(fileID,'%12.8f\n', Z_grid(:));
fclose(fileID);
