clear
close all
clc

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

m = 51;
%m = 71;
X = linspace(xmin,xmax,m);
Y = linspace(ymin,ymax,m);

[X_grid,Y_grid] = meshgrid(linspace(xmin,xmax,m),linspace(ymin,ymax,m));

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
title('Surface Height')

x_data = x(I);
y_data = y(I);
z_data = surface_height(I);
model = fit([x_data,y_data],z_data,'linearinterp');
K = find(x_data<-100);
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
title('Localized Surface Height')
local_surface_height = Z_grid;

% figure,
% Z_grid_nans = reshape(Z_grid_nans,m,m);
% surf(X_grid,Y_grid,Z_grid_nans)
% view(2)
% shading interp
% colorbar
% title('Localized Surface Height')

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

%% Averaging Code
% this is where we can smooth the thickness
s = Z_grid;
w = [1 0 1 0 1 0 1 0]/8 ; % averaging vector
for i = 1:m
    j =1;
    if i == 1
      stencil_values = [ s(i+1,j), s(i+1,j+1), s(i,j+1),0, ...
        0, 0 , 0 , 0];  
    s(i,j) = w(:)'*stencil_values(:);
    elseif i == m
        stencil_values = [ 0 , 0 , s(i,j+1),s(i-1,j+1), ...
        s(i-1,j), 0 , 0 , 0];
        s(i,j) = w(:)'*stencil_values(:);
    else 
    stencil_values = [ s(i+1,j), s(i+1,j+1), s(i,j+1),s(i-1,j+1), ...
        s(i-1,j), 0 , 0 , 0];
        s(i,j) = w(:)'*stencil_values(:);
    end
    for j = 2:m-1
        if i == 1
        stencil_values = [s(i+1,j), s(i+1,j+1),s(i,j+1), ...
            0, 0, 0, s(i,j-1), 0];    
        s(i,j) = w(:)'*stencil_values(:);
        elseif i == m
        stencil_values = [0, 0 ,s(i,j+1), ...
            s(i-1,j+1), s(i-1,j), s(i-1,j-1), s(i,j-1), s(i-1,j-1)];
        s(i,j) = w(:)'*stencil_values(:);
        else
        stencil_values = [s(i+1,j), s(i+1,j+1),s(i,j+1), ...
            s(i-1,j+1), s(i-1,j), s(i-1,j-1), s(i,j-1), s(i-1,j-1)];
        s(i,j) = w(:)'*stencil_values(:);
        end
    end
    j = m;
    if i == 1 
        stencil_values = [s(i+1,j), 0 , 0, 0, 0, 0, s(i,j-1)...
        , s(i+1,j-1)];
    s(i,j) = w(:)'*stencil_values(:);
    elseif i == m
        stencil_values = [0, 0 , 0, 0, s(i-1,j), s(i-1,j-1), s(i,j-1)...
        , 0];
    s(i,j) = w(:)'*stencil_values(:);
    else
    stencil_values = [s(i+1,j), 0 , 0, 0, s(i-1,j), s(i-1,j-1), s(i,j-1)...
        , s(i+1,j-1)];
    s(i,j) = w(:)'*stencil_values(:);
    end
end
figure,
surf(X_grid,Y_grid,s)
view(2)
shading interp
colorbar
title('Localized Thickness')

Z_grid = s;
local_thickness = Z_grid;
%%
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
title('Bedrock')

Z_grid = local_surface_height - local_thickness;
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Bedrock')

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
title('Log Basal Sliding')

x_data = x;
y_data = y;
z_data = log(basal_sliding);
model = fit([x_data,y_data],z_data,'linearinterp');
Z_grid = feval(model,[X_grid(:),Y_grid(:)]);
Z_grid = reshape(Z_grid,m,m);
figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Log Basal Sliding')

local_basal_sliding = Z_grid;

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
Z_grid = 0*Z_grid;
% 

figure,
surf(X_grid,Y_grid,Z_grid)
view(2)
shading interp
colorbar
title('Localized Forcing')

local_forcing = Z_grid;

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