% script to take in clean true state and subsample on a smaller mesh

% load in fine data
fine_data = load('clean_true_state.txt');
[nt_fine,n_fine] = size(fine_data);

t_fine = linspace(0,10,nt_fine);
% define indices for s,vx,vy
s_index_fine = 1:3:n_fine;
velx_index_fine = 2:3:n_fine;
vely_index_fine = 3:3:n_fine;

s_fine = fine_data(:,s_index_fine);
velx_fine = fine_data(:,velx_index_fine);
vely_fine = fine_data(:,vely_index_fine);

% define fine mesh
h_fine = sqrt(n_fine/3);
[X,Y] = meshgrid(linspace(0,1,h_fine)); % this is hard coded (need to improve)

nt = 61; % number of time steps for coarse problem
h = 71;
n = 3*h^2;

% define indices for s,vx,vy (coarse)
s_index_coarse = 1:3:n;
velx_index_coarse = 2:3:n;
vely_index_coarse = 3:3:n;

s_coarse = zeros(nt_fine,n/3); % matrices for coarse data
velx_coarse = zeros(nt_fine,n/3);
vely_coarse = zeros(nt_fine,n/3);

[Xq,Yq] = meshgrid(linspace(0,1,h)); % query points 

% interpolating  mesh at each time step
for i = 1:nt_fine
    Vs = reshape(s_fine(i,:),h_fine,h_fine);
    s = interp2(X,Y,Vs,Xq,Yq);
    s_coarse(i,:) = reshape(s,1,n/3);
    Vvelx = reshape(velx_fine(i,:),h_fine,h_fine);
    vx = interp2(X,Y,Vvelx,Xq,Yq);
    velx_coarse(i,:) = reshape(vx,1,n/3);
    Vvely = reshape(vely_fine(i,:),h_fine,h_fine);
    vy = interp2(X,Y,Vvely,Xq,Yq);
    vely_coarse(i,:) = reshape(vy,1,n/3);
end

t_coarse = linspace(0,10,nt);

% interpolating  time over mesh
s = interp1(t_fine(:),s_coarse,t_coarse);
vx = interp1(t_fine(:),velx_coarse,t_coarse);
vy = interp1(t_fine(:),vely_coarse,t_coarse);

% putting data back
coarse_data = zeros(nt,n);
for i = 1:nt
    coarse_data(i,s_index_coarse) = s(i,:);
    coarse_data(i,velx_index_coarse) = vx(i,:);
    coarse_data(i,vely_index_coarse) = vy(i,:);
end

% write to a file named coarse_data.txt
writematrix(coarse_data,'coarse_data.txt','Delimiter',' ')


