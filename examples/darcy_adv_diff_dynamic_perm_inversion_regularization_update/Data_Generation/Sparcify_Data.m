close all; clear; clc;

nt = 76;                    %% number of time steps
ni = 128;                   %% initial FE mesh dimension (double desired dimension)
nf = (ni/2)+1;              %% desired final dimension of data
m = nf^2;

nodes = load('nodes.txt');  %% load node coordinates
clean_data = load('clean_true_state.txt'); %% load clean data
clean_data = clean_data(1:2:end,:);

noisy_data_new = zeros(nt,nf*nf*2);

nodesx = reshape(nodes(:,1),ni+1,ni+1);
nodesy = reshape(nodes(:,2),ni+1,ni+1);
nodesx(:,2:2:end) = [];
nodesx(2:2:end,:) = [];
nodesy(:,2:2:end) = [];
nodesy(2:2:end,:) = [];

pressure_noise = 1 + .01*randn(m,1);

for k = 1:nt
    datap = clean_data(2,1:2:end);
    datap = reshape(datap,ni+1,ni+1);
    datap(:,2:2:end) = [];
    datap(2:2:end,:) = [];
    datap = reshape(datap,size(datap,1)*size(datap,2),1).*pressure_noise;
    
    datac = clean_data(k,2:2:end);
    datac = reshape(datac,ni+1,ni+1);
    datac(:,2:2:end) = [];
    datac(2:2:end,:) = [];
    datac = reshape(datac,size(datac,1)*size(datac,2),1).*(1 + .01*randn(m,1));
    
    noisy_data_new(k,1:2:end) = datap;
    noisy_data_new(k,2:2:end) = datac;
end 

writematrix(noisy_data_new,'data.txt','Delimiter','tab');



