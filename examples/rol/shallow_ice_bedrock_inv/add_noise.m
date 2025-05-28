% read in problem data 
% perturb it by noise and return data.txt

% load original data from 
% original_data = load('clean_true_state.txt');
% data = load('clean_true_state.txt'); 
original_data = load('coarse_data.txt'); % use coarse data to prevent an inverse crime
[m,n] = size(original_data);
% noise level
level = .05;
%level = 0;
% fix noise for repeatable experiments
rng(111);
% create new data as pertubation of old data
data = original_data.*(1+level.*randn(m,n));


% output new data 
writematrix(data,'data.txt','Delimiter',' ')
%  fileID = fopen('data.txt', 'w');
%  fprintf(fileID,'%12.8f\n', data);
%  fclose(fileID);
 
