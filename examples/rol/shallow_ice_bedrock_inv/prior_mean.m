% create the prior mean

m0 = zeros(441,1);
fileID = fopen('matlab_prior_mean.txt', 'w');
fprintf(fileID,'%12.8f\n', m0);
fclose(fileID);