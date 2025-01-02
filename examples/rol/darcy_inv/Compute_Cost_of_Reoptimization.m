function [total_pde_solves] = Compute_Cost_of_Reoptimization(filename)

fileID = fopen(filename, 'r');

% Step 2: Read the data
% Skip the first two lines (header and title)
headerLines = 3;
dataFormat = '%f %f %f %f %f %f %f %f %f %f'; % Define the format for the data
data = textscan(fileID, dataFormat, 'HeaderLines', headerLines, 'CollectOutput', true);

% Step 3: Close the file
fclose(fileID);

% Step 4: Process the data
% The data is stored in a cell array, we can convert it to a matrix
dataMatrix = data{1}; % Extract the numerical data

total_pde_solves = 2*nansum(dataMatrix(:,9)) + dataMatrix(end,6) + dataMatrix(end,7);
end