function [results] = Read_Results(filename)

% Open the file
fileID = fopen(filename, 'r');

results = struct();
% Initialize variables
results.num_B_vector_products = 0;
results.num_H_vector_products = 0;
results.num_gradient_evaluations = 0;
results.num_vectors_stored = 0;
results.solution_gradient_norm = 0;

% Read the file line by line
while ~feof(fileID)
    line = fgetl(fileID); % Read a line from the file

    % Check for each line and extract the corresponding value
    if contains(line, 'Number of B-vector products:')
        results.num_B_vector_products = str2double(extractAfter(line, ': '));
    elseif contains(line, 'Number of H-vector products:')
        results.num_H_vector_products = str2double(extractAfter(line, ': '));
    elseif contains(line, 'Number of gradient evaluations:')
        results.num_gradient_evaluations = str2double(extractAfter(line, ': '));
    elseif contains(line, 'Number of vectors stored for preconditioner:')
        results.num_vectors_stored = str2double(extractAfter(line, ': '));
    elseif contains(line, 'Solution gradient norm =')
        results.solution_gradient_norm = str2double(extractAfter(line, '= '));
    end
end

% Close the file
fclose(fileID);
end