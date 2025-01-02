function [params] = Read_Inputs_from_Xml()
% Specify the filename
filename = 'Sensitivity_input.xml'; % Replace with your actual file name

% Read the XML file
xmlDoc = xmlread(filename);

% Get the 'ParameterList' node that contains the 'Problem' parameters
parameterList = xmlDoc.getElementsByTagName('ParameterList');
problemParameters = parameterList.item(0).getElementsByTagName('Parameter');

% Initialize variables
params = struct();

% Loop through each 'Parameter' node and extract the attributes
for i = 0:problemParameters.getLength()-1
    paramNode = problemParameters.item(i);

    % Get the name and value attributes
    paramName = char(paramNode.getAttribute('name'));
    paramName = strrep(paramName, ' ', '_');
    paramValue = char(paramNode.getAttribute('value'));
    paramType = char(paramNode.getAttribute('type'));

    % Convert the value to the appropriate type
    switch paramType
        case 'double'
            params.(paramName) = str2double(paramValue);
        case 'int'
            params.(paramName) = str2double(paramValue);
        case 'bool'
            params.(paramName) = strcmp(paramValue, 'true');
        otherwise
            params.(paramName) = paramValue; % Fallback for other types
    end
end

end