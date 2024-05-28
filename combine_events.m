clc
clear

% Set the folder where CSV files are located
folderPath = 'D:\DD119_all\kinematics\1.75 m per s\sand\400 hz';

% Get a list of all CSV files in the folder
fileList = dir(fullfile(folderPath, '*.csv'));

% Initialize a cell array to store the data from each file
dataCellArray = cell(length(fileList), 1);

% Loop through each CSV file and read the data using readtable
for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
    currentTable = readtable(filePath);
    
    % Skip the first row and convert the table data to an array
    dataCellArray{i} = table2array(currentTable(2:end, :));
end

% Check if all files have the same size
sizes = cellfun(@size, dataCellArray, 'UniformOutput', false);
if ~all(cellfun(@(x) isequal(x, sizes{1}), sizes))
    error('All CSV files must have the same number of rows and columns.');
end

% Concatenate data along columns (assuming files have the same number of rows)
mergedData = horzcat(dataCellArray{:});

% Delete every third column starting from the 4th column
mergedData(:, 4:3:end) = [];

% Define the desired column names
desiredColumnNames = {'t1x', 't1y', 't2x', 't2y', 't3x', 't3y', 't4x', 't4y', ...
    't5x', 't5y', 't6x', 't6y', 't7x', 't7y', 't8x', 't8y', 'toetipx', 'toetipy', ...
    'footx', 'footy', 'anklex', 'ankley', 'kneex', 'kneey', 'hipx', 'hipy', ...
    'lower_backx', 'lower_backy', 'upper_backx', 'upper_backy', 'eyex', 'eyey', ...
    'nosex', 'nosey', 'time'};

% Assign column names to the merged data
mergedDataWithNames = array2table(mergedData, 'VariableNames', desiredColumnNames);

% Save the merged data to a new CSV file
mergedFilePath = 'D:\DD119_all\merged\merged_practice.csv';
writetable(mergedDataWithNames, mergedFilePath);

disp('Merging completed successfully.');
