clear
clc

% Set the folder where CSV files are located
folderPath = 'D:\DD119_all\kin and an\sand';

% Get a list of all CSV files in the folder
fileList = dir(fullfile(folderPath, '*.csv'));

% Initialize a cell array to store the data from each file
dataCellArray = cell(length(fileList), 1);

% Read the first file to get the original column names
firstFilePath = fullfile(folderPath, fileList(1).name);
firstFileTable = readtable(firstFilePath, 'Range', 'A1:AZ3'); % Explicitly read the first three rows
columnNames = {'frame','t1x','t1y','t2x','t2y','t3x','t3y','t4x','t4y',...
    't5x','t5y','t6x','t6y','t7x','t7y','t8x','t8y','toetipx','toetipy',...
    'footx','footy','anklex','ankley','kneex','kneey','hipx','hipy',...
    'lower_backx','lower_backy','upper_backx','upper_backy','eyex','eyey',...
    'nosex','nosey'};

% Loop through each CSV file and read the data using readtable
for i = 1:length(fileList)
    filePath = fullfile(folderPath, fileList(i).name);
    
    % Read the entire table for each file
    currentTable = readtable(filePath);
    
    % Convert the table data to an array
    dataCellArray{i} = table2array(currentTable);
end

% Check if all files have the same size
sizes = cellfun(@size, dataCellArray, 'UniformOutput', false);
if ~all(cellfun(@(x) isequal(x, sizes{1}), sizes))
    error('All CSV files must have the same number of rows and columns.');
end

% Concatenate data along rows (vertically)
mergedData = vertcat(dataCellArray{:});

% Delete every third column starting from the 4th column
mergedData(:, 4:3:end) = [];

% Assign column names to the merged data
mergedDataWithNames = array2table(mergedData, 'VariableNames', columnNames);

% Save the merged data to a new CSV file
mergedFilePath = 'D:\DD118B - invivo sand\kin and an\DD119_merged_sand.csv';
writetable(mergedDataWithNames, mergedFilePath);

disp('Merging completed successfully.');
