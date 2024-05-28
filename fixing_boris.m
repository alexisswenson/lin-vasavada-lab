clear
clc

% Import your data as a table
data = readtable("files for MATLAB/badboris/2019-10-15_15-09_Evt07-Camera1_No focal subject_No focal subject.xlsx");

% Define the desired final time (150s)
finalTime = 150;

% Calculate the number of rows to reach 150s
desiredRows = finalTime / 0.125 + 1;

% Create a table to store the extended data
extendedData = table('Size', [desiredRows, size(data, 2)], 'VariableTypes', varfun(@class, data, 'OutputFormat', 'cell'));

% Copy the column names from the data table to the extended data table
extendedData.Properties.VariableNames = data.Properties.VariableNames;

% Initialize the index for the extended data
currentIndex = 1;

% Copy the data and create new rows with the same data as the previous row
while currentIndex <= desiredRows
    extendedData(currentIndex, :) = data(ceil(currentIndex / 2), :);
    currentIndex = currentIndex + 1;
end

% Define the desired file path and name
file_path = 'D:\DD118B - invivo sand\BORIS\Solid\2019-10-15_15-09_Evt07-Camera1_No focal subject.csv';

% Save the updated data to an Excel file
writetable(extendedData, file_path);
