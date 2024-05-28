% Define file paths
borisFilePath = 'files for MATLAB/fixed boris/DD119/2019-08-29_09-23_Evt12-Camera1_No focal subject.csv';
kinematicsFilePath = "files for MATLAB/messy/DD119/2019-08-29_09-23_Evt12-Camera1DLC_mobnet_100_KRAT DATA WITH NEW DATAJun23shuffle1_700000.csv";
cleanedFilePath = 'C:\Users\Public\Documents\analysis\files for MATLAB\cleaned\DD119\2019-08-29_09-23_Evt12-Camera1DLC_mobnet_100_KRAT DATA CLEANED.csv';
mergedFilePath = 'C:\Users\Public\Documents\analysis\files for MATLAB/merged/DD119/DD119_event12_merged1.csv';

% Define time values from 0 to 75 seconds with 0.125-second increments
timeValues = 0:0.0625:75;

% Import your data as a table
all_data = readtable("files for MATLAB/tab/DD119 event 12.xlsx");

% Define the column names you want to keep
columnsToKeep = {'Time', 'ImageIndex', 'Behavior', 'BehaviorType'};
data = all_data(:, columnsToKeep);

% Create a list of unique behaviors
uniqueBehaviors = unique(data.Behavior);

% Define the frame indices
frameIndices = 1:1201; % Assuming you have 1201 frames

% Initialize a binary matrix
binaryMatrix = zeros(length(frameIndices), length(uniqueBehaviors));

% Iterate through the unique behaviors
for i = 1:length(uniqueBehaviors)
    behavior = uniqueBehaviors{i};

    % Convert the cell array "BehaviorType" to a string array for comparison
    behaviorTypes = string(data.BehaviorType);

    % Find the rows in the data table that correspond to the current behavior
    behaviorRows = strcmp(data.Behavior, behavior);

    % Find 'START,' 'STOP,' and 'POINT' events for the behavior
    startFrames = data.ImageIndex(behaviorRows & strcmp(behaviorTypes, 'START'));
    stopFrames = data.ImageIndex(behaviorRows & strcmp(behaviorTypes, 'STOP'));
    pointFrames = data.ImageIndex(behaviorRows & strcmp(behaviorTypes, 'POINT'));

    % Set 1s for 'START' and 'POINT' events
    binaryMatrix(ismember(frameIndices, startFrames), i) = 1;
    binaryMatrix(ismember(frameIndices, pointFrames), i) = 1;

    % Fill in between 'START' and 'STOP' with 1s
    for j = 1:length(startFrames)
        if j <= length(stopFrames)
            binaryMatrix(frameIndices >= startFrames(j) & frameIndices < stopFrames(j), i) = 1;
        end
    end
end

% Convert the binary matrix to a table with meaningful column names
binaryTable = array2table(binaryMatrix, 'VariableNames', uniqueBehaviors);

% Add the "FrameIndex" and "Time" columns to the binary table
binaryTable.FrameIndex = frameIndices(:);
binaryTable.Time = timeValues(:);

% Save the updated data to an Excel file
writetable(binaryTable, borisFilePath);

% Specify the 'VariableNamingRule' option as 'preserve'
opts = delimitedTextImportOptions('VariableNamingRule', 'preserve');

% Read the data from the CSV file into a table
data_kinematics = readtable(kinematicsFilePath, opts);

% Combine the body parts and coords into one name
row_names_kinematics = data_kinematics.Properties.RowNames;
combined_col_names_kinematics = strcat(data_kinematics{2, :}, data_kinematics{3, :});

% Remove any NaN values that might result from concatenation
combined_col_names_kinematics(ismissing(combined_col_names_kinematics)) = [];

% Convert string array to cell array of character vectors
combined_col_names_cell_kinematics = cellstr(combined_col_names_kinematics);

% Set the combined names as the new column names
data_kinematics.Properties.VariableNames = combined_col_names_cell_kinematics;

% Rename the column name from the first column
data_kinematics.Properties.VariableNames{1} = 'frame';

% Get the number of columns in the table
num_columns_kinematics = width(data_kinematics);

% Define the columns to remove (every third column starting from the fourth column)
columns_to_remove_kinematics = 4:3:num_columns_kinematics;

% Remove the 'likelihood' columns
data_kinematics(:, columns_to_remove_kinematics) = [];

% Remove first three rows (unnecessary labels)
data_kinematics(1:3, :) = [];

% Export the cleaned table to a CSV file
writetable(data_kinematics, cleanedFilePath);

% Read in kinematics and boris data files for the third part
kinematics_third = readtable(kinematicsFilePath);
boris_third = readtable(borisFilePath);

% % Removed last row from boris table because there was an extra one due to exporting as time vs frames
% rowIndex = 3233;
% boris_third(rowIndex,:) = [];

% Combine the files for the third part
merged_third = horzcat(kinematics_third, boris_third);

% Export the merged table to a CSV file
writetable(merged_third, mergedFilePath);
