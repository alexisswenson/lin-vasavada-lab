function binary_clean_merge(boris, kinematicsFilePath, cleanedFilePath, mergedFilePath)
   
% Specify the 'VariableNamingRule' option as 'preserve'
opts = delimitedTextImportOptions('VariableNamingRule', 'preserve');

% Define time values from 0 to 75 seconds with 0.125-second increments
    timeValues = 0:0.0625:75;

    % Script 1: Generate Binary Table   
% Import your data as a table
all_data = readtable("files for MATLAB/tab/DD119 event 12.xlsx");

% Define the column names you want to keep
columnsToKeep = {'Time', 'ImageIndex', 'Behavior', 'BehaviorType'};
data = all_data(:, columnsToKeep);

% Define the frame indices
frameIndices = 1:1201; % Assuming you have 1201 frames

% Define a list of all possible behaviors
allBehaviors = {'HandObstruction', 'Landing', 'NewVideo', 'Quadrupedal', 'Stance', 'SteadyHopping', 'Swing', 'Takeoff', 'Wedge_WedgeBack_', 'Wedge_WedgeFront_', 'Wedge'};

% Initialize a binary matrix with columns for all possible behaviors
binaryMatrix = zeros(length(frameIndices), length(allBehaviors));

% Iterate through all possible behaviors
for i = 1:length(allBehaviors)
    behavior = allBehaviors{i};

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
binaryTable = array2table(binaryMatrix, 'VariableNames', allBehaviors);

% Add the "FrameIndex" and "Time" columns to the binary table
binaryTable.FrameIndex = frameIndices(:);
binaryTable.Time = timeValues(:);

% Save the updated data to an Excel file
writetable(binaryTable, boris);


  
    % Script 2: Clean Kinematics Data
    % Specify the 'VariableNamingRule' option as 'preserve'
% opts = delimitedTextImportOptions('VariableNamingRule', 'preserve');

% Read the data from the CSV file into a table
file_path = "files for MATLAB/messy/DD119/2019-08-29_09-23_Evt12-Camera1DLC_mobnet_100_KRAT DATA WITH NEW DATAJun23shuffle1_700000.csv";
data = readtable(file_path, opts);

% combine the bodyparts and coords into one name
row_names = data.Properties.RowNames;
combined_col_names = strcat(data{2, :}, data{3, :});

% Remove any NaN values that might result from concatenation
combined_col_names(ismissing(combined_col_names)) = [];

% Set the combined names as the new column names
data.Properties.VariableNames = combined_col_names;

% Rename the column name from the first column
data.Properties.VariableNames{1} = 'frame';

% Get the number of columns in the table
num_columns = width(data);

% Define the columns to remove (every third column starting from the fourth column)
columns_to_remove = 5:3:num_columns;

% Remove the 'likelihood' columns
data(:, columns_to_remove) = [];

%remove first three rows (unnecesarry labels)
data(1:3, :) = [];

% Save the cleaned table to a CSV file
writetable(data, cleanedFilePath);

    % Script 3: Merge Kinematics and Boris Data
    kinematics_third = readtable(kinematicsFilePath);
    boris_third = readtable(boris);

    merged_third = horzcat(kinematics_third, boris_third);
    
    % Save the merged table to a CSV file
    writetable(merged_third, mergedFilePath);
end