clear
clc

% Define time values from 0 to 75 seconds with 0.125-second increments
timeValues = 0:0.0625:75;

% Import your data as a table
all_data = readtable("files for MATLAB/tab/DD119 event 5 fixed.xlsx");

% Define the column names you want to keep
columnsToKeep = {'Time', 'ImageIndex', 'Behavior', 'BehaviorType'};
data = all_data(:, columnsToKeep);

% Define the frame indices
frameIndices = 1:1201; % Assuming you have 1201 frames

% Define a list of all possible behaviors
% allBehaviors = {'HandObstruction', 'Landing', 'NewVideo', 'Quadrupedal', 'Stance', 'Steady Hopping', 'Swing', 'Takeoff', 'Wedge_WedgeBack_', 'Wedge_WedgeFront_', 'Wedge'};
allBehaviors = {'HandObstruction', 'Landing', 'New Video', 'Quadrupedal', 'Stance', 'SteadyHopping', 'Swing', 'Takeoff'};

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

% Define the desired file path and name
file_path = 'D:\DD119_all\BORIS\2019-08-29_09-00_Evt05-Camera1_No focal subject.csv';

% Save the updated data to an Excel file
writetable(binaryTable, file_path);

