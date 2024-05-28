clc
clear

% Set paths to the folders
emg_kinematics_path = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\merged_data';
boris_path = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\BORIS';

% Get a list of files from EMG kinematics folder
emgFileList = dir(fullfile(emg_kinematics_path, '*.csv'));

% Initialize cell array to store EMG data matrices
emgDataCellArray = cell(length(emgFileList), 1);

% Read data from each EMG file and store in the cell array
for i = 1:length(emgFileList)
    filePath = fullfile(emg_kinematics_path, emgFileList(i).name);
    data = readtable(filePath); % Assuming numeric data in CSV files
    emgDataCellArray{i} = data;
end

% Merge EMG data tables vertically
mergedEMGData = vertcat(emgDataCellArray{:});

% Identify columns to check for NaN values (starting from the fifth column)
nanColumns = mergedEMGData.Properties.VariableNames(5:end);

% Initialize logical mask to identify rows with NaN values in selected columns
nanMask = false(height(mergedEMGData), 1);

% Loop through selected columns and update nanMask
for i = 1:numel(nanColumns)
    columnData = mergedEMGData.(nanColumns{i});
    nanMask = nanMask | (isnan(columnData) & ~strcmp(nanColumns{i}, 'SonoLG')); % Exclude 'SonoLG' column
end

% Keep rows where nanMask is false
cleanedEMGData = mergedEMGData(~nanMask, :);



% Get a list of files from BORIS folder
borisFileList = dir(fullfile(boris_path, '*.csv'));

% Initialize cell array to store BORIS data matrices
borisDataCellArray = cell(length(borisFileList), 1);

% Read data from each BORIS file and store in the cell array
for i = 1:length(borisFileList)
    filePath = fullfile(boris_path, borisFileList(i).name);
    data = readtable(filePath); % Assuming numeric data in CSV files
    borisDataCellArray{i} = data;
end

% Merge BORIS data matrices vertically
mergedBorisData = vertcat(borisDataCellArray{:});

% Merge cleaned EMG data and BORIS data horizontally
mergedData = [cleanedEMGData, mergedBorisData];

% Define the old and new column names mapping
oldColumnNames = {'Hand_Obstruction', 'New_Video', 'Steady_Hopping'};
newColumnNames = {'HandObstruction', 'NewVideo', 'SteadyHopping'};

% Check if the columns exist in the mergedData table
existingColumns = intersect(oldColumnNames, mergedData.Properties.VariableNames);

% Rename the columns if they exist
if ~isempty(existingColumns)
    for i = 1:numel(existingColumns)
        mergedData.Properties.VariableNames{existingColumns{i}} = newColumnNames{i};
    end
end

% % Save merged data to a new CSV file
% mergedFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_merged_data.csv'; % Update with your desired file path
% writetable(mergedData, mergedFilePath);


% %DD119 solid (first 3 files (rows 1-3,603) and last 9 files (rows 22,820-33,628))
% %DD119 sand (files 3-18 (rows 3,604-22,819))
% 
% % Row indices for solid and sand sections
% solidRows = [1:3603, 22820:33628];
% sandRows = 3604:22819;
% 
% % Create solid and sand data tables
% solidData = mergedData(solidRows, :);
% sandData = mergedData(sandRows, :);
% 
% % Save solid and sand data to separate CSV files
% solidFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_solid_merged_data.csv';  % Specify the file path for solid data
% sandFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_sand_merged_data.csv';    % Specify the file path for sand data
% 
% writetable(solidData, solidFilePath);
% writetable(sandData, sandFilePath);


%DD121 solid (first 11 files (rows 1-13,211))
%DD121 sand (last 8 files (rows 13,212-22819))

% % Row indices for solid and sand sections
% solidRows = 1:13211;
% sandRows = 13212:22819;
% 
% % Create solid and sand data tables
% solidData = mergedData(solidRows, :);
% sandData = mergedData(sandRows, :);
% 
% % Save solid and sand data to separate CSV files
% solidFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD121\DD121_solid_merged_data.csv';  % Specify the file path for solid data
% sandFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD121\DD121_sand_merged_data.csv';    % Specify the file path for sand data
% 
% writetable(solidData, solidFilePath);
% writetable(sandData, sandFilePath);


% DD118A solid (first 17 files (rows 1-20,417))
% DD118A sand (next 10 files (rows 20,418-32,427))

% % Row indices for solid and sand sections
% solidRows = 1:20417;
% sandRows = 20418:32427;
% 
% % Create solid and sand data tables
% solidData = mergedData(solidRows, :);
% sandData = mergedData(sandRows, :);
% 
% % Save solid and sand data to separate CSV files
% solidFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118A\DD118A_solid_merged_data.csv';  % Specify the file path for solid data
% sandFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118A\DD118A_sand_merged_data.csv';    % Specify the file path for sand data
% 
% writetable(solidData, solidFilePath);
% writetable(sandData, sandFilePath);
% 
% disp('Solid and sand data files created successfully.');


% DD118B solid (first 14 files (rows 1-16,814)) I excluded evt 9
% DD118B sand (next 22 files (rows 16,815-43,236))

% Row indices for solid and sand sections
solidRows = 1:16814;
sandRows = 16815:43236;

% Create solid and sand data tables
solidData = mergedData(solidRows, :);
sandData = mergedData(sandRows, :);

% Save solid and sand data to separate CSV files
solidFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\DD118B_solid_merged_data.csv';  % Specify the file path for solid data
sandFilePath = 'C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\DD118B_sand_merged_data.csv';    % Specify the file path for sand data

writetable(solidData, solidFilePath);
writetable(sandData, sandFilePath);

disp('Solid and sand data files created successfully.');