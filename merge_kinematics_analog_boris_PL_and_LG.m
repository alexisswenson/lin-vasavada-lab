clear
clc

%%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters (All you need to change)
animal_number = 'SCDD15';
sand_or_solid = 'solid';
file_name = 'SCDD15_solid';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the folder where CSV files are located
folderPath = sprintf('D:\\%s\\kin and an\\%s', animal_number, sand_or_solid);

% Get a list of all CSV files in the folder
fileList = dir(fullfile(folderPath, '*.csv'));

% Initialize a cell array to store cleaned data
cleanedDataCell = cell(numel(fileList), 1);

% Loop through each file in the folder
for i = 1:numel(fileList)
    % Read the current CSV file
    filePath = fullfile(folderPath, fileList(i).name);
    data = readtable(filePath);
    
    % Step 2: Perform filtering and data processing
    Combined_Tendon = data.CombinedTendon;
    LG_Sono = data.LGSono;
    PL_Sono = data.PLSono;
    LG_EMG = data.LGEMG;
    PL_EMG = data.PLEMG;

    % Remove mean from EMG_1, EMG_2, Buckle
    LG_EMG = detrend(LG_EMG);
    PL_EMG = detrend(PL_EMG);

    % Rectify absolute values
    LG_EMG = abs(LG_EMG);
    PL_EMG = abs(PL_EMG);

    % Check for NaN values in SONO
    nanIndices = isnan(LG_Sono);
    % Check for infinite values in SONO
    infIndices = isinf(LG_Sono);
    % Replace NaN and infinite values with zeros
    LG_Sono(nanIndices) = 0;
    LG_Sono(infIndices) = 0;

    % Check for NaN values in EMG_1
    nanIndices = isnan(LG_EMG);
    % Check for infinite values in EMG_1
    infIndices = isinf(LG_EMG);
    % Replace NaN and infinite values with zeros
    LG_EMG(nanIndices) = 0;
    LG_EMG(infIndices) = 0;

   % Lowpass and highpass filter information
    fs = 4000; % Sampling frequency in Hz
    f_lowpass = 50; % Desired cutoff frequency for lowpass filter in Hz
    f_highpass = 0.5; % Desired cutoff frequency for highpass filter in Hz
    order = 2; % Filter order

    % Design lowpass and highpass Butterworth filters
    [b_lowpass, a_lowpass] = butter(order, f_lowpass / (fs / 2), 'low');
    [b_highpass, a_highpass] = butter(order, f_highpass / (fs / 2), 'high');

    % Apply the lowpass filter to EMG 1 and 2, SONO, and Buckle
    filtered_Combined_Tendon_lowpass = filtfilt(b_lowpass, a_lowpass, Combined_Tendon);
    filtered_LG_EMG_lowpass = filtfilt(b_lowpass, a_lowpass, LG_EMG);
    filtered_PL_EMG_lowpass = filtfilt(b_lowpass, a_lowpass, PL_EMG);
    filtered_LG_Sono_lowpass = filtfilt(b_lowpass, a_lowpass, LG_Sono);
    filtered_PL_Sono_lowpass = filtfilt(b_lowpass, a_lowpass, PL_Sono);

    % Apply the highpass filter to the lowpass filtered signals
    filtered_Combined_Tendon_highpass = filtfilt(b_highpass, a_highpass, filtered_Combined_Tendon_lowpass);
    filtered_LG_EMG_highpass = filtfilt(b_highpass, a_highpass, filtered_LG_EMG_lowpass);
    filtered_PL_EMG_highpass = filtfilt(b_highpass, a_highpass, filtered_PL_EMG_lowpass);
    filtered_LG_Sono_highpass = filtfilt(b_highpass, a_highpass, filtered_LG_Sono_lowpass);
    filtered_PL_Sono_highpass = filtfilt(b_highpass, a_highpass, filtered_PL_Sono_lowpass);

    % Add the filtered columns to the table
    data.filtered_Combined_Tendon = filtered_Combined_Tendon_highpass;
    data.filtered_LG_EMG = filtered_LG_EMG_highpass;
    data.filtered_PL_EMG = filtered_PL_EMG_highpass;
    data.filtered_LG_Sono = filtered_LG_Sono_highpass;
    data.filtered_PL_Sono = filtered_PL_Sono_highpass;

    % Step 3: Store the cleaned and filtered data
    cleanedDataCell{i} = data;
end

% Loop through each file in the folder
for i = 1:numel(fileList)
    % Read the current CSV file
    cleanedData = cleanedDataCell{i};
    
    % Remove NaN values from all columns except 'Sono LG'
    nanColumns = cleanedData.Properties.VariableNames(5:end);
    sononColumn = strcmp(nanColumns, 'PLSono');
    sononColumnIndex = find(sononColumn);
    
    for j = 1:numel(nanColumns)
        if j ~= sononColumnIndex
            columnData = cleanedData.(nanColumns{j});
            nanMask = isnan(columnData);
            cleanedData = cleanedData(~nanMask, :);
        end
    end

    % Store the cleaned data in the cell array
    cleanedDataCell{i} = cleanedData;
end

% Concatenate all cleaned data together vertically
mergedData = vertcat(cleanedDataCell{:});

% Write the merged data to a new CSV file
mergedFilePath = sprintf('D:\\%s\\kin and an\\%s\\%s.csv', animal_number, sand_or_solid, file_name);
writetable(mergedData, mergedFilePath);

% Specify the folders containing the sand and solid files
folder_path_sand = sprintf('D:\\%s\\BORIS\\%s\\400 hz', animal_number, sand_or_solid);

% Get a list of files in the sand folder
file_list_sand = dir(fullfile(folder_path_sand, '*.xlsx'));

% Initialize cell arrays to store data from each folder
data_cell_sand = cell(length(file_list_sand), 1);

% Read and store data from files in the sand folder
for i = 1:length(file_list_sand)
    file_path = fullfile(folder_path_sand, file_list_sand(i).name);
    data = readtable(file_path);
    data_cell_sand{i} = data;
end

% Vertically concatenate the data from all sand files
concatenated_data_sand = vertcat(data_cell_sand{:});

% Define the file paths for the concatenated sand and solid files
output_file_path_sand = sprintf('D:\\%s\\BORIS\\%s\\%s.csv', animal_number, sand_or_solid, file_name);

% Write the concatenated sand data to a new CSV file
writetable(concatenated_data_sand, output_file_path_sand);

% read in kinematics and boris data files
kinematics = readtable(sprintf('D:\\%s\\kin and an\\%s\\%s.csv', animal_number, sand_or_solid, file_name));
boris = readtable(sprintf('D:\\%s\\BORIS\\%s\\%s.csv', animal_number, sand_or_solid, file_name));

% combine the files
merged = horzcat(kinematics, boris);

% Define the file name and path where you want to save the CSV file
file_name = sprintf('merged_%s.csv', file_name);
file_path = sprintf('D:\\%s\\merged', animal_number);

% Combine the file name and path
full_file_path = fullfile(file_path, file_name);

% Export the merged table to a CSV file
 writetable(merged, full_file_path);

disp('Merging and saving completed successfully.');