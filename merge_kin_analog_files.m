clc
clear

% Set the folder where CSV files are located
folderPath = 'D:\DD121 - InVivo sand\kin and an\solid';

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
    EMG_1 = data.EMG1LG;
    EMG_2 = data.EMG2LG;
    SONO = data.SonoLG;
    Buckle = data.Buckle;

    % Remove mean from EMG_1, EMG_2, Buckle
    EMG_1 = detrend(EMG_1);
    EMG_2 = detrend(EMG_2);

    % Rectify absolute values
    EMG_1 = abs(EMG_1);
    EMG_2 = abs(EMG_2);

    % Invert the Buckle data
    Buckle = -Buckle;

    % Check for NaN values in SONO
    nanIndices = isnan(SONO);
    % Check for infinite values in SONO
    infIndices = isinf(SONO);
    % Replace NaN and infinite values with zeros
    SONO(nanIndices) = 0;
    SONO(infIndices) = 0;

    % Check for NaN values in EMG_1
    nanIndices = isnan(EMG_1);
    % Check for infinite values in EMG_1
    infIndices = isinf(EMG_1);
    % Replace NaN and infinite values with zeros
    EMG_1(nanIndices) = 0;
    EMG_1(infIndices) = 0;

   % Lowpass and highpass filter information
    fs = 4000; % Sampling frequency in Hz
    f_lowpass = 50; % Desired cutoff frequency for lowpass filter in Hz
    f_highpass = 0.5; % Desired cutoff frequency for highpass filter in Hz
    order = 2; % Filter order

    % Design lowpass and highpass Butterworth filters
    [b_lowpass, a_lowpass] = butter(order, f_lowpass / (fs / 2), 'low');
    [b_highpass, a_highpass] = butter(order, f_highpass / (fs / 2), 'high');

    % Apply the lowpass filter to EMG 1 and 2, SONO, and Buckle
    filtered_EMG_1_lowpass = filtfilt(b_lowpass, a_lowpass, EMG_1);
    filtered_EMG_2_lowpass = filtfilt(b_lowpass, a_lowpass, EMG_2);
    filtered_SONO_lowpass = filtfilt(b_lowpass, a_lowpass, SONO);
    filtered_Buckle_lowpass = filtfilt(b_lowpass, a_lowpass, Buckle);

    % Apply the highpass filter to the lowpass filtered signals
    filtered_EMG_1_highpass = filtfilt(b_highpass, a_highpass, filtered_EMG_1_lowpass);
    filtered_EMG_2_highpass = filtfilt(b_highpass, a_highpass, filtered_EMG_2_lowpass);
    filtered_SONO_highpass = filtfilt(b_highpass, a_highpass, filtered_SONO_lowpass);
    filtered_Buckle_highpass = filtfilt(b_highpass, a_highpass, filtered_Buckle_lowpass);

    % Add the filtered columns to the table
    data.filtered_EMG_1 = filtered_EMG_1_highpass;
    data.filtered_EMG_2 = filtered_EMG_2_highpass;
    data.filtered_SONO = filtered_SONO_highpass;
    data.filtered_Buckle = filtered_Buckle_highpass;

% Step 3: Store the cleaned and filtered data
cleanedDataCell{i} = data;

end

% Loop through each file in the folder
for i = 1:numel(fileList)
    % Read the current CSV file
    cleanedData = cleanedDataCell{i};
    
    % Remove NaN values from all columns except 'Sono LG'
    nanColumns = cleanedData.Properties.VariableNames(5:end);
    sononColumn = strcmp(nanColumns, 'SonoLG');
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
mergedFilePath = 'D:\DD121 - InVivo sand\kin and an\DD121_solid!.csv';
writetable(mergedData, mergedFilePath);

disp('Merging and saving completed successfully.');
