function [EMG_back_w_cycles] = EMG_back_w_cycles(data_nan)

% remove all NAN cells
data=rmmissing(data_nan); 

% data that needs to be filtered
EMG_1 = data.EMG_1_LG;
EMG_2 = data.EMG_2_LG;
SONO = data.SONO_LG;
Buckle = data.Buckle;

% remove mean from from EMG_1, EMG_2
EMG_1 = detrend(EMG_1);
EMG_2 = detrend(EMG_2);

% recity absolute values
EMG_1 = abs(EMG_1);
EMG_2 = abs(EMG_2);

% lowpass filter information
fs = 500; % Sampling frequency in Hz
f = 30; % Desired cutoff frequency in Hz

% Design a lowpass Butterworth filter
order = 4 ; % Filter order
[b, a] = butter(order, f/(fs/2), 'low');

% Apply the low pass filter to EMG 1 and 2, SONO, and Buckle
filtered_EMG_1 = filtfilt(b, a, EMG_1);
filtered_EMG_2 = filtfilt(b, a, EMG_2);
filtered_SONO = filtfilt(b, a, SONO);
filtered_Buckle = filtfilt(b, a, Buckle);

% assign column names and convert from array2table
filt_EMG_1_table = array2table(filtered_EMG_1, 'VariableNames', {'filtered_EMG_1'});
filt_EMG_2_table = array2table(filtered_EMG_2, 'VariableNames', {'filtered_EMG_2'});
filt_SONO_table = array2table(filtered_SONO, 'VariableNames', {'filtered_SONO'});
filt_Buckle_table = array2table(filtered_Buckle, 'VariableNames', {'filtered_Buckle'});

% Concatenate the original table (DD139) with the new tables containing filtered EMG data
data_updated = [data, filt_SONO_table, filt_Buckle_table, filt_EMG_1_table, filt_EMG_2_table];

if ismember('Wedge', data_updated.Properties.VariableNames)
filtered_data = data_updated(data_updated.SteadyHopping == 1 & (data_updated.Wedge_WedgeBack_ == 1 | data_updated.Wedge == 1), :);
else
filtered_data = data_updated(data_updated.SteadyHopping == 1 & (data_updated.Wedge_WedgeBack_ == 1), :);
end

% Find the indices where the column equals 1
indices = find(filtered_data.Landing(:, 1) == 1);

% Split the data into cycles
cycles_all = cell(length(indices)-1, 1);
% cycles_all = 1
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = filtered_data(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(0, 100, numRows)';
    a.time = time;
    cycles_all{i} = a;
end

% Create a cell array to store the new tables
EMG_back_w_cycles = cell(size(cycles_all));

% Define the relevant columns
selectedColumns = {'time', 'SONO_LG', 'Buckle', 'EMG_1_LG', 'EMG_2_LG', 'filtered_SONO', 'filtered_Buckle', 'filtered_EMG_1', 'filtered_EMG_2'};

% Iterate through each table in the 'cycles_all' cell array
for i = 1:numel(cycles_all)
    % Select the relevant columns from the current table
    current_table = cycles_all{i};
    EMG_back_w_cycles{i} = current_table(:, selectedColumns);
end
