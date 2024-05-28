clear
clc

% read in data from correct filepath
 data_nan = readtable("files for MATLAB/merged/DD_old/DD139_slope_all trials_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD133_tread_slopes_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD133_slopes_merged.csv");

% data_nan = readtable("files for MATLAB/merged/DD130_slopes_01_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD130_slopes_02_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD135_slopes_02_merged.csv");

% data_nan = readtable("files for MATLAB/merged/SCDD03_hard.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD03_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand2.csv");

% data_nan = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\DD118B_sand_merged_data.csv");

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

%high pass filter information
fs = 500; % Sampling frequency in Hz
f_signal = 2; % Frequency of the signal in Hz

% Design a highpass Butterworth filter
order = 4; % Filter order
f_cutoff = 30; % Cutoff frequency in Hz
[d, c] = butter(order, f_cutoff/(fs/2), 'high');

% Apply the high and low pass filter to EMG 1 and 2
filt_EMG_1 = filtfilt(b, a, EMG_1);
filt_EMG_2 = filtfilt(b, a, EMG_2);
filtered_EMG_1 = filtfilt(d, c, filt_EMG_1);
filtered_EMG_2 = filtfilt(d, c, filt_EMG_2);

% apply the low pass filter to SONO and Buckle
filtered_SONO = filtfilt(b, a, SONO);
filtered_Buckle = filtfilt(b, a, Buckle);

% assign column names and convert from array2table
filt_EMG_1_table = array2table(filtered_EMG_1, 'VariableNames', {'filtered_EMG_1'});
filt_EMG_2_table = array2table(filtered_EMG_2, 'VariableNames', {'filtered_EMG_2'});
filt_SONO_table = array2table(filtered_SONO, 'VariableNames', {'filtered_SONO'});
filt_Buckle_table = array2table(filtered_Buckle, 'VariableNames', {'filtered_Buckle'});

% Concatenate the original table (DD139) with the new tables containing filtered EMG data
data_updated = [data, filt_SONO_table, filt_Buckle_table, filt_EMG_1_table, filt_EMG_2_table];

% filter data to get only good hops (=1) and to get jumps not affected by wedge (=0)
filtered_data = data_updated(data_updated.SteadyHopping == 1 & data_updated.Wedge_WedgeBack_ == 0 & data_updated.Wedge_WedgeFront_ == 0, :);

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
    time= linspace(-10, 110, numRows)';
    a.time = time;
    cycles_all{i} = a;
end

% Create a cell array to store the new tables
cycles = cell(size(cycles_all));

% Define the relevant columns
selectedColumns = {'time', 'SONO_LG', 'Buckle', 'EMG_1_LG', 'EMG_2_LG', 'filtered_SONO', 'filtered_Buckle', 'filtered_EMG_1', 'filtered_EMG_2'};

% Iterate through each table in the 'cycles_all' cell array
for i = 1:numel(cycles_all)
    % Select the relevant columns from the current table
    current_table = cycles_all{i};
    cycles{i} = current_table(:, selectedColumns);
end

% gets number of tables in cell array
 num_hop_cycles = numel(cycles);
% num_hop_cycles = 3;

% Plot EMG and force plate data with respect to pseudotime
figure(1)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.filtered_EMG_1;
    
    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
    
    xlabel("HOP CYCLE PERCENT")
    ylabel("EMG 1 LG (Volts)")
    title("Filtered EMG 1 vs Hop Cycle Percent")
end
xlim([-10, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(2)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.filtered_EMG_2;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("EMG 2 LG (Volts)")
    title("Filtered EMG 2 vs Hop Cycle Percent")
  
    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(3)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.filtered_SONO;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("SONO LG (Volts)")
    title("Filtered SONO LG vs Hop Cycle Percent")
    hold on

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(4)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.filtered_Buckle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Buckle (Volts)")
    title("Filtered Buckle vs Hop Cycle Percent")
    hold on;

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(5)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.EMG_1_LG;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("EMG 1 LG")
    title("Unfiltered EMG 1 vs Hop Cycle Percent")
    hold on

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(6)
hold on
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.EMG_2_LG;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("EMG 2 LG")
    title("Unfiltered EMG 2 vs Hop Cycle Percent")
    hold on

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(7)
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.SONO_LG;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("SONO LG (Volts)")
    title("Unfiltered SONO LG vs Hop Cycle Percent")
    hold on

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

figure(8)
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.Buckle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Buckle (Volts)")
    title("Unfiltered Buckle vs Hop Cycle Percent")
    hold on;

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
xlim([0, 100]);
hold off

% mean and std of EMG data
% Step 1: Compute the mean pseudotime for each cycle
pseudotime_mean_per_cycle = cellfun(@(data) mean(data.time), cycles);

% Step 2: Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 9; % Assuming you have 9 columns of interest
aligned_data_matrix = zeros(num_pseudotime_points, num_columns);

for i = 1:num_hop_cycles
    current_data = cycles{i};
    time_vector = current_data{:, 'time'};
    columns_data = table2array(current_data(:, 1:end));

    % Resample the data to align with a common pseudotime range
    pseudotime_range = linspace(time_vector(1), time_vector(end), num_pseudotime_points);
    aligned_data = interp1(time_vector, columns_data, pseudotime_range, 'spline');

    % Add the aligned data to the combined matrix
    aligned_data_matrix = aligned_data_matrix + aligned_data;
end

% Step 3: Compute the mean of each column over all cycles at corresponding pseudotime points
mean_values_array = aligned_data_matrix / num_hop_cycles;
mean_values=array2table(mean_values_array, 'VariableNames', {'time', 'SONO_LG', 'Buckle', 'EMG_1_LG', 'EMG_2_LG', 'filtered_SONO', 'filtered_Buckle', 'filtered_EMG_1', 'filtered_EMG_2'});

% Calculating the standard deviation of each column across all cycles at corresponding pseudotime points
std_values = std(mean_values, 0, 1);

% plotting mean and standard deviations of EMG data
figure(9)
errorbar(mean_values.time, mean_values.filtered_SONO, std_values.filtered_SONO, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_SONO, -std_values.filtered_SONO, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("SONO LG (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered SONO LG vs Hop Cycle Percent')
xlim([0, 100]);

figure(10)
errorbar(mean_values.time, mean_values.filtered_Buckle, std_values.filtered_Buckle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_Buckle, -std_values.filtered_Buckle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Buckle (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered Buckle vs Hop Cycle Percent')
xlim([0, 100]);

figure(11)
errorbar(mean_values.time, mean_values.filtered_EMG_1, std_values.filtered_EMG_1, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_EMG_1, -std_values.filtered_EMG_1, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("EMG 1 LG (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered EMG 1 vs Hop Cycle Percent')
xlim([0, 100]);

figure(12)
errorbar(mean_values.time, mean_values.filtered_EMG_2, std_values.filtered_EMG_2, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_EMG_2, -std_values.filtered_EMG_2, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("EMG 2 LG (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered EMG 2 vs Hop Cycle Percent')
xlim([0, 100]);

% Define the callback function for the line click event
function lineClicked(src, ~)
    cycle_num_str = get(src, 'DisplayName');
    cycle_num = str2double(regexp(cycle_num_str, '\d+', 'match'));
    if ~isnan(cycle_num)
        disp(['Clicked on cycle number: ' num2str(cycle_num)]);
    end
end