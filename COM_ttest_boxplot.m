clc
clear

% Sand
%  data_sand = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118A\DD118A_sand_merged_data.csv");
%  data_sand = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD121\DD121_sand_merged_data.csv");
%  data_sand = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_sand_merged_data.csv");
%  data_sand = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\DD118B_sand_merged_data.csv");

%  data_sand = readtable('D:\DD119_all\merged\merged_DD119_sand.csv');
%  data_sand = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_sand.csv');
%  data_sand = readtable('D:\DD118B - invivo sand\merged\merged_DD118B_sand.csv');
%  data_sand = readtable('D:\DD118A_all\merged\merged_DD118A_sand.csv');

% data_sand = readtable('D:\DD118B - invivo sand\merged\merged_DD118B_sand!.csv');
% data_sand = readtable('D:\DD118A_all\merged\merged_DD118A_sand!.csv');
% data_sand = readtable('D:\DD119_all\merged\merged_DD119_sand.csv');
data_sand = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_sand5!.csv');

% Solid
%  data_solid = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118A\DD118A_solid_merged_data.csv");
%  data_solid = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD121\DD121_solid_merged_data.csv");
%  data_solid = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_solid_merged_data.csv");
%  data_solid = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD118B\DD118B_solid_merged_data.csv");

%  data_solid = readtable('D:\DD119_all\merged\merged_DD119_solid.csv');
%  data_solid = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_solid.csv');
%  data_solid = readtable('D:\DD118B - invivo sand\merged\merged_DD118B_solid.csv');
%  data_solid = readtable('D:\DD118A_all\merged\merged_DD118A_solid.csv');

% data_solid = readtable('D:\DD118B - invivo sand\merged\merged_DD118B_solid!2.csv');
% data_solid = readtable('D:\DD118A_all\merged\merged_DD118A_solid!.csv');
% data_solid = readtable('D:\DD119_all\merged\merged_DD119_solid!.csv');
 data_solid = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_solid!.csv');

[data_sand] = kin_filt(data_sand);
[data_solid] = kin_filt(data_solid);

% indexing
[upperback_sand, hip_sand, knee_sand, ankle_sand, foot_sand, toetip_sand, eye_sand] = filt_krat_table(data_sand);
[upperback_solid, hip_solid, knee_solid, ankle_solid, foot_solid, toetip_solid, eye_solid] = filt_krat_table(data_solid);

% nameing columns x and y for eye amd hip
eye_sand = array2table(eye_sand, 'VariableNames', {'x_eye', 'y_eye'});
hip_sand = array2table(hip_sand, 'VariableNames', {'x_hip', 'y_hip'});

eye_solid = array2table(eye_solid, 'VariableNames', {'x_eye', 'y_eye'});
hip_solid = array2table(hip_solid, 'VariableNames', {'x_hip', 'y_hip'});

% Scale factor of 2.5 mm/pixel
scale_factor = 2.5;

% Divide the x and y coordinates of eye by the scale factor
eye_sand.x_eye = eye_sand.x_eye / scale_factor;
eye_sand.y_eye = eye_sand.y_eye / scale_factor;

eye_solid.x_eye = eye_solid.x_eye / scale_factor;
eye_solid.y_eye = eye_solid.y_eye / scale_factor;

% Divide the x and y coordinates of hip by the scale factor
hip_sand.x_hip = hip_sand.x_hip / scale_factor;
hip_sand.y_hip = hip_sand.y_hip / scale_factor;

hip_solid.x_hip = hip_solid.x_hip / scale_factor;
hip_solid.y_hip = hip_solid.y_hip / scale_factor;

% Calculate the COM coordinates for the filtered data
com_x_sand = (eye_sand.x_eye + hip_sand.x_hip) / 2;
com_y_sand = (eye_sand.y_eye + hip_sand.y_hip) / 2;

com_x_solid = (eye_solid.x_eye + hip_solid.x_hip) / 2;
com_y_solid = (eye_solid.y_eye + hip_solid.y_hip) / 2;

% Initialize an array to store the internal angles
pitch_angle_sand = zeros(height(hip_sand), 1);
pitch_angle_solid = zeros(height(hip_solid), 1);

% Calculate internal angles for each pair of coordinates
for i = 1:height(hip_sand)
    x_eye = eye_sand.x_eye(i);
    y_eye = eye_sand.y_eye(i);
    
    x_hip = hip_sand.x_hip(i);
    y_hip = hip_sand.y_hip(i);
    
    % Calculate the angle in radians using the arctangent function (atan2)
    pitch_angle_sand(i) = atan2(y_hip - y_eye, x_hip - x_eye);
end

for i = 1:height(hip_solid)
    x_eye = eye_solid.x_eye(i);
    y_eye = eye_solid.y_eye(i);
    
    x_hip = hip_solid.x_hip(i);
    y_hip = hip_solid.y_hip(i);
    
    % Calculate the angle in radians using the arctangent function (atan2)
    pitch_angle_solid(i) = atan2(y_hip - y_eye, x_hip - x_eye);
end

% Convert the angles to degrees if needed
angles_array_sand = rad2deg(pitch_angle_sand);
angles_array_solid = rad2deg(pitch_angle_solid);

% Create a table for the COM coordinates
com_table_sand = table(com_x_sand, com_y_sand, angles_array_sand, 'VariableNames', {'com_x', 'com_y', 'pitch_angle'});
com_table_solid = table(com_x_solid, com_y_solid, angles_array_solid, 'VariableNames', {'com_x', 'com_y', 'pitch_angle'});

% Negate the 'com_y' and 'pitch_angle' columns
com_table_sand.pitch_angle = -com_table_sand.pitch_angle;
com_table_sand.com_y = -com_table_sand.com_y;

com_table_solid.pitch_angle = -com_table_solid.pitch_angle;
com_table_solid.com_y = -com_table_solid.com_y;

% Convert the 'Landing' and 'SteadyHopping' variables to tables
landing_column_sand = table(data_sand.Landing, 'VariableNames', {'Landing'});
steadyhopping_column_sand = table(data_sand.SteadyHopping, 'VariableNames', {'SteadyHopping'});

landing_column_solid = table(data_solid.Landing, 'VariableNames', {'Landing'});
steadyhopping_column_solid = table(data_solid.SteadyHopping, 'VariableNames', {'SteadyHopping'});

% Add the 'Landing' and 'SteadyHopping' columns to com_table
com_table_sand = [com_table_sand, landing_column_sand, steadyhopping_column_sand];
com_table_solid = [com_table_solid, landing_column_solid, steadyhopping_column_solid];

% Find the indices where the column equals 1
indices_sand = find(com_table_sand.Landing(:, 1) == 1);
indices_solid = find(com_table_solid.Landing(:, 1) == 1);

% Split the data into cycles
cycles_sand = cell(length(indices_sand)-1, 1);
for i = 1:length(indices_sand)-1
    startIdx = indices_sand(i);
    endIdx = indices_sand(i+1)-1;
    a = com_table_sand(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles_sand{i} = a;
end

cycles_solid = cell(length(indices_solid)-1, 1);
for i = 1:length(indices_solid)-1
    startIdx = indices_solid(i);
    endIdx = indices_solid(i+1)-1;
    a = com_table_solid(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles_solid{i} = a;
end

% Define column names for indexing
time_column_name_sand = 'time';
com_y_column_name_sand = 'com_y';
pitch_angle_column_name_sand = 'pitch_angle';

time_column_name_solid = 'time';
com_y_column_name_solid = 'com_y';
pitch_angle_column_name_solid = 'pitch_angle';

% Group the column names in a cell array
column_names_sand = {time_column_name_sand, com_y_column_name_sand, pitch_angle_column_name_sand};
column_names_solid = {time_column_name_solid, com_y_column_name_solid, pitch_angle_column_name_solid};

% Iterate over columns to get angles and pseudotime
% Grouped each cycle with this data in cells
cell_array_sand = cycles_sand;
num_tables_sand = numel(cell_array_sand);
num_columns_sand = numel(column_names_sand);
data_cell_sand = cell(num_tables_sand, 1); 

cell_array_solid = cycles_solid;
num_tables_solid = numel(cell_array_solid);
num_columns_solid = numel(column_names_solid);
data_cell_solid = cell(num_tables_solid, 1); 

for i = 1:num_tables_sand
    current_table_sand = cell_array_sand{i};
    columns_data_sand = cell(1, num_columns_sand);
    for j = 1:num_columns_sand
        a_column_data_sand = current_table_sand.(column_names_sand{j});
        columns_data_sand{j} = a_column_data_sand;
    end
    data_cell_sand{i} = [columns_data_sand{:}];
end

for i = 1:num_tables_solid
    current_table_solid = cell_array_solid{i};
    columns_data_solid = cell(1, num_columns_solid);
    for j = 1:num_columns_solid
        a_column_data_solid = current_table_solid.(column_names_solid{j});
        columns_data_solid{j} = a_column_data_solid;
    end
    data_cell_solid{i} = [columns_data_solid{:}];
end

% defining and plotting the COM
cycles_array_sand = data_cell_sand;
num_hop_cycles_sand = height(cycles_array_sand);

cycles_array_solid = data_cell_solid;
num_hop_cycles_solid = height(cycles_array_solid);

% Iterate over each cell in cycles_array
for i = 1:numel(cycles_array_sand)
    current_data_sand = cycles_array_sand{i};  % Get the current cell

    % Choose the column to subtract the mean from (e.g., column 1)
    column_to_subtract_mean_from = 2;

    % Calculate the mean of the selected column
    mean_value_sand = mean(current_data_sand(:, column_to_subtract_mean_from));

    % Subtract the mean from every point in the selected column
    current_data_sand(:, column_to_subtract_mean_from) = current_data_sand(:, column_to_subtract_mean_from) - mean_value_sand;

    % Update the current cell in cycles_array
    cycles_array_sand{i} = current_data_sand;
end

for i = 1:numel(cycles_array_solid)
    current_data_solid = cycles_array_solid{i};  % Get the current cell

    % Choose the column to subtract the mean from (e.g., column 1)
    column_to_subtract_mean_from = 2;

    % Calculate the mean of the selected column
    mean_value_solid = mean(current_data_solid(:, column_to_subtract_mean_from));

    % Subtract the mean from every point in the selected column
    current_data_solid(:, column_to_subtract_mean_from) = current_data_solid(:, column_to_subtract_mean_from) - mean_value_solid;

    % Update the current cell in cycles_array
    cycles_array_solid{i} = current_data_solid;
end

% Angular momentum calculations
% Define the sampling rate
sampling_rate = 100;

% Loop over each double in the cell array
for i = 1:numel(cycles_array_solid)
    % Get the current double
    current_double = cycles_array_solid{i};
    
    % Extract the third column (assuming it's the column you want to differentiate)
    column_to_differentiate = current_double(:, 3);
    
    % Compute the derivative using diff and divide by the sampling rate
    derivative_column = diff(column_to_differentiate) / sampling_rate;
    
    % Add NaN as the first element to maintain alignment
    derivative_column = [NaN; derivative_column];
    
    % Add the derivative column as the fourth column in the current double
    current_double(:, 4) = derivative_column;
    
    % Update the current double in the cell array
    cycles_array_solid{i} = current_double;
end

% Loop over each double in the cell array
for i = 1:numel(cycles_array_sand)
    % Get the current double
    current_double = cycles_array_sand{i};
    
    % Extract the third column (assuming it's the column you want to differentiate)
    column_to_differentiate = current_double(:, 3);
    
    % Compute the derivative using diff and divide by the sampling rate
    derivative_column = diff(column_to_differentiate) / sampling_rate;
    
    % Add NaN as the first element to maintain alignment
    derivative_column = [NaN; derivative_column];
    
    % Add the derivative column as the fourth column in the current double
    current_double(:, 4) = derivative_column;
    
    % Update the current double in the cell array
    cycles_array_sand{i} = current_double;
end

% Filter angular momentum column
% Define the lowpass filter parameters
fs = 100; % Sampling frequency in Hz
fc = 5;   % Cutoff frequency in Hz

% Design a lowpass Butterworth filter
order = 4; % Filter order
[b, a] = butter(order, fc / (fs / 2), 'low');

% Loop over each double in the cell array (solid data)
for i = 1:numel(cycles_array_solid)
    % Get the current double
    current_double = cycles_array_solid{i};
    
    % Extract the fourth column (assumed to be the angular momentum column)
    angular_momentum_column = current_double(:, 4);
    
    % Check if the input data contains non-finite values
    if any(~isfinite(angular_momentum_column))
        continue;
    end
    
    % Apply lowpass filtering
    filtered_angular_momentum = filtfilt(b, a, angular_momentum_column);
    
    % Update the fourth column with the filtered angular momentum
    current_double(:, 4) = filtered_angular_momentum;
    
    % Update the current double in the cell array
    cycles_array_solid{i} = current_double;
end

% Loop over each double in the cell array (sand data)
for i = 1:numel(cycles_array_sand)
    % Get the current double
    current_double = cycles_array_sand{i};
    
    % Extract the fourth column (assumed to be the angular momentum column)
    angular_momentum_column = current_double(:, 4);
    
    % Check if the input data contains non-finite values
    if any(~isfinite(angular_momentum_column))
        continue;
    end
    
    % Apply lowpass filtering
    filtered_angular_momentum = filtfilt(b, a, angular_momentum_column);
    
    % Update the fourth column with the filtered angular momentum
    current_double(:, 4) = filtered_angular_momentum;
    
    % Update the current double in the cell array
    cycles_array_sand{i} = current_double;
end



% Define the column names for each table (struct) within a_data_cell
column_names_sand = {'time', 'com', 'pitch_angle', 'angular_momentum'};
column_names_solid = {'time', 'com', 'pitch_angle', 'angular_momentum'};

% Loop through each table (struct) in a_data_cell and set the column names
for i = 1:numel(cycles_array_sand)
    data_cell_sand{i} = array2table(cycles_array_sand{i}, 'VariableNames', column_names_sand);
end

for i = 1:numel(cycles_array_solid)
    data_cell_solid{i} = array2table(cycles_array_solid{i}, 'VariableNames', column_names_solid);
end


% Step 1: Compute the mean pseudotime for each cycle for sand data
pseudotime_mean_per_cycle_sand = cellfun(@(data_sand) mean(data_sand), data_cell_sand, 'UniformOutput', false);

% Step 2: Use the mean pseudotime as a reference to align the sand data
num_pseudotime_points = 100;
num_columns = 4; % Assuming you have 4 columns of interest (excluding the time column)
aligned_data_matrix_sand = zeros(num_pseudotime_points, num_columns);

for i = 1:num_tables_sand
    current_data_sand = data_cell_sand{i};
    time_vector_sand = current_data_sand{:, 'time'};
    columns_data_sand = table2array(current_data_sand(:, 1:end));

    % Resample the sand data to align with a common pseudotime range
    pseudotime_range_sand = linspace(time_vector_sand(1), time_vector_sand(end), num_pseudotime_points);
    aligned_data_sand = interp1(time_vector_sand, columns_data_sand, pseudotime_range_sand, 'spline');

    % Add the aligned sand data to the combined matrix
    aligned_data_matrix_sand = aligned_data_matrix_sand + aligned_data_sand;
end

% Step 3: Compute the mean of each column over all cycles at corresponding pseudotime points for sand data
mean_values_array_sand = aligned_data_matrix_sand / num_tables_sand;
mean_values_sand = array2table(mean_values_array_sand, 'VariableNames', {'Time', 'COM', 'pitch_angle', 'angular_momentum'});

% Duplicated for Solid Data

% Step 1: Compute the mean pseudotime for each cycle for solid data
pseudotime_mean_per_cycle_solid = cellfun(@(data_solid) mean(data_solid), data_cell_solid, 'UniformOutput', false);

% Step 2: Use the mean pseudotime as a reference to align the solid data
num_pseudotime_points = 100;
num_columns = 4; % Assuming you have 4 columns of interest (excluding the time column)
aligned_data_matrix_solid = zeros(num_pseudotime_points, num_columns);

for i = 1:num_tables_solid
    current_data_solid = data_cell_solid{i};
    time_vector_solid = current_data_solid{:, 'time'};
    columns_data_solid = table2array(current_data_solid(:, 1:end));

    % Resample the solid data to align with a common pseudotime range
    pseudotime_range_solid = linspace(time_vector_solid(1), time_vector_solid(end), num_pseudotime_points);
    aligned_data_solid = interp1(time_vector_solid, columns_data_solid, pseudotime_range_solid, 'spline');

    % Add the aligned solid data to the combined matrix
    aligned_data_matrix_solid = aligned_data_matrix_solid + aligned_data_solid;
end

% Step 3: Compute the mean of each column over all cycles at corresponding pseudotime points for solid data
mean_values_array_solid = aligned_data_matrix_solid / num_tables_solid;
mean_values_solid = array2table(mean_values_array_solid, 'VariableNames', {'Time', 'COM', 'pitch_angle', 'angular_momentum'});

% Calculating the standard deviation of each column across all cycles at corresponding pseudotime points
std_values_sand = std(mean_values_sand, 0, 1);
std_values_solid = std(mean_values_solid, 0, 1);

% Calculate the average displacement across all cycles
average_displacement_sand = mean(mean_values_sand.COM);
average_angle_sand = mean(mean_values_sand.pitch_angle);

average_displacement_solid = mean(mean_values_solid.COM);
average_angle_solid = mean(mean_values_solid.pitch_angle);

% Calculate the maximum and minimum displacement values for "Solid" and "Sand"
max_displacement_solid = max(mean_values_solid.COM) - min(mean_values_solid.COM);
max_displacement_sand = max(mean_values_sand.COM) - min(mean_values_sand.COM);

% Initialize cell arrays to store max and min values
max_values_sand_com = cell(length(data_cell_sand), 1);
min_values_sand_com = cell(length(data_cell_sand), 1);
max_values_sand_angle = cell(length(data_cell_sand), 1);
min_values_sand_angle = cell(length(data_cell_sand), 1);

max_values_solid_com = cell(length(data_cell_solid), 1);
min_values_solid_com = cell(length(data_cell_solid), 1);
max_values_solid_angle = cell(length(data_cell_solid), 1);
min_values_solid_angle = cell(length(data_cell_solid), 1);

% Calculate max and min values for "Sand" data
for i = 1:length(data_cell_sand)
    current_cycle = data_cell_sand{i};
    max_values_sand_com{i} = max(current_cycle.com);  
    min_values_sand_com{i} = min(current_cycle.com);  
    max_values_sand_angle{i} = max(current_cycle.pitch_angle);  
    min_values_sand_angle{i} = min(current_cycle.pitch_angle); 
end

% Calculate max and min values for "Solid" data
for i = 1:length(data_cell_solid)
    current_cycle = data_cell_solid{i};
    max_values_solid_com{i} = max(current_cycle.com);  
    min_values_solid_com{i} = min(current_cycle.com); 
    max_values_solid_angle{i} = max(current_cycle.pitch_angle);  
    min_values_solid_angle{i} = min(current_cycle.pitch_angle); 
end

% Choose the number of cycles to compare
num_cycles_to_compare = min(length(max_values_sand_com), length(max_values_solid_com));

% Take the first N cycles from both datasets
max_values_sand_com = max_values_sand_com(1:num_cycles_to_compare);
min_values_sand_com = min_values_sand_com(1:num_cycles_to_compare);
max_values_solid_com = max_values_solid_com(1:num_cycles_to_compare);
min_values_solid_com = min_values_solid_com(1:num_cycles_to_compare);

max_values_sand_angle = max_values_sand_angle(1:num_cycles_to_compare);
min_values_sand_angle = min_values_sand_angle(1:num_cycles_to_compare);
max_values_solid_angle = max_values_solid_angle(1:num_cycles_to_compare);
min_values_solid_angle = min_values_solid_angle(1:num_cycles_to_compare);

% Convert cell arrays to numerical arrays
max_values_sand_com = cell2mat(max_values_sand_com);
min_values_sand_com = cell2mat(min_values_sand_com);
max_values_solid_com = cell2mat(max_values_solid_com);
min_values_solid_com = cell2mat(min_values_solid_com);

max_values_sand_angle = cell2mat(max_values_sand_angle);
min_values_sand_angle = cell2mat(min_values_sand_angle);
max_values_solid_angle = cell2mat(max_values_solid_angle);
min_values_solid_angle = cell2mat(min_values_solid_angle);

% Perform paired t-tests for max-min displacement and max-min pitch angle differences
[h_maxmin_com, p_maxmin_com] = ttest2(max_values_solid_com - min_values_solid_com, max_values_sand_com - min_values_sand_com);
[h_maxmin_angle, p_maxmin_angle] = ttest2(max_values_solid_angle - min_values_solid_angle, max_values_sand_angle - min_values_sand_angle);

fprintf('P-value for Maximum-Minimum Displacement (Solid vs. Sand): %.5f\n', p_maxmin_com);
fprintf('P-value for Maximum-Minimum Pitch Angle Difference (Solid vs. Sand): %.5f\n', p_maxmin_angle);

if p_maxmin_com < 0.05
    fprintf('Result for Maximum-Minimum Displacement is statistically significant (reject the null hypothesis)\n');
else
    fprintf('Result for Maximum-Minimum Displacement is not statistically significant (fail to reject the null hypothesis)\n');
end

if p_maxmin_angle < 0.05
    fprintf('Result for Maximum-Minimum Pitch Angle Difference is statistically significant (reject the null hypothesis)\n');
else
    fprintf('Result for Maximum-Minimum Pitch Angle Difference is not statistically significant (fail to reject the null hypothesis)\n');
end

% one_cycle = cycles_array_solid{1};
% one_cycle =array2table(one_cycle);
% one_cycle = renamevars(one_cycle, {'one_cycle1','one_cycle2','one_cycle3'}, {'time', 'com_y', 'pitch_angle'});

% figure(11)
% plot(one_cycle.time, one_cycle.com_y, 'LineWidth', 4)
% xlabel("HOP CYCLE PERCENT")
% ylabel("Center of Mass Displacement [mm]")
% title('Center of Mass Displacement on a Solid Surface vs Hop Cycle Percent')

 % num_hop_cycles_solid = 100;
 % num_hop_cycles_sand = 30;


figure(1)
for i = 1:num_hop_cycles_solid
 % Skip cycle 44
    if i == 44
        continue;
    end

    current_array = cycles_array_solid{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    x_values = current_table.time;
    y_values = current_table.com_y;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Center of Mass Displacement [mm]")
    title('Center of Mass Displacement on a Solid Surface vs Hop Cycle Percent')
    hold on;

     % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off


figure(2)
errorbar(mean_values_solid.Time, mean_values_solid.COM, std_values_solid.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values_solid.Time, mean_values_solid.COM, -std_values_solid.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Center of Mass Displacement [mm]")
legend('Mean', 'Std Dev')
title('Average Center of Mass Displacement on a Solid Surface vs Hop Cycle Percent')


figure(3)
for i = 1:num_hop_cycles_solid
 % Skip cycle 44
   if i == 36 || i == 25
        continue;
    end

    current_array = cycles_array_solid{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    x_values = current_table.time;
    y_values = current_table.pitch_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Pitch Angle [deg]")
    title('Pitch Angle on a Solid Surface vs Hop Cycle Percent')
    hold on;
% Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end

% Set the scale for the y-axis
    ylim([min(-10), max(25)]);

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off

figure(4)
errorbar(mean_values_solid.Time, mean_values_solid.pitch_angle, std_values_solid.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values_solid.Time, mean_values_solid.pitch_angle, -std_values_solid.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Pitch Angle [deg]")
legend('Mean', 'Std Dev')
title('Average Pitch Angle on a Solid Surface vs Hop Cycle Percent')


figure(5)
for i = 1:num_hop_cycles_sand
    % Skip cycle 29, 209
    if i == 29
        continue;
    end
    
    current_array = cycles_array_sand{i};
    current_table_wrong_names = array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    
    x_values = current_table.time;
    y_values = current_table.com_y;
    
    plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Center of Mass Displacement [mm]")
    title('Center of Mass Displacement on a Sand Surface vs Hop Cycle Percent')
    hold on;
end

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off

figure(6)
errorbar(mean_values_sand.Time, mean_values_sand.COM, std_values_sand.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values_sand.Time, mean_values_sand.COM, -std_values_sand.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Center of Mass Displacement [mm]")
legend('Mean', 'Std Dev')
title('Average Center of Mass Displacement on a Sand Surface vs Hop Cycle Percent')


figure(7)
for i = 1:num_hop_cycles_sand
 % Skip cycle 29, 209
    if i == 29 || i == 48 
        continue;
    end

    current_array = cycles_array_sand{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    x_values = current_table.time;
    y_values = current_table.pitch_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Pitch Angle [deg]")
    title('Pitch Angle on a Sand Surface vs Hop Cycle Percent')
    hold on;
% Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end

% Set the scale for the y-axis
    ylim([min(-10), max(25)]);

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off

figure(8)
errorbar(mean_values_sand.Time, mean_values_sand.pitch_angle, std_values_sand.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values_sand.Time, mean_values_sand.pitch_angle, -std_values_sand.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Pitch Angle [deg]")
legend('Mean', 'Std Dev')
title('Average Pitch Angle on a Sand Surface vs Hop Cycle Percent')


% Combine max-min displacement data for sand and solid into a single matrix
max_min_displacement_data = [max_values_sand_com - min_values_sand_com, max_values_solid_com - min_values_solid_com];

% Create a table and label the columns
max_min_displacement_table = array2table(max_min_displacement_data, 'VariableNames', {'sand_m_m', 'solid_m_m'});

mean_max_min_disp_sand = mean(max_min_displacement_table.sand_m_m);
mean_max_min_disp_solid = mean(max_min_displacement_table.solid_m_m);


%Combine max-min pitch angle data for sand and solid into a single matrix
max_min_pitch_angle_data = [max_values_sand_angle - min_values_sand_angle, max_values_solid_angle - min_values_solid_angle];

% Create a table and label the columns
max_min_angle_table = array2table(max_min_pitch_angle_data, 'VariableNames', {'sand_m_m', 'solid_m_m'});

mean_max_min_angle_sand = mean(max_min_angle_table.sand_m_m);
mean_max_min_angle_solid = mean(max_min_angle_table.solid_m_m);


figure(9);
violin(max_min_displacement_data, 'Labels', {'Sand', 'Solid'});
title('Max-Min Displacement Violin Plot (Sand vs. Solid)');
ylabel('Max-Min Displacement [mm]');
xlabel('Sand                                    Solid')
xticks([]);


figure(10);
violin(max_min_pitch_angle_data, 'Labels', {'Sand', 'Solid'});
title('Max-Min Pitch Angle Violin Plot (Sand vs. Solid)');
ylabel('Max-Min Pitch Angle [degrees]');
xlabel('Sand                                    Solid')
xticks([]);

figure(11)
for i = 1:num_hop_cycles_sand
 % Skip cycle 29, 209
    if i == 29 || i == 48 
        continue;
    end

    current_array = cycles_array_sand{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    x_values = current_table.time;
    y_values = current_table.angular_momentum;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Velocity [deg/hop cycle percent]")
    title('Velocity on a Sand Surface vs Hop Cycle Percent')
    hold on;
% Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off

figure(12)
for i = 1:num_hop_cycles_solid
 % Skip cycle 29, 209
    if i == 36 || i == 25 
        continue;
    end

    current_array = cycles_array_solid{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3', 'current_array4'}, {'time', 'com_y', 'pitch_angle', 'angular_momentum'});
    x_values = current_table.time;
    y_values = current_table.angular_momentum;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Velocity [deg/hop cycle percent]")
    title('Velocity on a Solid Surface vs Hop Cycle Percent')
    hold on;
% Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);
hold off

% Specify the directory where you want to save the figures
% saveDir = 'C:\Users\Public\Documents\analysis\figures\COM\DD118A';
% 
% Create the directory if it doesn't exist
% if ~isfolder(saveDir)
%     mkdir(saveDir);
% end
% 
% Store the figures in a cell array
% figures = {figure(1), figure(2), figure(3), figure(4), figure(5), figure(6), figure(7), figure(8), figure(9), figure(10), figure(11), figure(12)};
% 
% Save each figure with its original title
% for i = 1:numel(figures)
%     Set the current figure
%     figure(figures{i});
% 
%     Extract the original title
%     figTitle = get(get(gca, 'title'), 'string');
% 
%     Save the figure with its original title
%     saveas(figures{i}, fullfile(saveDir, [figTitle, '.png'])); % Save in the specified directory
% end
 
% Define the callback function for the line click event
function lineClicked(src, ~)
    cycle_num_str = get(src, 'DisplayName');
    cycle_num = str2double(regexp(cycle_num_str, '\d+', 'match'));
    if ~isnan(cycle_num)
        disp(['Clicked on cycle number: ' num2str(cycle_num)]);
    end
end

