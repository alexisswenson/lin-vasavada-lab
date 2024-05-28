clear
clc

% Sand
%data = readtable("C:\Users\Public\Documents\scrape_experimental_txt_files\scrape_experimental_txt_files\no_trigger_offsets_align_data\DD119\DD119_sand_merged_data.csv");
 
% Solid
 data = readtable("files for MATLAB/merged/DD118B/DD118B_event07.csv");
 

[data_updated] = kin_filt(data);

% indexing
[upperback, hip, knee, ankle, foot, toetip, eye] = filt_krat_table(data_updated);

% nameing columns x and y for eye amd hip
eye = array2table(eye, 'VariableNames', {'x_eye', 'y_eye'});
hip = array2table(hip, 'VariableNames', {'x_hip', 'y_hip'});

% Scale factor of 2.5 mm/pixel
scale_factor = 2.5;

% Divide the x and y coordinates of eye by the scale factor
eye.x_eye = eye.x_eye / scale_factor;
eye.y_eye = eye.y_eye / scale_factor;

% Divide the x and y coordinates of hip by the scale factor
hip.x_hip = hip.x_hip / scale_factor;
hip.y_hip = hip.y_hip / scale_factor;

% Calculate the COM coordinates for the filtered data
com_x = (eye.x_eye + hip.x_hip) / 2;
com_y = (eye.y_eye + hip.y_hip) / 2;

% Initialize an array to store the internal angles
pitch_angle = zeros(height(hip), 1);

% Calculate internal angles for each pair of coordinates
for i = 1:height(hip)
    x_eye = eye.x_eye(i);
    y_eye = eye.y_eye(i);
    
    x_hip = hip.x_hip(i);
    y_hip = hip.y_hip(i);
    
    % Calculate the angle in radians using the arctangent function (atan2)
    pitch_angle(i) = atan2(y_hip - y_eye, x_hip - x_eye);
end

% Convert the angles to degrees if needed
angles_array = rad2deg(pitch_angle);

% Create a table for the COM coordinates
com_table = table(com_x, com_y, angles_array, 'VariableNames', {'com_x', 'com_y', 'pitch_angle'});

% Negate the 'com_y' and 'pitch_angle' columns
com_table.pitch_angle = -com_table.pitch_angle;
com_table.com_y = -com_table.com_y;

% Convert the 'Landing' and 'SteadyHopping' variables to tables
landing_column = table(data_updated.Landing, 'VariableNames', {'Landing'});
steadyhopping_column = table(data_updated.SteadyHopping, 'VariableNames', {'SteadyHopping'});

% Add the 'Landing' and 'SteadyHopping' columns to com_table
com_table = [com_table, landing_column, steadyhopping_column];

% Find the indices where the column equals 1
indices = find(com_table.Landing(:, 1) == 1);

% Split the data into cycles
cycles = cell(length(indices)-1, 1);
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = com_table(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles{i} = a;
end

% Define column names for indexing
time_column_name = 'time';
com_y_column_name = 'com_y';
pitch_angle_column_name = 'pitch_angle';

% Group the column names in a cell array
column_names = {time_column_name, com_y_column_name, pitch_angle_column_name};

% Iterate over columns to get angles and pseudotime
% Grouped each cycle with this data in cells
cell_array = cycles;
num_tables = numel(cell_array);
num_columns = numel(column_names);
data_cell = cell(num_tables, 1); 

for i = 1:num_tables
    current_table = cell_array{i};
    columns_data = cell(1, num_columns);
    for j = 1:num_columns
        a_column_data = current_table.(column_names{j});
        columns_data{j} = a_column_data;
    end
    data_cell{i} = [columns_data{:}];
end

% defining and plotting the COM
cycles_array = data_cell;
num_hop_cycles = height(cycles_array);

% Iterate over each cell in cycles_array
for i = 1:numel(cycles_array)
    current_data = cycles_array{i};  % Get the current cell

    % Choose the column to subtract the mean from (e.g., column 1)
    column_to_subtract_mean_from = 2;

    % Calculate the mean of the selected column
    mean_value = mean(current_data(:, column_to_subtract_mean_from));

    % Subtract the mean from every point in the selected column
    current_data(:, column_to_subtract_mean_from) = current_data(:, column_to_subtract_mean_from) - mean_value;

    % Update the current cell in cycles_array
    cycles_array{i} = current_data;
end


% Define the column names for each table (struct) within a_data_cell
column_names = {'time', 'com', 'pitch_angle'};

% Loop through each table (struct) in a_data_cell and set the column names
for i = 1:numel(cycles_array)
    data_cell{i} = array2table(cycles_array{i}, 'VariableNames', column_names);
end

% Finding the means across each cycle with the corresponding psuedotime
% Step 1: Compute the mean pseudotime for each cycle
pseudotime_mean_per_cycle = cellfun(@(data) mean(data), data_cell, 'UniformOutput', false);

% Step 2: Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 3; % Assuming you have 2 columns of interest (excluding the time column)
aligned_data_matrix = zeros(num_pseudotime_points, num_columns);

for i = 1:num_tables
    current_data = data_cell{i};
    time_vector = current_data{:, 'time'};
    columns_data = table2array(current_data(:, 1:end));

    % Resample the data to align with a common pseudotime range
    pseudotime_range = linspace(time_vector(1), time_vector(end), num_pseudotime_points);
    aligned_data = interp1(time_vector, columns_data, pseudotime_range, 'spline');

    % Add the aligned data to the combined matrix
    aligned_data_matrix = aligned_data_matrix + aligned_data;
end

% Step 3: Compute the mean of each column over all cycles at corresponding pseudotime points
mean_values_array = aligned_data_matrix / num_tables;
mean_values=array2table(mean_values_array, 'VariableNames', {'Time', 'COM', 'pitch_angle'});

% Calculating the standard deviation of each column across all cycles at corresponding pseudotime points
std_values = std(mean_values, 0, 1);

% Calculate the average displacement across all cycles
average_displacement = mean(mean_values.COM);
average_angle = mean(mean_values.pitch_angle);

figure(1)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'com_y', 'pitch_angle'});
    x_values = current_table.time;
    % current_table.com_y =
    y_values = current_table.com_y;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Center of Mass Displacement [mm]")
    title('Center of Mass Displacement on a Sand Surface vs Hop Cycle Percent')
    hold on;
end
hold off


figure(2)
errorbar(mean_values.Time, mean_values.COM, std_values.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.COM, -std_values.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Center of Mass Displacement [mm]")
legend('Mean', 'Std Dev')
title('Average Center of Mass Displacement on a Sand Surface vs Hop Cycle Percent')


figure(3)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'com_y', 'pitch_angle'});
    x_values = current_table.time;
    y_values = current_table.pitch_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Pitch Angle [deg]")
    title('Pitch Angle on a Sand Surface vs Hop Cycle Percent')
    hold on;
end
hold off

figure(4)
errorbar(mean_values.Time, mean_values.pitch_angle, std_values.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.pitch_angle, -std_values.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Pitch Angle [deg]")
legend('Mean', 'Std Dev')
title('Average Pitch Angle on a Sand Surface vs Hop Cycle Percent')



% figure(1)
% for i = 1:num_hop_cycles
%     current_array = cycles_array{i};
%     current_table_wrong_names=array2table(current_array);
%     current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'com_y', 'pitch_angle'});
%     x_values = current_table.time;
%     y_values = current_table.com_y;
%     plot(x_values, y_values);
%     xlabel("HOP CYCLE PERCENT")
%     ylabel("Center of Mass Displacement [mm]")
%     title('Center of Mass Displacement on a Solid Surface vs Hop Cycle Percent')
%     hold on;
% end
% hold off
% 
% 
% figure(2)
% errorbar(mean_values.Time, mean_values.COM, std_values.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
% hold on;
% errorbar(mean_values.Time, mean_values.COM, -std_values.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
% hold off;
% xlabel("HOP CYCLE PERCENT")
% ylabel("Center of Mass Displacement [mm]")
% legend('Mean', 'Std Dev')
% title('Average Center of Mass Displacement on a Solid Surface vs Hop Cycle Percent')
% 
% 
% figure(3)
% for i = 1:num_hop_cycles
%     current_array = cycles_array{i};
%     current_table_wrong_names=array2table(current_array);
%     current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'com_y', 'pitch_angle'});
%     x_values = current_table.time;
%     y_values = current_table.pitch_angle;
%     plot(x_values, y_values);
%     xlabel("HOP CYCLE PERCENT")
%     ylabel("Pitch Angle [deg]")
%     title('Pitch Angle on a Solid Surface vs Hop Cycle Percent')
%     hold on;
% end
% hold off
% 
% figure(4)
% errorbar(mean_values.Time, mean_values.pitch_angle, std_values.COM, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
% hold on;
% errorbar(mean_values.Time, mean_values.pitch_angle, -std_values.COM, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
% hold off;
% xlabel("HOP CYCLE PERCENT")
% ylabel("Pitch Angle [deg]")
% legend('Mean', 'Std Dev')
% title('Average Pitch Angle on a Solid Surface vs Hop Cycle Percent')
