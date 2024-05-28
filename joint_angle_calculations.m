clear
clc

% read in data from correct filepath
 data_nan = readtable("files for MATLAB/merged/DD139_slope_all trials_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD133_tread_slopes_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD133_slopes_merged.csv");

% data_nan = readtable("files for MATLAB/merged/DD130_slopes_01_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD130_slopes_02_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD135_slopes_02_merged.csv");

% data_nan = readtable("files for MATLAB/merged/SCDD03_hard.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD03_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand2.csv");

% remove all NAN cells
data = rmmissing(data_nan);

% indexing
[upperback, hip, knee, ankle, foot, toetip] = krat_table(data);

% defining t 
t=0:height(data)-1;

% calculating joint angles
hipang_ = seg_dist(knee,hip,upperback);
kneeang_ = seg_dist_knee(ankle,knee,hip);
ankleang_ = seg_dist(foot,ankle,knee);
footang_ = seg_dist(toetip,foot,ankle);

% creating table with angle values
angles = [hipang_, kneeang_, ankleang_, footang_];

% plot first 100 frames
figure(1)
plot(t(1:100),hipang_(1:100))
title('Hip, Knee, Ankle, and Foot Angle vs Time')
xlabel('Time')
ylabel('Angle (degrees)')

hold on
plot(t(1:100),kneeang_(1:100))
plot(t(1:100),ankleang_(1:100))
plot(t(1:100),footang_(1:100))
legend('Hip Angle','Knee Angle','Ankle Angle', 'Foot Angle')
hold off

% add landing, steadyhopping, and wedges column to angles
landing_column = table2array(data(:,'Landing'));
angles = [angles, landing_column];

steadyhopping_column =  table2array(data(:,'SteadyHopping'));
angles = [angles, steadyhopping_column];

wedgeback_column = table2array(data(:,'Wedge_WedgeBack_'));
angles = [angles, wedgeback_column];

wedgefront_column = table2array(data(:,'Wedge_WedgeFront_'));
angles = [angles, wedgefront_column];

% Check if the 'wedge' column exists in the original data table
if ismember('Wedge', data.Properties.VariableNames)
    wedge_column = table2array(data(:,'Wedge'));
    angles = [angles, wedge_column];

    % Get the original angles variable names
    original_angle_names = {'angles1', 'angles2', 'angles3', 'angles4', 'angles5', 'angles6', 'angles7', 'angles8', 'angles9'};
    new_angle_names = {'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle', 'landing', 'steadyhopping', 'wedge_back', 'wedge_front', 'wedge'};

    angles_with_time_table = array2table(angles);

    % Rename the angles variable names
    angles_table = renamevars(angles_with_time_table, original_angle_names, new_angle_names);
else
    % Get the original angles variable names
    original_angle_names = {'angles1', 'angles2', 'angles3', 'angles4', 'angles5', 'angles6', 'angles7', 'angles8'};
    new_angle_names = {'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle', 'landing', 'steadyhopping', 'wedge_back', 'wedge_front'};
    
    angles_with_time_table = array2table(angles);

    % Rename the angles variable names
    angles_table = renamevars(angles_with_time_table, original_angle_names, new_angle_names);
end

% Check if the 'wedge' column exists in the original data table
if ismember('Wedge', angles_table.Properties.VariableNames)
    filtered_data = angles_table(angles_table.steadyhopping == 1 & angles_table.wedge_back == 0 & angles_table.wedge_front == 0 & angles_table.wedge == 0, :);
else
    filtered_data = angles_table(angles_table.steadyhopping == 1 & angles_table.wedge_back == 0 & angles_table.wedge_front == 0, :);
end
   
% Find the indices where the column equals 1
indices = find(filtered_data.landing(:, 1) == 1);

% Split the data into cycles
cycles = cell(length(indices)-1, 1);
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = filtered_data(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles{i} = a;
end

% Define column names for indexing
time_column_name = 'time';
hip_angle_column_name = 'hip_angle';
knee_angle_column_name = 'knee_angle';
ankle_angle_column_name = 'ankle_angle';
foot_angle_column_name = 'foot_angle';

% Group the column names in a cell array
column_names = {time_column_name, hip_angle_column_name, knee_angle_column_name, ankle_angle_column_name, foot_angle_column_name};

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

% defining and plotting the joint angles
cycles_array = data_cell;
% num_hop_cycles = 10;
num_hop_cycles = height(cycles_array);

% plot angles across all good hop cycles
figure(2)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', 'current_array5'}, {'time', 'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle'});
    x_values = current_table.time;
    y_values = current_table.hip_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Hip Angle")
    title('Hip Angles vs Hop Cycle Percent')
    hold on;
end
hold off

figure(3)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', 'current_array5'}, {'time', 'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle'});
    x_values = current_table.time;
    y_values = current_table.knee_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Knee Angle")
    title('Knee Angles vs Hop Cycle Percent')
    hold on;
end
hold off

figure(4)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', 'current_array5'}, {'time', 'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle'});
    x_values = current_table.time;
    y_values = current_table.ankle_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Ankle Angle")
    title('Ankle Angles vs Hop Cycle Percent')
    hold on;
end
hold off;

figure(5)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', 'current_array5'}, {'time', 'hip_angle', 'knee_angle', 'ankle_angle', 'foot_angle'});
    x_values = current_table.time;
    y_values = current_table.foot_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Foot Angle")
    title('Foot Angles vs Hop Cycle Percent')
    hold on;
end
hold off;

% Define the column names for each table (struct) within a_data_cell
column_names = {'time', 'Hip_Angle', 'Knee_Angle', 'Ankle_Angle', 'Foot_Angle'};

% Loop through each table (struct) in a_data_cell and set the column names
for i = 1:numel(data_cell)
    data_cell{i} = array2table(data_cell{i}, 'VariableNames', column_names);
end

% Finding the means across each cycle with the corresponding psuedotime

% Step 1: Compute the mean pseudotime for each cycle
pseudotime_mean_per_cycle = cellfun(@(data) mean(data.time), data_cell);

% Step 2: Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 5; % Assuming you have 5 columns of interest (excluding the time column)
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
mean_values=array2table(mean_values_array, 'VariableNames', {'Time', 'Hip_Angle', 'Knee_Angle', 'Ankle_Angle', 'Foot_Angle'});

% Calculating the standard deviation of each column across all cycles at corresponding pseudotime points
std_values = std(mean_values, 0, 1);

% axis_limits = [0, 100, -50, 125];
% axis(axis_limits)

% Plotting the mean joint angles with error bars for +/- 1 standard deviation
figure(6)
errorbar(mean_values.Time, mean_values.Hip_Angle, std_values.Hip_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.Hip_Angle, -std_values.Hip_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

% % Add filled area between the standard deviation lines
% x = [mean_values.Time; flipud(mean_values.Time)];
% y = [mean_values.Hip_Angle + std_values.Hip_Angle; flipud(mean_values.Hip_Angle - std_values.Hip_Angle)];
% fill(x, y, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% hold off;

xlabel("HOP CYCLE PERCENT")
ylabel("Hip Angle")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Hip Angle vs Hop Cycle Percent')

figure(7)
errorbar(mean_values.Time, mean_values.Knee_Angle, std_values.Knee_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.Knee_Angle, -std_values.Knee_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Knee Angle")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Knee Angle vs Hop Cycle Percent')

figure(8)
errorbar(mean_values.Time, mean_values.Ankle_Angle, std_values.Ankle_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.Ankle_Angle, -std_values.Ankle_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Ankle Angle")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Ankle Angle vs Hop Cycle Percent')

figure(9)
errorbar(mean_values.Time, mean_values.Foot_Angle, std_values.Foot_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.Foot_Angle, -std_values.Foot_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Foot Angle")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Foot Angle vs Hop Cycle Percent')