clear
clc

%DD139_nan = readtable("files for MATLAB/merged/DD139_slope_all trials_merged.csv");
DD139_nan = readtable("files for MATLAB/merged/DD_old/DD133_tread_slopes_merged.csv");
%DD139_nan = readtable("files for MATLAB/merged/DD133_slopes_merged.csv");

%DD139_nan = readtable("files for MATLAB/merged/DD133_tread_slopes_merged.csv");

DD139=rmmissing(DD139_nan);

krat_table(DD139)

%indexing
[upperback, hip, knee, ankle, foot, toetip] = krat_table(DD139);

%defining t 
t=0:height(DD139)-1;

%calculating joint angles
hipang_=seg_dist(knee,hip,upperback);
kneeang_=seg_dist_knee(ankle,knee,hip);
ankleang_=seg_dist(toetip,ankle,knee);

alexis_angles = [hipang_, kneeang_, ankleang_];

% plot first 100 frames
figure(1)
plot(t(1:100),hipang_(1:100))
title('Hip, Knee, and Ankle Angle vs Time')
xlabel('Time')
ylabel('Angle (degrees)')

hold on
plot(t(1:100),kneeang_(1:100))
plot(t(1:100),ankleang_(1:100))
legend('Hip Angle','Knee Angle','Ankle Angle')
hold off

%Nick's angles
A_hip=DD139.A_hip;
A_knee=DD139.A_knee;
A_ankle=DD139.A_ankle;

%mine and nicks angles
figure(2)
plot(t(1:100),hipang_(1:100))
title('Hip, Knee, and Ankle Angle vs Time')
xlabel('Time')
ylabel('Angle (degrees)')

hold on
plot(t(1:100),kneeang_(1:100))
plot(t(1:100),ankleang_(1:100))
plot(t(1:100),A_hip(1:100))
plot(t(1:100),A_knee(1:100))
plot(t(1:100),A_ankle(1:100))
legend('Hip Angle','Knee Angle','Ankle Angle','Nick Hip Angle','Nick Knee Angle','Nick Ankle Angle')
hold off

% Get the names of all columns in the table
DD139_column_names = DD139.Properties.VariableNames;

% Define the range of columns to index through
start_column_name = 't1y';
end_column_name = 'nosey_v';

% Find the indices of the start and end columns
start_index = find(strcmp(DD139_column_names, start_column_name));
end_index = find(strcmp(DD139_column_names, end_column_name));

% Extract the selected columns and create a new table
DD139_noEMG = DD139(:, start_index:end_index);

%filter data to get only good hops (=1) and to get jumps not affected by wedge (=0)
filteredDD139 = DD139_noEMG(DD139_noEMG.SteadyHopping == 1 & DD139_noEMG.Wedge_WedgeBack_ == 0 & DD139_noEMG.Wedge_WedgeFront_ == 0, :);
% filteredDD139 = DD139_noEMG(DD139_noEMG.Wedge == 0,:)

% Find the indices where the column equals 1
indices = find(filteredDD139.Landing(:, 1) == 1);

% Split the data into cycles
cycles = cell(length(indices)-1, 1);
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = filteredDD139(startIdx:endIdx, :);
    
    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles{i} = a;
end

%iterate over columns to get nicks angles and pseudotime
%grouped each cycle with this data in cells
cell_array = cycles;
column_indices = [35, 49, 50, 51];

num_hop_cycles = numel(cell_array);
num_columns = numel(column_indices);
data_cell = cell(num_hop_cycles, 1); 

for i = 1:num_hop_cycles
    current_table = cell_array{i};
    columns_data = cell(1, num_columns);
    for j = 1:num_columns
        column_data = current_table{:, column_indices(j)};
        columns_data{j} = column_data;
    end
    data_cell{i} = cat(2, columns_data{:});
end


%add landing, steadyhopping, and wedges column to alexis_angles_with_time
landing_column = table2array(DD139(:,'Landing'));
alexis_angles = [alexis_angles, landing_column];

steadyhopping_column =  table2array(DD139(:,'SteadyHopping'));
alexis_angles = [alexis_angles, steadyhopping_column];

wedgeback_column = table2array(DD139(:,'Wedge_WedgeBack_'));
alexis_angles = [alexis_angles, wedgeback_column];

wedgefront_column = table2array(DD139(:,'Wedge_WedgeFront_'));
alexis_angles = [alexis_angles, wedgefront_column];

%check for a wedge columns
hasWedge = ismember('Wedge', DD139_noEMG.Properties.VariableNames);
hasWedgeFront = ismember('Wedge_WedgeFront', DD139_noEMG.Properties.VariableNames);
hasWedgeBack = ismember('Wedge_WedgeBack', DD139_noEMG.Properties.VariableNames);

% Filter the data based on the conditions, including the presence or absence of wedge data
if hasWedge && hasWedgeBack && hasWedgeFront
    wedge_column = table2array(DD139(:,'Wedge'));
    alexis_angles = [alexis_angles, wedge_column];
    filteredDD139 = DD139_noEMG(DD139_noEMG.SteadyHopping == 1 & DD139_noEMG.Wedge_WedgeBack_ == 0 & DD139_noEMG.Wedge_WedgeFront_ == 0, :);
    alexis_angles_with_time_table = array2table(alexis_angles);
    alexis_table = renamevars(alexis_angles_with_time_table, {'alexis_angles1', 'alexis_angles2', 'alexis_angles3','alexis_angles4', 'alexis_angles5','alexis_angles6','alexis_angles7','alexis_angles8'}, {'hip_angle', 'knee_angle', 'ankle_angle', 'landing', 'steadyhopping', 'wedge', 'wedge_back', 'wedge_front'});
    %filter data to get only good hops (=1) and to get jumps not affected by wedge (=0)
    filtered_alexis_table = alexis_table(alexis_table.steadyhopping == 1 & alexis_table.wedge == 0 & alexis_table.wedge_back == 0 & alexis_table.wedge_front == 0, :);
    
else
   alexis_angles_with_time_table = array2table(alexis_angles);
   alexis_table = renamevars(alexis_angles_with_time_table, {'alexis_angles1', 'alexis_angles2', 'alexis_angles3','alexis_angles4', 'alexis_angles5','alexis_angles6','alexis_angles7'}, {'hip_angle', 'knee_angle', 'ankle_angle', 'landing', 'steadyhopping', 'wedge_back', 'wedge_front'}); 
   %filter data to get only good hops (=1) and to get jumps not affected by wedge (=0)
   filtered_alexis_table = alexis_table(alexis_table.steadyhopping == 1 & alexis_table.wedge_back == 0 & alexis_table.wedge_front == 0, :);
end

% Find the indices where the column equals 1
indices_a = find(filtered_alexis_table.landing(:, 1) == 1);

% Split the data into cycles
alexis_cycles = cell(length(indices_a)-1, 1);
for i = 1:length(indices_a)-1
    startIdx = indices_a(i);
    endIdx = indices_a(i+1)-1;
    a = filtered_alexis_table(startIdx:endIdx, :);

    % make pseudo-time
    numRows = size(a, 1);
    time= linspace(1, 100, numRows)';
    a.time = time;
    alexis_cycles{i} = a;
end

% Define column names for indexing
time_column_name = 'time';
hip_angle_column_name = 'hip_angle';
knee_angle_column_name = 'knee_angle';
ankle_angle_column_name = 'ankle_angle';

% Group the column names in a cell array
a_column_names = {time_column_name, hip_angle_column_name, knee_angle_column_name, ankle_angle_column_name};

% Iterate over columns to get Alexis' angles and pseudotime
% Grouped each cycle with this data in cells
a_cell_array = alexis_cycles;
a_num_tables = numel(a_cell_array);
a_num_columns = numel(a_column_names);
a_data_cell = cell(a_num_tables, 1); 

for i = 1:a_num_tables
    a_current_table = a_cell_array{i};
    a_columns_data = cell(1, a_num_columns);
    for j = 1:a_num_columns
        a_column_data = a_current_table.(a_column_names{j});
        a_columns_data{j} = a_column_data;
    end
    a_data_cell{i} = [a_columns_data{:}];
end

%defining and plotting the joint angles (alexis)
alexis_array = a_data_cell;
% num_hop_cycles = 10;
num_hop_cycles = height(alexis_array);

figure(3)
for i = 1:num_hop_cycles
    current_array = alexis_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', }, {'time', 'hip_angle', 'knee_angle', 'ankle_angle'});
    x_values = current_table.time;
    y_values = current_table.hip_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Hip Angle")
    hold on;
end

figure(4)
for i = 1:num_hop_cycles
    current_array = alexis_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', }, {'time', 'hip_angle', 'knee_angle', 'ankle_angle'});
    x_values = current_table.time;
    y_values = current_table.knee_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Knee Angle")
    hold on;
end

figure(5)
for i = 1:num_hop_cycles
    current_array = alexis_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3','current_array4', }, {'time', 'hip_angle', 'knee_angle', 'ankle_angle'});
    x_values = current_table.time;
    y_values = current_table.ankle_angle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Ankle Angle")
    hold on;
end

hold off;

% Define the column names for each table (struct) within a_data_cell
column_names = {'time', 'Hip_Angle', 'Knee_Angle', 'Ankle_Angle'};

% Loop through each table (struct) in a_data_cell and set the column names
for i = 1:numel(a_data_cell)
    a_data_cell{i} = array2table(a_data_cell{i}, 'VariableNames', column_names);
end

%Finding the means across each cycle with the corresponding psuedotime

% Step 1: Compute the mean pseudotime for each cycle
pseudotime_mean_per_cycle = cellfun(@(data) mean(data.time), a_data_cell);

% Step 2: Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 4; % Assuming you have 4 columns of interest (excluding the time column)
aligned_data_matrix = zeros(num_pseudotime_points, num_columns);

for i = 1:a_num_tables
    current_data = a_data_cell{i};
    time_vector = current_data{:, 'time'};
    columns_data = table2array(current_data(:, 1:end));

    % Resample the data to align with a common pseudotime range
    pseudotime_range = linspace(time_vector(1), time_vector(end), num_pseudotime_points);
    aligned_data = interp1(time_vector, columns_data, pseudotime_range, 'spline');

    % Add the aligned data to the combined matrix
    aligned_data_matrix = aligned_data_matrix + aligned_data;
end

% Step 3: Compute the mean of each column over all cycles at corresponding pseudotime points
mean_values_array = aligned_data_matrix / a_num_tables;
mean_values=array2table(mean_values_array, 'VariableNames', {'Time', 'Hip_Angle', 'Knee_Angle', 'Ankle_Angle'});


% Calculating standard deviation and plotting it with the means

% Calculating the standard deviation of each column across all cycles at corresponding pseudotime points
std_values = std(mean_values, 0, 1);

% Plotting the mean joint angles (alexis) with error bars for +/- 1 standard deviation
figure(6)
errorbar(mean_values.Time, mean_values.Hip_Angle, std_values.Hip_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.Time, mean_values.Hip_Angle, -std_values.Hip_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
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