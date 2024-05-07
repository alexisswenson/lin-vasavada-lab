function [back_w_cycles] = back_wedge_cycles(data_nan)

% remove all NAN cells
data = rmmissing(data_nan);

[data_updated] = kin_filt(data);

% indexing
[upperback, hip, knee, ankle, foot, toetip] = filt_krat_table(data_updated);

% calculating joint angles
hipang_= seg_dist(knee,hip,upperback);
kneeang_= seg_dist_knee(ankle,knee,hip);
ankleang_ = seg_dist(foot,ankle,knee);
footang_ = seg_dist(toetip,foot,ankle);

% creating table with angle values
angles = [hipang_, kneeang_, ankleang_, footang_];

% add landing, steadyhopping, and wedges column to angles
landing_column = table2array(data_updated(:,'Landing'));
angles = [angles, landing_column];

steadyhopping_column =  table2array(data_updated(:,'SteadyHopping'));
angles = [angles, steadyhopping_column];

wedgeback_column = table2array(data_updated(:,'Wedge_WedgeBack_'));
angles = [angles, wedgeback_column];

wedgefront_column = table2array(data_updated(:,'Wedge_WedgeFront_'));
angles = [angles, wedgefront_column];

% Check if the 'wedge' column exists in the original data table
if ismember('Wedge', data_updated.Properties.VariableNames)
    wedge_column = table2array(data_updated(:,'Wedge'));
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

if ismember('Wedge', angles_table.Properties.VariableNames)
filtered_angles_table = angles_table(angles_table.steadyhopping == 1 & (angles_table.wedge_back == 1 | angles_table.wedge == 1), :);
else
filtered_angles_table = angles_table(angles_table.steadyhopping == 1 & angles_table.wedge_back == 1, :);
end
   
% Find the indices where the column equals 1
indices = find(filtered_angles_table.landing(:, 1) == 1);

% Split the data into cycles
cycles = cell(length(indices)-1, 1);
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = filtered_angles_table(startIdx:endIdx, :);

    % make pseudo-time
    numRows = size(a, 1);
    time = linspace(1, 100, numRows)';
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

% Iterate over cycles to get angles and pseudotime
% Group each cycle data in separate tables and store them in a cell array
cell_array = cycles;
num_tables = numel(cell_array);
num_columns = numel(column_names);
back_w_cycles = cell(num_tables, 1);

for i = 1:num_tables
    current_table = cell_array{i};
    columns_data = cell(1, num_columns);
    for j = 1:num_columns
        a_column_data = current_table.(column_names{j});
        columns_data{j} = a_column_data;
    end
    back_w_cycles{i} = table(columns_data{:}, 'VariableNames', column_names);
end