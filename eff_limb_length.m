clc
clear

data_nan = readtable("files for MATLAB/merged/DD_old/DD139_slope_all trials_merged.csv");

% remove all NAN cells
data = rmmissing(data_nan);

[data_updated] = kin_filt(data);

% indexing
[upperback, hip, knee, ankle, foot, toetip] = filt_krat_table(data_updated);

% nameing columns x and y fpr toe amd hip
toetip = array2table(toetip, 'VariableNames', {'x_toe', 'y_toe'});
hip = array2table(hip, 'VariableNames', {'x_hip', 'y_hip'});

% differences between x and y coords
dx = hip.x_hip - toetip.x_toe;
dy = hip.y_hip - toetip.y_toe;

% length of line segment between hip and toetip
length_table = sqrt(dx.^2 + dy.^2);

% Initialize an array to store the internal angles
internalAngles = zeros(height(hip), 1);

% Calculate internal angles for each pair of coordinates
for i = 1:height(hip)
    x_toe = toetip.x_toe(i);
    y_toe = toetip.y_toe(i);
    
    x_hip = hip.x_hip(i);
    y_hip = hip.y_hip(i);
    
    % Calculate the angle in radians using the arctangent function (atan2)
    internalAngles(i) = atan2(y_hip - y_toe, x_hip - x_toe);
end

% Convert the angles to degrees if needed
angles_array = rad2deg(internalAngles);
angles_column_name = 'angles';
angles = array2table(angles_array, 'VariableNames', {angles_column_name});

% add angles, landing, steadyhopping, and wedges column to length_table
angles_column = table2array(angles(:,'angles'));
length_table = [length_table, angles_column];

landing_column = table2array(data_updated(:,'Landing'));
length_table = [length_table, landing_column];

steadyhopping_column =  table2array(data_updated(:,'SteadyHopping'));
length_table = [length_table, steadyhopping_column];

wedgeback_column = table2array(data_updated(:,'Wedge_WedgeBack_'));
length_table = [length_table, wedgeback_column];

wedgefront_column = table2array(data_updated(:,'Wedge_WedgeFront_'));
length_table = [length_table, wedgefront_column];

% Check if the 'wedge' column exists in the original data table
if ismember('Wedge', data_updated.Properties.VariableNames)
    wedge_column = table2array(data_updated(:,'Wedge'));
    length_table = [angles, wedge_column];

    % Get the original variable names
    original_names = {'length_table1', 'length_table2', 'length_table3', 'length_table4', 'length_table5', 'length_table6', 'length_table7'};
    new_names = {'length', 'angles', 'landing', 'steadyhopping', 'wedge_back', 'wedge_front', 'wedge'};

    length_table = array2table(length_table);

    % Rename the variable names
    length_table = renamevars(length_table, original_names, new_names);
else
    % Get the original variable names
    original_names = {'length_table1', 'length_table2', 'length_table3', 'length_table4', 'length_table5', 'length_table6'};
    new_names = {'length', 'angles', 'landing', 'steadyhopping', 'wedge_back', 'wedge_front'};
    
    length_table = array2table(length_table);

    % Rename the variable names
    length_table = renamevars(length_table, original_names, new_names);
end

% Check if the 'wedge' column exists in the original data table
if ismember('Wedge', length_table.Properties.VariableNames)
    filtered_length_table = length_table(length_table.steadyhopping == 1 & length_table.wedge_back == 0 & length_table.wedge_front == 0 & length_table.wedge == 0, :);
else
    filtered_length_table = length_table(length_table.steadyhopping == 1 & length_table.wedge_back == 0 & length_table.wedge_front == 0, :);
end

% Find the indices where the column equals 1
indices = find(filtered_length_table.landing(:, 1) == 1);

% Split the data into cycles
cycles = cell(length(indices)-1, 1);
for i = 1:length(indices)-1
    startIdx = indices(i);
    endIdx = indices(i+1)-1;
    a = filtered_length_table(startIdx:endIdx, :);

    % make pseudo-time
    numRows = size(a, 1);
    time = linspace(1, 100, numRows)';
    a.time = time;
    cycles{i} = a;
end

% Define column names for indexing
time_column_name = 'time';
length_column_name = 'length';
angles_column_name = 'angles';

% Group the column names in a cell array
column_names = {time_column_name, length_column_name, angles_column_name};

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

% defining and plotting the effective limb length
cycles_array = data_cell;
num_hop_cycles = height(cycles_array);

figure(1)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'length', 'angles'});
    x_values = current_table.time;
    y_values = current_table.length;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Effective Limb Length")
    title('Effective Limb Length vs Hop Cycle Percent')
    hold on;
end
hold off

figure(2)
for i = 1:num_hop_cycles
    current_array = cycles_array{i};
    current_table_wrong_names=array2table(current_array);
    current_table = renamevars(current_table_wrong_names, {'current_array1','current_array2','current_array3'}, {'time', 'length', 'angles'});
    x_values = current_table.time;
    y_values = current_table.angles;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Effective Limb Length Angle")
    title('Effective Limb Length Angle vs Hop Cycle Percent')
    hold on;
end
hold off
