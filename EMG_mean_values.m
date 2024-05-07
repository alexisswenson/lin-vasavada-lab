function [mean_values] = EMG_mean_values(EMG_cycles)

num_tables = numel(EMG_cycles);

% Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 9; % Assuming you have 9 columns of interest (excluding the time column)
aligned_data_matrix = zeros(num_pseudotime_points, num_columns);

for i = 1:num_tables
    current_data = EMG_cycles{i};
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
mean_values=array2table(mean_values_array, 'VariableNames', {'time', 'SONO_LG', 'Buckle', 'EMG_1_LG', 'EMG_2_LG', 'filtered_SONO', 'filtered_Buckle', 'filtered_EMG_1', 'filtered_EMG_2'});
