function [mean_values] = mean_wedge_cycles(w_cycles)

num_tables = numel(w_cycles);

% Use the mean pseudotime as a reference to align the data
num_pseudotime_points = 100;
num_columns = 5; % Assuming you have 5 columns of interest (excluding the time column)
aligned_data_matrix = zeros(num_pseudotime_points, num_columns);

for i = 1:num_tables
    current_data = w_cycles{i};
    time_vector = current_data{:, 'time'};
    columns_data = table2array(current_data(:, 1:end));

    % Resample the data to align with a common pseudotime range
    pseudotime_range = linspace(time_vector(1), time_vector(end), num_pseudotime_points);
    aligned_data = interp1(time_vector, columns_data, pseudotime_range, 'spline');

    % Add the aligned data to the combined matrix
    aligned_data_matrix = aligned_data_matrix + aligned_data;
end

% Compute the mean of each column over all cycles at corresponding pseudotime points
mean_values_array = aligned_data_matrix / num_tables;
mean_values=array2table(mean_values_array, 'VariableNames', {'Time', 'Hip_Angle', 'Knee_Angle', 'Ankle_Angle', 'Foot_Angle'});

end