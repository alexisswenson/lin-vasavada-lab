clc
clear

%%%%%%%%%%%%%%%%%%%%%%%%%% Input parameters (All you need to change)
animal_number = 'SCDD15';
sand_or_solid = 'sand';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters
input_boris_folder = sprintf('D:\\%s\\BORIS\\%s', animal_number, sand_or_solid);
output_boris_folder = sprintf('D:\\%s\\BORIS\\%s\\400 hz', animal_number, sand_or_solid);
input_data_folder = sprintf('D:\\%s\\kinematics\\1.75 m per s\\%s', animal_number, sand_or_solid);
output_data_folder = sprintf('D:\\%s\\kinematics\\1.75 m per s\\%s\\400 hz', animal_number, sand_or_solid);

% Resample and interpolate BORIS data
boris_files = dir(fullfile(input_boris_folder, '*.xlsx'));
for i = 1:numel(boris_files)
    % Read BORIS data
    boris_data = readtable(fullfile(input_boris_folder, boris_files(i).name));

    % Define the current and desired sampling rates
    current_sampling_rate = 300;  % Hz
    desired_sampling_rate = 400;  % Hz

    % Calculate the resampling factor
    resampling_factor = desired_sampling_rate / current_sampling_rate;

    % Create a new time vector with the desired number of samples
    interpolated_time = linspace(0, (height(boris_data) - 1) / current_sampling_rate, round(resampling_factor * height(boris_data)));

    % Preallocate interpolated data matrix
    interpolated_data = zeros(length(interpolated_time), size(boris_data, 2) - 1);  % Exclude the 'Time' column

    % Resample each column of the data matrix (excluding the 'Time' column)
    for col = 2:size(boris_data, 2)  % Iterate over columns, starting from 2 to exclude 'Time'
        % Interpolate the column using linear interpolation
        interpolated_column = interp1((0:height(boris_data) - 1) / current_sampling_rate, boris_data{:, col}, interpolated_time, 'linear');

        % Round the interpolated values to 0 or 1
        interpolated_column(interpolated_column > 0.62) = 1;
        interpolated_column(interpolated_column < 0.75) = 0;

        % Assign the interpolated column to the resampled data matrix
        interpolated_data(:, col - 1) = interpolated_column;
    end

    % Create a new table with the resampled data
    resampled_boris_data = array2table(interpolated_data, 'VariableNames', boris_data.Properties.VariableNames(2:end));

    % Add a new column for the resampled time at the beginning
    resampled_boris_data = addvars(resampled_boris_data, interpolated_time(:), 'Before', 1, 'NewVariableNames', 'Time');

    % Save the resampled BORIS data to the output folder
    output_filename = fullfile(output_boris_folder, strcat('resampled_', boris_files(i).name));
    writetable(resampled_boris_data, output_filename);
end

% Resample and interpolate other data
data_files = dir(fullfile(input_data_folder, '*.csv'));
for i = 1:numel(data_files)
    % Read data
    data = readtable(fullfile(input_data_folder, data_files(i).name));

    % Define the current and desired sampling rates
    current_sampling_rate = 300;  % Hz
    desired_sampling_rate = 400;  % Hz

    % Calculate the resampling factor
    resampling_factor = desired_sampling_rate / current_sampling_rate;

    % Create a new time vector with the desired number of samples
    interpolated_time = linspace(0, (height(data) - 1) / current_sampling_rate, round(resampling_factor * height(data)));

    % Preallocate interpolated data matrix
    interpolated_data = zeros(length(interpolated_time), size(data, 2));  % Preallocate interpolated data matrix

    % Resample each column of the data matrix
    for col = 1:size(data, 2)  % Iterate over columns
        % Interpolate the column using linear interpolation
        interpolated_data(:, col) = interp1((0:height(data) - 1) / current_sampling_rate, data{:, col}, interpolated_time, 'linear');
    end

    % Create a new table with the resampled data
    resampled_data = array2table(interpolated_data, 'VariableNames', data.Properties.VariableNames);

    % Add a new column for the resampled time at the beginning
    resampled_data = addvars(resampled_data, interpolated_time(:), 'Before', 1, 'NewVariableNames', 'Time');

    % Save the resampled data to the output folder
    output_filename = fullfile(output_data_folder, strcat('resampled_', data_files(i).name));
    writetable(resampled_data, output_filename);
end

disp('Merging and saving completed successfully.');
