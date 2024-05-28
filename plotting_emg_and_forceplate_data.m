clear
clc

% read in data from correct filepath
% data_nan = readtable("D:\DD118B - invivo sand\merged\merged_DD118B_solid!.csv");       %filtered before downsampling
% data_nan = readtable("D:\DD118B - invivo sand\merged\merged_DD118B_sand!.csv");
 data_nan = readtable('D:\DD118A_all\merged\merged_DD118A_solid!.csv');
% data_nan = readtable('D:\DD118A_all\merged\merged_DD118A_sand!.csv');
% data_nan = readtable('D:\DD119_all\merged\merged_DD119_sand!.csv');
% data_nan = readtable('D:\DD119_all\merged\merged_DD119_solid!.csv');
% data_nan = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_solid!.csv');
% data_nan = readtable('D:\DD121 - InVivo sand\merged\merged_DD121_sand5!.csv'); 

% filter data to get only good hops (=1)
filtered_data = data_nan(data_nan.SteadyHopping == 1, :);

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
    time= linspace(1, 100, numRows)';
    a.time = time;
    cycles_all{i} = a;
end

% Create a cell array to store the new tables
cycles = cell(size(cycles_all));

% Define the relevant columns
selectedColumns = {'time', 'SonoLG', 'Buckle', 'EMG1LG', 'EMG2LG', 'filtered_SONO', 'filtered_Buckle', 'filtered_EMG_1', 'filtered_EMG_2'};

% Iterate through each table in the 'cycles_all' cell array
for i = 1:numel(cycles_all)
    % Select the relevant columns from the current table
    current_table = cycles_all{i};
    cycles{i} = current_table(:, selectedColumns);
end

% gets number of tables in cell array
 num_hop_cycles = numel(cycles);
% num_hop_cycles = 1;

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

hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);


figure(2)
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
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(3)
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
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);

figure(4)
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.EMG1LG;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("EMG 1 LG")
    title("Unfiltered EMG 1 vs Hop Cycle Percent")
    hold on

    % Create a separate line plot for each cycle
    line_handle = plot(x_values, y_values, 'DisplayName', ['Cycle ' num2str(i)]);
end
hold off

% Set the callback function for all the line plots
set(findobj(gca, 'Type', 'Line'), 'ButtonDownFcn', @lineClicked);


figure(5)
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.SonoLG;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("SONO LG (Volts)")
    title("Unfiltered SONO LG vs Hop Cycle Percent")
    hold on
end
hold off

figure(6)
for i = 1:num_hop_cycles
    current_table = cycles{i}; 
    x_values = current_table.time; 
    y_values = current_table.Buckle;
    plot(x_values, y_values);
    xlabel("HOP CYCLE PERCENT")
    ylabel("Buckle (Volts)")
    title("Unfiltered Buckle vs Hop Cycle Percent")
    hold on;
end
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
figure(7)
errorbar(mean_values.time, mean_values.filtered_SONO, std_values.filtered_SONO, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_SONO, -std_values.filtered_SONO, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("SONO LG (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered SONO LG vs Hop Cycle Percent')

figure(8)
errorbar(mean_values.time, mean_values.filtered_Buckle, std_values.filtered_Buckle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_Buckle, -std_values.filtered_Buckle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("Buckle (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered Buckle vs Hop Cycle Percent')

figure(9)
errorbar(mean_values.time, mean_values.filtered_EMG_1, std_values.filtered_EMG_1, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_values.time, mean_values.filtered_EMG_1, -std_values.filtered_EMG_1, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold off;
xlabel("HOP CYCLE PERCENT")
ylabel("EMG 1 LG (Volts)")
legend('Mean', '+1 Std Dev', '-1 Std Dev')
title('Average Filtered EMG 1 vs Hop Cycle Percent')



% % Specify the directory where you want to save the figures
% saveDir = 'C:\Users\Public\Documents\analysis\figures\EMG data\DD118B';
% 
% % Create the directory if it doesn't exist
% if ~isfolder(saveDir)
%     mkdir(saveDir);
% end
% 
% % Store the figures in a cell array
% figures = {figure(1), figure(2), figure(3), figure(4), figure(5), figure(6), figure(7), figure(8), figure(9)};
% 
% % Save each figure with its original title
% for i = 1:numel(figures)
%     % Set the current figure
%     figure(figures{i});
% 
%     % Extract the original title
%     figTitle = get(get(gca, 'title'), 'string');
% 
%     % Save the figure with its original title
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

