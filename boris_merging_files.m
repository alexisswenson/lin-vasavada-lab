clear
clc

% Specify the folders containing the sand and solid files
folder_path_sand = 'D:\DD121 - InVivo sand\BORIS\Solid\400 hz';
%folder_path_solid = 'D:\DD119_all\BORIS\solid';

% Get a list of files in each folder
file_list_sand = dir(fullfile(folder_path_sand, '*.csv'));
%file_list_solid = dir(fullfile(folder_path_solid, '*.csv'));

% Initialize cell arrays to store data from each folder
data_cell_sand = cell(length(file_list_sand), 1);
%data_cell_solid = cell(length(file_list_solid), 1);

% Read and store data from files in the sand folder
for i = 1:length(file_list_sand)
    file_path = fullfile(folder_path_sand, file_list_sand(i).name);
    data = readtable(file_path);
    data_cell_sand{i} = data;
end

% % Read and store data from files in the solid folder
% for i = 1:length(file_list_solid)
%     file_path = fullfile(folder_path_solid, file_list_solid(i).name);
%     data = readtable(file_path);
%     data_cell_solid{i} = data;
% end

% Vertically concatenate the data from all sand files
concatenated_data_sand = vertcat(data_cell_sand{:});

 % Vertically concatenate the data from all solid files
 % concatenated_data_solid = vertcat(data_cell_solid{:});

% Define the file paths for the concatenated sand and solid files
 output_file_path_sand = 'D:\DD121 - InVivo sand\BORIS\DD121_solid!.csv';
% output_file_path_solid = 'D:\DD119_all\BORIS\DD119_merged_solid.csv';

% % Write the concatenated sand data to a new CSV file
 writetable(concatenated_data_sand, output_file_path_sand);

% Write the concatenated solid data to a new CSV file
% writetable(concatenated_data_solid, output_file_path_solid);
