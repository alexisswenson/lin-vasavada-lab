clear
clc

% read in kinematics and boris data files
kinematics = readtable("D:\DD121 - InVivo sand\kin and an\DD121_solid!.csv");
boris = readtable("D:\DD121 - InVivo sand\BORIS\DD121_solid!.csv");

% combine the files
merged = horzcat(kinematics, boris);

% Define the file name and path where you want to save the CSV file
file_name = 'merged_DD121_solid4000!.csv';
file_path = 'D:\DD121 - InVivo sand\merged';

% Combine the file name and path
full_file_path = fullfile(file_path, file_name);

% Export the merged table to a CSV file
writetable(merged, full_file_path);