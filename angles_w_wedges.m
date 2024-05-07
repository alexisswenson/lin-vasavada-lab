clear
clc

%has steadyhopping with wedges
DD139 = readtable("files for MATLAB/merged/DD_old/DD139_slope_all trials_merged.csv");
[w_cycles_DD139] = wedge_cycles(DD139);

DD133_tread = readtable("files for MATLAB/merged/DD_old/DD133_tread_slopes_merged.csv");
[w_cycles_DD133_tread] = wedge_cycles(DD133_tread);

DD133 = readtable("files for MATLAB/merged/DD_old/DD133_slopes_merged.csv");
[w_cycles_DD133] = wedge_cycles(DD133);

DD135_slopes_02 = readtable("files for MATLAB/merged/DD_old/DD135_slopes_02_merged.csv");
[w_cycles_DD135_slopes_02] = wedge_cycles(DD135_slopes_02);


%does not have steady hopping with wedges
DD130_slopes_01 = readtable("files for MATLAB/merged/DD_old/DD130_slopes_01_merged.csv");
[w_cycles_DD130_slopes_01] = wedge_cycles(DD130_slopes_01);

DD130_slopes_02 = readtable("files for MATLAB/merged/DD_old/DD130_slopes_02_merged.csv");
[w_cycles_DD130_slopes_02] = wedge_cycles(DD130_slopes_01);

SCDD03_hard = readtable("files for MATLAB/merged/SCD/SCDD03_hard.csv");
[w_cycles_SCDD03_hard] = wedge_cycles(SCDD03_hard);

SCDD03_sand = readtable("files for MATLAB/merged/SCD/SCDD03_sand.csv");
[w_cycles_SCDD03_sand] = wedge_cycles(SCDD03_sand);

SCDD10_sand = readtable("files for MATLAB/merged/SCD/SCDD10_sand.csv");
[w_cycles_SCDD10_sand] = wedge_cycles(SCDD10_sand);

SCDD10_sand2 = readtable("files for MATLAB/merged/SCD/SCDD10_sand2.csv");
[w_cycles_SCDD10_sand2] = wedge_cycles(SCDD10_sand2);


% Calculate the mean and std of the cycles
[mean_DD139] = mean_wedge_cycles(w_cycles_DD139);
std_DD139 = std(mean_DD139, 0, 1);

[mean_DD133_tread] = mean_wedge_cycles(w_cycles_DD133_tread);
std_DD133_tread = std(mean_DD133_tread, 0, 1);

[mean_DD133] = mean_wedge_cycles(w_cycles_DD133);
std_DD133 = std(mean_DD133, 0, 1);

[mean_DD135_slopes_02] = mean_wedge_cycles(w_cycles_DD135_slopes_02);
std_DD135_slopes_02 = std(mean_DD135_slopes_02, 0, 1);

% Create a cell array to store the names of each dataset
dataset_names = {'DD139', 'DD133 tread', 'DD133', 'DD135 slopes 02'};

% Create a cell array to store the mean and std values for each dataset
mean_values_all = {mean_DD139, mean_DD133_tread, mean_DD133, mean_DD135_slopes_02};
std_values_all = {std_DD139, std_DD133_tread, std_DD133, std_DD135_slopes_02};

% Define colors and markers for each dataset
colors = {'b', 'r', 'g', 'm'};
markers = {'^', '+', '*', 'x'};

% Iterate through each dataset and plot the mean with error bars (std)
figure(1)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.Time;
y = mean_values.Hip_Angle;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("Hip Angle")
legend(dataset_names, 'Location', 'best')
title('Average Hip Angle vs Hop Cycle Percent')
hold off

figure(2)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.Time;
y = mean_values.Knee_Angle;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("Knee Angle")
legend(dataset_names, 'Location', 'best')
title('Average Knee Angle vs Hop Cycle Percent')
hold off

figure(3)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.Time;
y = mean_values.Ankle_Angle;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("Ankle Angle")
legend(dataset_names, 'Location', 'best')
title('Average Ankle Angle vs Hop Cycle Percent')
hold off


% figure()
% for i = 1:numel(dataset_names)
%     hold on;
%     mean_values = mean_values_all{i};
%     std_values = std_values_all{i};
% 
%     errorbar(mean_values.Time, mean_values.Knee_Angle, std_values.Knee_Angle,  markers{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', 5);
%     errorbar(mean_values.Time, mean_values.Knee_Angle, -std_values.Knee_Angle,  markers{i}, 'MarkerFaceColor', colors{i}, 'MarkerSize', 5);
%     xlabel("HOP CYCLE PERCENT")
%     ylabel("Knee Angle") 
%     legend('Mean', '+1 Std Dev', '-1 Std Dev')
%     title('Average Knee Angle vs Hop Cycle Percent')  
%     hold off
% end


% Now onto EMG data with steadyhopping and wedge cycles
[mean_EMG_DD139] = EMG_mean_values(DD139);
std_EMG_DD139 = std(mean_EMG_DD139, 0, 1);

[mean_EMG_DD133_tread] = EMG_mean_values(DD133_tread);
std_EMG_DD133_tread = std(mean_EMG_DD133_tread, 0, 1);

[mean_EMG_DD133] = EMG_mean_values(DD133);
std_EMG_DD133 = std(mean_EMG_DD133, 0, 1);

% Create a cell array to store the names of each dataset
dataset_names = {'DD139', 'DD133 tread', 'DD133'};

% Create a cell array to store the mean and std values for each dataset
mean_values_all = {mean_EMG_DD139, mean_EMG_DD133_tread, mean_EMG_DD133};
std_values_all = {std_EMG_DD139, std_EMG_DD133_tread, std_EMG_DD133};

% Define colors and markers for each dataset
colors = {'b', 'r', 'g', 'm'};
markers = {'^', '+', '*', 'x'};

% Iterate through each dataset and plot the mean with error bars (std)
figure(4)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.time;
y = mean_values.filtered_SONO;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("SONO LG")
legend(dataset_names, 'Location', 'best')
title('Average SONO LG vs Hop Cycle Percent')
hold off

figure(5)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.time;
y = mean_values.filtered_Buckle;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("Buckle")
legend(dataset_names, 'Location', 'best')
title('Average Buckle vs Hop Cycle Percent')
hold off

figure(6)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.time;
y = mean_values.filtered_EMG_1;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 1")
legend(dataset_names, 'Location', 'best')
title('Average EMG 1 vs Hop Cycle Percent')
hold off

figure(7)
hold on
for i = 1:numel(dataset_names)
    mean_values = mean_values_all{i};
    std_values = std_values_all{i};
    
x = mean_values.time;
y = mean_values.filtered_EMG_2;
plot(x,y)
end

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 2")
legend(dataset_names, 'Location', 'best')
title('Average EMG 2 vs Hop Cycle Percent')
hold off