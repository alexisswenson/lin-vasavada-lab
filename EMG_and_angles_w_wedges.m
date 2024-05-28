clear
clc

% read in data from correct filepath
% data_nan = readtable("files for MATLAB/merged/DD139_slope_all trials_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD139_slope_all trials_merged(old).csv");
 data_nan = readtable("files for MATLAB/merged/DD_old/DD133_tread_slopes_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD133_slopes_merged.csv");

% data_nan = readtable("files for MATLAB/merged/DD130_slopes_01_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD130_slopes_02_merged.csv");
% data_nan = readtable("files for MATLAB/merged/DD135_slopes_02_merged.csv");

% data_nan = readtable("files for MATLAB/merged/SCDD03_hard.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD03_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand.csv");
% data_nan = readtable("files for MATLAB/merged/SCDD10_sand2.csv");

% cycles with steady hopping and wedges
[front_w_cycles] = front_wedge_cycles(data_nan);
[front_mean_w_cycles] = mean_wedge_cycles(front_w_cycles);
std_front_w_cycles = std(front_mean_w_cycles, 0, 1);

[back_w_cycles] = back_wedge_cycles(data_nan);
[back_mean_w_cycles] = mean_wedge_cycles(back_w_cycles);
std_back_w_cycles = std(back_mean_w_cycles, 0, 1);

% [w_cycles] = wedge_cycles(data_nan);

[EMG_front_w_cycles] = EMG_front_w_cycles(data_nan);
[front_mean_w_EMG] = EMG_mean_values(EMG_front_w_cycles);
std_front_w_emg = std(front_mean_w_EMG, 0, 1);

[EMG_back_w_cycles] = EMG_back_w_cycles(data_nan);
[back_mean_w_EMG] = EMG_mean_values(EMG_back_w_cycles);
std_back_w_emg = std(back_mean_w_EMG, 0, 1);

%cycles with steady hopping and no wedges
[cycles] = cycles(data_nan);
[mean_cycles] = mean_wedge_cycles(cycles);
std_cycles = std(mean_cycles, 0, 1);

[EMG_cycles] = EMG_cycles(data_nan);
[mean_EMG] = EMG_mean_values(EMG_cycles);
std_EMG = std(mean_EMG, 0, 1);


% Plotting the mean joint angles with error bars for +/- 1 standard deviation
% overlay individual steady hopping and wedge cycles
figure(1)
errorbar(mean_cycles.Time, mean_cycles.Hip_Angle, std_cycles.Hip_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Hip_Angle, -std_cycles.Hip_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(front_w_cycles)
    table_data = front_w_cycles{i};
    plot(table_data.time, table_data.hip_angle, 'b');
end

for i = 1:numel(back_w_cycles)
    table_data = back_w_cycles{i};
    plot(table_data.time, table_data.hip_angle, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("Hip Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Hip Angle vs Hop Cycle Percent')
hold off


figure(2)
errorbar(mean_cycles.Time, mean_cycles.Knee_Angle, std_cycles.Knee_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Knee_Angle, -std_cycles.Knee_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(front_w_cycles)
    table_data = front_w_cycles{i};
    plot(table_data.time, table_data.knee_angle, 'b');
end

for i = 1:numel(back_w_cycles)
    table_data = back_w_cycles{i};
    plot(table_data.time, table_data.knee_angle, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("Knee Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Knee Angle vs Hop Cycle Percent')
hold off

figure(3)
errorbar(mean_cycles.Time, mean_cycles.Ankle_Angle, std_cycles.Ankle_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Ankle_Angle, -std_cycles.Ankle_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(front_w_cycles)
    table_data = front_w_cycles{i};
    plot(table_data.time, table_data.ankle_angle, 'b');
end

for i = 1:numel(back_w_cycles)
    table_data = back_w_cycles{i};
    plot(table_data.time, table_data.ankle_angle, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("Ankle Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Ankle Angle vs Hop Cycle Percent')


figure(4)
errorbar(mean_cycles.Time, mean_cycles.Foot_Angle, std_cycles.Foot_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Foot_Angle, -std_cycles.Foot_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(front_w_cycles)
    table_data = front_w_cycles{i};
    plot(table_data.time, table_data.foot_angle, 'b');
end

for i = 1:numel(back_w_cycles)
    table_data = back_w_cycles{i};
    plot(table_data.time, table_data.foot_angle, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("Foot Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Foot Angle vs Hop Cycle Percent')



% plotting means of steady hopping and wedge cycles
figure(5)
errorbar(mean_cycles.Time, mean_cycles.Hip_Angle, std_cycles.Hip_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Hip_Angle, -std_cycles.Hip_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_cycles.Time, front_mean_w_cycles.Hip_Angle, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_cycles.Time, back_mean_w_cycles.Hip_Angle, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("Hip Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Hip Angle vs Hop Cycle Percent')


figure(6)
errorbar(mean_cycles.Time, mean_cycles.Knee_Angle, std_cycles.Knee_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Knee_Angle, -std_cycles.Knee_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_cycles.Time, front_mean_w_cycles.Knee_Angle, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_cycles.Time, back_mean_w_cycles.Knee_Angle, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("Knee Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Knee Angle vs Hop Cycle Percent')


figure(7)
errorbar(mean_cycles.Time, mean_cycles.Ankle_Angle, std_cycles.Ankle_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Ankle_Angle, -std_cycles.Ankle_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_cycles.Time, front_mean_w_cycles.Ankle_Angle, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_cycles.Time, back_mean_w_cycles.Ankle_Angle, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("Ankle Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Ankle Angle vs Hop Cycle Percent')


figure(8)
errorbar(mean_cycles.Time, mean_cycles.Foot_Angle, std_cycles.Foot_Angle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_cycles.Time, mean_cycles.Foot_Angle, -std_cycles.Foot_Angle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_cycles.Time, front_mean_w_cycles.Foot_Angle, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_cycles.Time, back_mean_w_cycles.Foot_Angle, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("Foot Angle")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Foot Angle vs Hop Cycle Percent')



% plotting mean and standard deviations of EMG data
figure(9)
errorbar(mean_EMG.time, mean_EMG.filtered_SONO, std_EMG.filtered_SONO, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_SONO, -std_EMG.filtered_SONO, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(EMG_front_w_cycles)
    emg_table_data = EMG_front_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_SONO, 'b');
end

for i = 1:numel(EMG_back_w_cycles)
    emg_table_data = EMG_back_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_SONO, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("SONO LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered SONO LG vs Hop Cycle Percent')

figure(10)
errorbar(mean_EMG.time, mean_EMG.filtered_Buckle, std_EMG.filtered_Buckle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_Buckle, -std_EMG.filtered_Buckle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(EMG_front_w_cycles)
    emg_table_data = EMG_front_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_Buckle, 'b');
end

for i = 1:numel(EMG_back_w_cycles)
    emg_table_data = EMG_back_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_Buckle, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("Buckle (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered Buckle vs Hop Cycle Percent')

figure(11)
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_1, std_EMG.filtered_EMG_1, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_1, -std_EMG.filtered_EMG_1, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(EMG_front_w_cycles)
    emg_table_data = EMG_front_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_EMG_1, 'b');
end

for i = 1:numel(EMG_back_w_cycles)
    emg_table_data = EMG_back_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_EMG_1, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 1 LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered EMG 1 vs Hop Cycle Percent')

figure(12)
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_2, std_EMG.filtered_EMG_2, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_2, -std_EMG.filtered_EMG_2, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

for i = 1:numel(EMG_front_w_cycles)
    emg_table_data = EMG_front_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_EMG_2, 'b');
end

for i = 1:numel(EMG_back_w_cycles)
    emg_table_data = EMG_back_w_cycles{i};
    plot(emg_table_data.time, emg_table_data.filtered_EMG_2, 'r');
end

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 2 LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered EMG 2 vs Hop Cycle Percent')


% plotting mean and standard deviations of EMG data
figure(13)
errorbar(mean_EMG.time, mean_EMG.filtered_SONO, std_EMG.filtered_SONO, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_SONO, -std_EMG.filtered_SONO, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_EMG.time, front_mean_w_EMG.filtered_SONO, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_EMG.time, back_mean_w_EMG.filtered_SONO, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("SONO LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered SONO LG vs Hop Cycle Percent')

figure(14)
errorbar(mean_EMG.time, mean_EMG.filtered_Buckle, std_EMG.filtered_Buckle, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_Buckle, -std_EMG.filtered_Buckle, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_EMG.time, front_mean_w_EMG.filtered_Buckle, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_EMG.time, back_mean_w_EMG.filtered_Buckle, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("Buckle (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered Buckle vs Hop Cycle Percent')

figure(15)
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_1, std_EMG.filtered_EMG_1, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_1, -std_EMG.filtered_EMG_1, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_EMG.time, front_mean_w_EMG.filtered_EMG_1, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_EMG.time, back_mean_w_EMG.filtered_EMG_1, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 1 LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered EMG 1 vs Hop Cycle Percent')

figure(16)
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_2, std_EMG.filtered_EMG_2, '-o', 'MarkerFaceColor', 'r', 'MarkerSize', 5);
hold on;
errorbar(mean_EMG.time, mean_EMG.filtered_EMG_2, -std_EMG.filtered_EMG_2, 'o', 'MarkerFaceColor', 'b', 'MarkerSize', 5);

plot(front_mean_w_EMG.time, front_mean_w_EMG.filtered_EMG_2, 'LineWidth', 3, 'Color', 'b');
plot(back_mean_w_EMG.time, back_mean_w_EMG.filtered_EMG_2, 'LineWidth', 3, 'Color', 'r');

xlabel("HOP CYCLE PERCENT")
ylabel("EMG 2 LG (Volts)")
legend('Front Wedge', 'Back Wedge', 'Location', 'best')
title('Average Filtered EMG 2 vs Hop Cycle Percent')