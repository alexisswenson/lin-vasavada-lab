clc
clear

data = readtable("D:\DD118B - invivo sand\raw_sensor_data\DD118B_tr03.txt");  % solid

% %for one hop (0.2 seconds)
% EMG1 = data.EMG1LG(77540:78050, :);
% Buckle = data.Buckle(77540:78050, :);
% Sono = data.SonoLG(77540:78050, :);
% time = data.ChannelTitle_(77540:78050, :);

EMG1_all = data.EMG1LG;
Buckle_all = data.Buckle;
Sono_all = data.SonoLG;
time_all = data.ChannelTitle_;

% for multiple hops (1 second)
EMG1 = data.EMG1LG(78540:83050, :);
Buckle = data.Buckle(78540:83050, :);
Sono = data.SonoLG(78540:83050, :);

% Lowpass filter information
    fs = 4000; % Sampling frequency in Hz
    f = 50; % Desired cutoff frequency in Hz
    order = 2; % Filter order

    % Design a lowpass Butterworth filter
    [b, a] = butter(order, f/(fs/2), 'low');

    % Apply the lowpass filter to EMG 1, SONO, and Buckle
    filtered_EMG_1 = filtfilt(b, a, EMG1_all);
    filtered_SONO = filtfilt(b, a, Sono_all);
    filtered_Buckle = filtfilt(b, a, Buckle_all);

% Notch filter specifications
fs = 4000;   % Sampling frequency in Hz
f0 = 60;     % Center frequency of the notch filter in Hz
Q = 30;      % Quality factor of the notch filter

% Design the notch filter
d = fdesign.notch('N,F0,Q', 4, f0, Q, fs);
notchFilter = design(d);

% Apply the notch filter to filtered_EMG_1
filtered_EMG_1_notch = filtfilt(notchFilter.sosMatrix, notchFilter.ScaleValues, filtered_EMG_1);

% for multiple hops (1 second)
filtered_EMG_1 = filtered_EMG_1(78540:83050, :);
filtered_EMG_1_notch = filtered_EMG_1_notch(78540:83050, :);
time = data.ChannelTitle_(78540:83050, :);

figure(20)
plot(time, EMG1);
xlabel('Time [seconds]'); % Label for x-axis
ylabel('EMG1LG [Volts]'); % Label for y-axis
title('Raw EMG1LG vs. Time'); % Title of the plot

figure(21)
plot(time, filtered_EMG_1);
xlabel('Time [seconds]'); % Label for x-axis
ylabel('Low Pass EMG1LG [Volts]'); % Label for y-axis
title('Low Pass EMG1LG vs. Time'); % Title of the plot

figure(22)
plot(time, filtered_EMG_1_notch);
xlabel('Time [seconds]'); % Label for x-axis
ylabel('EMG1LG [Volts]'); % Label for y-axis
title('Low Pass and Notch EMG1LG vs. Time'); % Title of the plot

% figure(23)
% plot(time, filtered_Buckle);
% xlabel('Time [seconds]'); % Label for x-axis
% ylabel('Buckle [Volts]'); % Label for y-axis
% title('Buckle vs. Time'); % Title of the plot
% 
% figure(24)
% plot(time, filtered_SONO);
% xlabel('Time [seconds]'); % Label for x-axis
% ylabel('Sono [Volts]'); % Label for y-axis
% title('Sono vs. Time'); % Title of the plot
% 
% figure(25)
% plot(time, Buckle);
% xlabel('Time [seconds]'); % Label for x-axis
% ylabel('Buckle [Volts]'); % Label for y-axis
% title('Buckle vs. Time'); % Title of the plot
% 
% figure(26)
% plot(time, Sono);
% xlabel('Time [seconds]'); % Label for x-axis
% ylabel('Sono [Volts]'); % Label for y-axis
% title('Sono vs. Time'); % Title of the plot