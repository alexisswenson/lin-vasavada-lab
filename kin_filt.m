function [data_updated] = kin_filt(data_nan)

% data that needs to be filtered
upper_backx = data_nan.upper_backx;
upper_backy = data_nan.upper_backy;
hipx = data_nan.hipx;
hipy = data_nan.hipy;
kneex = data_nan.kneex;
kneey = data_nan.kneey;
anklex = data_nan.anklex;
ankley = data_nan.ankley;
footx = data_nan.footx;
footy = data_nan.footy;
toetipx = data_nan.toetipx;
toetipy = data_nan.toetipy;
eyex = data_nan.eyex;
eyey = data_nan.eyey;


% lowpass filter information
fs = 500; % Sampling frequency in Hz
f = 30; % Desired cutoff frequency in Hz
% t = 0:1/fs:1; % Time vector

% Design a lowpass Butterworth filter
order = 4 ; % Filter order
[b, a] = butter(order, f/(fs/2), 'low');

% Apply the filter to the signal
filtered_upperbackx = filtfilt(b, a, upper_backx);
filtered_upperbacky = filtfilt(b, a, upper_backy);
filtered_hipx = filtfilt(b, a, hipx);
filtered_hipy = filtfilt(b, a, hipy);
filtered_kneex = filtfilt(b, a, kneex);
filtered_kneey = filtfilt(b, a, kneey);
filtered_anklex = filtfilt(b, a, anklex);
filtered_ankley = filtfilt(b, a, ankley);
filtered_footx = filtfilt(b, a, footx);
filtered_footy = filtfilt(b, a, footy);
filtered_toetipx = filtfilt(b, a, toetipx);
filtered_toetipy = filtfilt(b, a, toetipy);
filtered_eyex = filtfilt(b, a, eyex);
filtered_eyey = filtfilt(b, a, eyey);

% assign column names and convert from array2table
filtered_upperbackx_table = array2table(filtered_upperbackx, 'VariableNames', {'filtered_upperbackx'});
filtered_upperbacky_table = array2table(filtered_upperbacky, 'VariableNames', {'filtered_upperbacky'});
filtered_hipx_table = array2table(filtered_hipx, 'VariableNames', {'filtered_hipx'});
filtered_hipy_table = array2table(filtered_hipy, 'VariableNames', {'filtered_hipy'});
filtered_kneex_table = array2table(filtered_kneex, 'VariableNames', {'filtered_kneex'});
filtered_kneey_table = array2table(filtered_kneey, 'VariableNames', {'filtered_kneey'});
filtered_anklex_table = array2table(filtered_anklex, 'VariableNames', {'filtered_anklex'});
filtered_ankley_table = array2table(filtered_ankley, 'VariableNames', {'filtered_ankley'});
filtered_footx_table = array2table(filtered_footx, 'VariableNames', {'filtered_footx'});
filtered_footy_table = array2table(filtered_footy, 'VariableNames', {'filtered_footy'});
filtered_toetipx_table = array2table(filtered_toetipx, 'VariableNames', {'filtered_toetipx'});
filtered_toetipy_table = array2table(filtered_toetipy, 'VariableNames', {'filtered_toetipy'});
filtered_eyex_table = array2table(filtered_eyex, 'VariableNames', {'filtered_eyex'});
filtered_eyey_table = array2table(filtered_eyey, 'VariableNames', {'filtered_eyey'});

 % Concatenate the original table with the new tables containing filtered EMG data
    data_updated = [data_nan, filtered_eyey_table, filtered_eyex_table, filtered_toetipx_table, filtered_toetipy_table, filtered_footx_table, filtered_footy_table, filtered_anklex_table, filtered_ankley_table, filtered_kneex_table, filtered_kneey_table, filtered_hipx_table, filtered_hipy_table, filtered_upperbackx_table, filtered_upperbacky_table];