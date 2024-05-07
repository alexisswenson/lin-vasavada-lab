function [DD139_updated] = lowpass_filt(DD139)

%data that needs to be filtered
upperbackx = DD139.upper_backx;
upperbacky = DD139.upper_backy;
hipx = DD139.hipx;
hipy = DD139.hipy;
kneex = DD139.kneex;
kneey = DD139.kneey;
anklex = DD139.anklex;
ankley = DD139.ankley;
footx = DD139.footx;
footy = DD139.footy;
toetipx = DD139.toetipx;
toetipy = DD139.toetipy;

%lowpass filter information
fs = 500; % Sampling frequency in Hz
t = 0:1/fs:1; % Time vector
f = 25; % Desired cutoff frequency in Hz

% Design a lowpass Butterworth filter
order = 4 ; % Filter order
[b, a] = butter(order, f/(fs/2), 'low');

% Apply the filter to the signal
filtered_upperbackx = filtfilt(b, a, upperbackx);
filtered_upperbacky = filtfilt(b, a, upperbacky);
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

%assign column names and convert from array2table
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

% Concatenate the original table (DD139) with the new tables containing filtered EMG data
DD139_updated = [DD139, filtered_upperbackx_table, filtered_upperbacky_table, filtered_hipx_table, filtered_hipy_table, filtered_kneex_table, filtered_kneey_table, filtered_anklex_table, filtered_ankley_table, filtered_footx_table, filtered_footy_table, filtered_toetipx_table, filtered_toetipy_table];