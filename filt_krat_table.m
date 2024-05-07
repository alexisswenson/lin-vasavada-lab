function [upperback, hip, knee, ankle, foot, toetip, eye] = filt_krat_table(data_filtered)

% indexing by column name
upperback_xy = {'filtered_upperbackx', 'filtered_upperbacky'};
hip_xy = {'filtered_hipx', 'filtered_hipy'};
knee_xy = {'kneex', 'kneey'};
ankle_xy = {'filtered_anklex', 'filtered_ankley'};
foot_xy = {'filtered_footx', 'filtered_footy'};
toetip_xy = {'filtered_toetipx', 'filtered_toetipy'};
eye_xy = {'filtered_eyex', 'filtered_eyey'};

% kinematic data
upperback_table=[data_filtered(:,upperback_xy)];
hip_table=[data_filtered(:,hip_xy)];
knee_table=[data_filtered(:, knee_xy)];
ankle_table=[data_filtered(:,ankle_xy)];
kratfoot_table=[data_filtered(:,foot_xy)];
toetip_table=[data_filtered(:,toetip_xy)];
eye_table=[data_filtered(:,eye_xy)];

% table to array
upperback = table2array(upperback_table);
hip = table2array(hip_table);
knee = table2array(knee_table);
ankle = table2array(ankle_table);
foot = table2array(kratfoot_table);
toetip = table2array(toetip_table);
eye = table2array(eye_table);