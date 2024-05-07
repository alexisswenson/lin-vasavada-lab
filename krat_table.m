function [upperback, hip, knee, ankle, foot, toetip] = krat_table(data)

% indexing by column name
upperback_xy = {'upper_backx', 'upper_backy'};
hip_xy = {'hipx', 'hipy'};
knee_xy = {'kneex', 'kneey'};
ankle_xy = {'anklex', 'ankley'};
foot_xy = {'footx', 'footy'};
toetip_xy = {'toetipx', 'toetipy'};

% kinematic data
upperback_table=[data(:,upperback_xy)];
hip_table=[data(:,hip_xy)];
knee_table=[data(:, knee_xy)];
ankle_table=[data(:,ankle_xy)];
kratfoot_table=[data(:,foot_xy)];
toetip_table=[data(:,toetip_xy)];

% table to array
upperback = table2array(upperback_table);
hip = table2array(hip_table);
knee = table2array(knee_table);
ankle = table2array(ankle_table);
foot = table2array(kratfoot_table);
toetip = table2array(toetip_table);