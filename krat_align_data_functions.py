# -*- coding: utf-8 -*-
"""
Functions to merge the data together
Created on Thu Jun 15 13:58:09 2023

@author: nicko
"""
import pandas as pd
import numpy as np
import os

# Trim the DeepLabCut outputted kinematics file #


def dlc_to_xy(path):
    ''' 
    Takes the path to the DeepLabCut file, and returns a pandas dataframe that
    only has the xy columns.
    '''
    # Read dataframe in
    df = pd.read_csv(path, header=[1, 2], dtype=np.float64)
    # Drop likelihood columns
    df_xy = df[df.columns.drop(list(df.filter(regex='likelihood')))]
    # columns in form markernamex and markernamey for all markers
    df_xy.columns = ['{}{}'.format(x[0], x[1]) for x in df_xy.columns]
    return df_xy.drop('bodypartscoords', axis=1)  # drop this extra column

# Correct BORIS labels


def correct_BORIS_labels(gt_gait, keys):
    '''
    Problem: the labels from BORIS are slightly off. There is a blank frame between stance
    and takeoff, and swing and landing. So this function just find value changes from 0 to 1
    and then if 0, it updates that index to be 1.

    This was necessary for our DD139 labels. I am unsure if this will be the 
    case for everyvideo. I believe that it is an artifact of BORIS exporting 
    w.r.t time instead of frames. 

    gt_gait is the BORIS labels file
    keys is the list of strings, where strings are the columns you need to 
    correct (should just be stance and/or swing)

    This function:
        1. iterates over keys provided
        2. iterates over rows in that key's column

    Since this is a script that looks at current and previous values, we have
    an exception for the first index, and we are looking at the difference
    between p and c. 
    '''
    df = gt_gait.copy()
    # Iterate over columns you want to correct
    for key in keys:
        # iterate over values in that column
        for idx, c in enumerate(df[key]):
            # if first iteration, p = c
            if idx == 0:
                p = c
            # if p not equal to c
            if p != c:
                # if c == 0, assign to current index 1
                if c == 0:
                    # print(f"Changing idx:{idx}")
                    df.loc[idx, key] = 1
                if c == 1:
                    # print(f"Changing idx: {idx}")
                    df.loc[idx-1, key] = 1
            # previous for next iteration is this one's current
            p = c
    return df

##############################################################
########### These are for scraping the .txt files ############


def split_data_files(d):
    '''
    It appears that the analog data files sometimes have multiple headers
    in the same .txt data file
    . This function splits these into seperate dfs
    and we will put those in a list for convenience.
    '''

    # Make copy
    df = d.copy()
    # Get column name of the Timeformat block
    bcolname = df.columns[0]
    # Find all rows with Timeformat, indicating a new trial of data collection

    idxs = df.loc[df[bcolname] == "TimeFormat=	StartOfBlock"]
    idxs = idxs.index.values

    idxs = [i - 1 for i in idxs]
    idxs.append(df.index[-1])
    dfs = []
    for i in range(len(idxs)-1):
        dfs.append(df.iloc[idxs[i]:idxs[i+1]-1, :])
    return dfs


def scrape_metadata(d):
    '''
    Getting data from the header, correcting where need be.

    When we read it in, the data is in form (n_rows,1), and each entry is
    delimited by a tab.
    So we take the entry,
    make sure it is a string d-type, (.str)
    split that string (it interprets the \t as the auto delimiter), (.split())
    convert that to a list of strings, (.to_list())
    and take the first string in that list. ([0])
    '''

    # Interval is interpreted as the column
    Interval = d.columns.str.split().to_list()[0]

    # ExcelDateTime is the first row
    ExcelDateTime = d.iloc[0].str.split().to_list()[0]

    # TimeFormat is the second row
    TimeFormat = d.iloc[1].str.split().to_list()[0]

    # DateFormat is the third row
    DateFormat = d.iloc[2].str.split().to_list()[0]

    # ChannelTitle is the 3rd row
    ChannelTitle = d.iloc[3].str.split('\t').to_list()[0]
    print('yeet')
    return Interval, ExcelDateTime, TimeFormat, DateFormat, ChannelTitle


def fix_ChannelTitle(ct):
    ''' Fix the ChannelTitle. Probably not robust to new data.'''

    first_entry = ct.pop(0)
    trigger = ct.pop(-1)
    new_cts = []
    p = 0
    for idx, i in enumerate(ct):
        # print(i,new_cts)
        if i != 'LG':
            new_cts.append(i)
        if i == 'LG':
            new_cts.remove(p)
            new_cts.append(p+'_'+i)
        p = i
    new_cts.insert(0, first_entry)
    new_cts.append(trigger)
    # print(new_cts)
    return new_cts


def fix_Range(r):
    '''
    Fix the range. Not sure how robust this is. All this is doing is since
    the 10 and V are different list entries, this just concatenates them

    '''
    first_entry = r.pop(0)
    new_r = []
    p = 0
    # enumerate list entries from range
    for idx, i in enumerate(r):
        if i != 'V':
            # if not "v" (then its a number, we append to our new_r)
            new_r.append(i)
        if i == 'V':
            # if "V" remove the previous entry (number) and make it of form:
            # "10V" for example
            new_r.remove(p)
            new_r.append(p+i)
        p = i
    # our loop requires looking at previous value, so the first value needs
    # to be accounted for. we do that here
    new_r.insert(0, first_entry)
    return new_r


#######################################################################
### Cleaning Trigger Offset Data ###
#######################################################################

def get_triggered_times(datas, threshold, animal_n):
    trigged = {}
    for k, v in datas.items():
        start_trigs = []
        if animal_n in k:
            df = pd.to_numeric(v['Trigger'])
            # 5 is a good threshold, usually trigger is at 6
            # SC videos go to zero when triggered
            # other videos go above 5 when triggered
            triggered = df[df > threshold]
            for i, d in enumerate(triggered.index):
                # if start of loop, we take that as the first trigger
                if i == 0:
                    start_trigs.append(d)
                    p = d
                # If the jump is more than 3 (this can be lower, mainly just seeing
                # if there is a jump)
                if abs(d - p) > 30:  # maybe?
                    start_trigs.append(d)
                p = d
            trigged[k] = start_trigs
    return trigged


def TO_sign_change(TO):
    '''
    This will iterate over TO's Trigger Offset column (2nd column). This
    accomplishes this:
        - Identifying frames/trigger offsets that could be the start/end of an
          window

    It only looks at the instance where the current value is negative, and the
    previous value is positive because this case will provide us both values
    we want.

    pre_trig will have all values where the current value is negative, and the
    previous value is positive.

    post_trig will have all values where the current value is positive, and the
    previous value is negative.

    We clean these arrays in the next functions called.

    abbreviations:
        p = previous value
        c = current value
    '''
    p = 0
    pre_trig = []
    post_trig = []
    for i, c in enumerate(TO.iloc[:, 1]):
        # if positive, and previous is 0 or neg, we append current to pre_trig
        # and append previous to post_trig. breaks at first use but we can pop it
        if c < 0 and p >= 0:

            pre_trig.append(TO.iloc[i, :].to_list())
            post_trig.append(TO.iloc[i-1, :].to_list())
        # if 0, append to both
        if c == 0:
            pre_trig.append(TO.iloc[i, :].to_list())
            post_trig.append(TO.iloc[i, :].to_list())
        # previous for next iter is current for this iter
        p = c

    # some formatting for triggers
    pre_trig = np.vstack(pre_trig)
    post_trig = np.vstack(post_trig)
    return pre_trig, post_trig


def clean_pre_trigs(pre_trig):
    '''
    Pre_trigs may include values that are not the start of the event window,
    and here we clean them up. We iterate over the length of pre_trig, and if
    the current value is negative and the previous is zero, then we append that
    to the cleaned list of pre_trigs.
    '''
    pre_trigs = []
    p = 0
    # this loop is to get rid of the false datas from the ocr messing up
    for i in range(len(pre_trig)):
        c = pre_trig[i, 1]
        # if current pos and previous neg or 0, append
        if c < 0 and p >= 0:
            pre_trigs.append(pre_trig[i, :])
        p = c
    pre_trigs = np.vstack(pre_trigs)
    return pre_trigs


def clean_post_trigs(post_trig, TO):
    '''
    Post_trigs may include values that are not the end of the event window,
    and here we clean them up. We iterate over the length of the post_trig, and
    if the current value is positive and the previous value is zero, then we
    append that to the cleaned list of post_trigs.
    '''
    # some formatting for post triggers
    post_trigs = []
    p = 0
    # this loop is to get rid of the false datas from the ocr messing up
    for i in range(len(post_trig)):
        c = post_trig[i, 1]
        if c > 0 and p <= 0:
            post_trigs.append(post_trig[i, :])
        p = c

    # Loop misses last value, append that row here
    post_trigs.append(TO.iloc[-1, :].to_numpy())
    # list of numpy arrays --> numpy array
    post_trigs = np.vstack(post_trigs)
    # post_trig has last entry in first row, so this removes it
    post_trigs = post_trigs[1:, :]
    return post_trigs


def combine_pre_post_triggers(pre_trigs, post_trigs):
    '''
    This combines the pre and post triggers we got. It also does a little
    organizing, cleaning, and getting it from ms to s so we can use it for
    aligning the data in seconds.
    '''
    # Concatenate
    trig_indices = np.vstack([pre_trigs,  post_trigs])
    # Sort by frame
    trig_indices = trig_indices[trig_indices[:, 0].argsort()]
    # Round (important to not have forever floats, and gets rid of the 0.1
    # issue when it displays on the video)
    trig_indices = trig_indices.round(decimals=0)
    # going from ms to s, to get same format as the sensor data
    trig_indices[:, 1] = trig_indices[:, 1] / 1000
    return trig_indices


def setup_for_alignment(df_to_merge, TO, TO_0, video_sampling_rate, analog_sampling_rate, channels, key):
    '''
    This takes in our data, and sets it up for alignment.
        - The first thing we do is get the frames of the video,
        - Subtracts 1 to start at 0,
        - Multiply those frame numbers by the video_sampling_rate to get
          the video time with analog sampling rate
        - Make sure these values are rounded (or else merging
          the data doesn't work when we have a long float),
        - Add this column to the video-derived data that we want to align with
          the analog data.

    Currently we have 8 columns necessary. These are (in order):
        1. Video time in analog sampling rate
        2. .txt filename
        3. t_real from the analog file (to go back and validate)
        4. SONO_LG
        5. Buckle
        6. EMG_1_LG
        7. EMG_2_LG
        8. Trigger

    So we create an (n,8) matrix, where n is the length of the frames of the
    video data, multiplied by the sampling conversion factor. This should
    be robust to adding a new channel, as the # of columns is calculated by
    the formula columns = 2 + # of channels in raw data

    I rounded to 5 decimal points, which is the minimum for the current analog
    samping frequency, but if you get equipment that has a sampling rate that
    requires more, you'll have to change it.
    '''
    # Create a time column to add to video-derived data

    video_time = TO['Frame'].to_numpy() - 1       # -1 to have 0 start

    video_time = video_time * video_sampling_rate  # multiply frames by sampling rate

    # round the floats for consistency
    video_time = [round(x, 5) for x in video_time]

    # Iterate over df you want to merge
    for key, df in df_to_merge.items():
        df_to_merge[key]['t'] = video_time
    # Get sampling conversion rate
    sampling_conv = int(video_sampling_rate/analog_sampling_rate)

    # Get the necessary length of an array to hold analog data that corresponds
    # to the video data

    video_len = len(df_to_merge["Kinematics"])
    video_len_in_analog_sampling_rate = video_len * sampling_conv

    # create array we will assign the analog values to
    analog = np.zeros(
        (video_len_in_analog_sampling_rate, 2 + len(channels)),
        dtype=object)         # object dtype for mixed data types

    # create an arrange of length video frames * sampling conversion rate
    # then convert this to time by multiplying by analog_sampling-rate
    t_f = np.arange(video_len_in_analog_sampling_rate) * analog_sampling_rate
    t_f = [round(x, 5) for x in t_f]
    analog[:, 0] = t_f                           # assign the time to analog

    # This just subtracts one from the frame number
    # it was written when I wanted it to be robust to videos that didn't
    # start at frame 1 (DD133_tread_slopes), but thats all this line does
    TO_0[:, 0] = TO_0[:, 0] - TO.iloc[0, 0]

    # Trigger offsets that line up with the frames in analog
    to_0 = TO_0 * sampling_conv
    return analog, df_to_merge, to_0


def setup_for_alignment_no_TO(df_to_merge, video_sampling_rate, analog_sampling_rate, channels, key):
    '''
    This takes in our data, and sets it up for alignment.
        - The first thing we do is get the frames of the video,
        - Subtracts 1 to start at 0,
        - Multiply those frame numbers by the video_sampling_rate to get
          the video time with analog sampling rate
        - Make sure these values are rounded (or else merging
          the data doesn't work when we have a long float),
        - Add this column to the video-derived data that we want to align with
          the analog data.

    Currently we have 8 columns necessary. These are (in order):
        1. Video time in analog sampling rate
        2. .txt filename
        3. t_real from the analog file (to go back and validate)
        4. SONO_LG
        5. Buckle
        6. EMG_1_LG
        7. EMG_2_LG
        8. Trigger

    So we create an (n,8) matrix, where n is the length of the frames of the
    video data, multiplied by the sampling conversion factor. This should
    be robust to adding a new channel, as the # of columns is calculated by
    the formula columns = 2 + # of channels in raw data

    I rounded to 5 decimal points, which is the minimum for the current analog
    samping frequency, but if you get equipment that has a sampling rate that
    requires more, you'll have to change it.
    '''
    # Create a time column to add to video-derived data
    keys = [k for k in df_to_merge["Kinematics"].keys()]

    video_time = TO['Frame'].to_numpy() - 1       # -1 to have 0 start
    video_time = df_to_merge["Kinematics"]

    video_time = video_time * video_sampling_rate  # multiply frames by sampling rate

    # round the floats for consistency
    video_time = [round(x, 5) for x in video_time]

    # Iterate over df you want to merge
    for key, df in df_to_merge.items():
        df_to_merge[key]['t'] = video_time
    # Get sampling conversion rate
    sampling_conv = int(video_sampling_rate/analog_sampling_rate)

    # Get the necessary length of an array to hold analog data that corresponds
    # to the video data

    video_len = len(df_to_merge["Kinematics"])
    video_len_in_analog_sampling_rate = video_len * sampling_conv

    # create array we will assign the analog values to
    analog = np.zeros(
        (video_len_in_analog_sampling_rate, 2 + len(channels)),
        dtype=object)         # object dtype for mixed data types

    # create an arrange of length video frames * sampling conversion rate
    # then convert this to time by multiplying by analog_sampling-rate
    t_f = np.arange(video_len_in_analog_sampling_rate) * analog_sampling_rate
    t_f = [round(x, 5) for x in t_f]
    analog[:, 0] = t_f                           # assign the time to analog

    # This just subtracts one from the frame number
    # it was written when I wanted it to be robust to videos that didn't
    # start at frame 1 (DD133_tread_slopes), but thats all this line does
    TO_0[:, 0] = TO_0[:, 0] - TO.iloc[0, 0]

    # Trigger offsets that line up with the frames in analog
    to_0 = TO_0 * sampling_conv
    return analog, df_to_merge, to_0


def add_TO_start_to_a(analog, to_0):
    '''
    This adds 1's to the second column of analog in the rows where the trigger
    is first triggered. So we add these in, then in align_data we use these
    values as the center of the events.

    analog is our (n,8) numpy array that has the video time in the analog
    sampling rate in the first column. The second column will have the 1's in
    the centers.

    to_0 is the rows of analog that the trigger is first triggered (i.e. center
    of our data)

    This function:
        - iterates over the length of analog rows
        - if the index is in to_0 ()

    '''
    # iterate over length of analog
    for i in range(len(analog[:, 0])):
        # if i (frame) is in to_0 we add this to the second column of analog
        if i in to_0[:, 0]:
            analog[i, 1] = 1
    return analog


def align_data(trigged, analog, datas, trig_indices, analog_sampling_rate, sampling_conversion_factor):
    '''
    This script assigns the analog data to our analog array. We pull out
    windows of data that correspond to the window of the event window. This is
    variable, as you can see in trig_indices. The start of the event window is
    usually -882, but not always. The end of the window is always 926 from my
    observations. But the above code *should* be robust to this number changing
    if it does in the future.

    I have print statements commented out for sanity, but go ahead and 
    uncomment them if you want to see what is happening or need to debug.

    Input variables:
        - datas is the dictionary of all analog.txt data
        - trigged is the dictionary of indices where analog data first
          goes above 5
        - analog is our analog array. it is inputted with only time and 1's
          in the first column where the trigger first is triggered.
        - trig_indices has the time values for the event window

    This is what the script does:
        1. Find rows where analog's first column == 1 (we center the trigger
           events around this value)
        2. Iterate over the analog data files from the analog data (trigged)
        3. Pull the dataframe associated with the analog data file
        4. Iterate over the triggered events for this data file
        5. The center is the value in analog_1. We then calculate the number
           of indices we need to go back (f_pre) and go ahead (f_post) from the
           df of the analog data.

    Frames_pre = time (s) before trigger offset == 0 / analog sampling rate
    Frames_pre = time (s) after trigger offset == 0 / analog sampling rate
    (i.e., s / frame * s^-1 == frame)

    We then find the f_pre for analog, and for the df from datas.
    The difference is that the analog variable uses cent for centering the
    assignment, and the df variable uses idx for centering.

    Additionally, we add the sampling_conversion_factor to pull that # of 
    frames after the post event, or else we have that # of frames missing from
    analog when we assign.

    '''
    # Find frames to center assignment to
    analog_1 = np.where(analog[:, 1] == 1)[0]

    # Initialize counter variable
    cnt = 0
    # Iterate over analog data files for this animal
    for k, v in trigged.items():
        # print(k, v)
        # pull data associate with this analog data file
        df = datas[k].astype('float64')
        # Iterate over the trigger events in this analog data file
        for idx in v:
            # print(a_1[cnt], cnt)

            # idx is the center inside df that we are pulling data from
            # cent is the center inside analog that we are assigning data to
            cent = analog_1[cnt]

            # calculating frames to go back (according to trigger offsets)
            f_pre = trig_indices[2*cnt, 1] / analog_sampling_rate

            # different frame numbers depending on which data
            f_pre_a = int(cent + f_pre)
            f_pre_d = int(idx + f_pre)

            # calculating frames to go forward (according to trigger offsets)
            f_post = trig_indices[2*cnt+1, 1] / analog_sampling_rate

            # different frame numbers depending on which data
            f_post_a = int(cent + f_post + sampling_conversion_factor)
            f_post_d = int(idx + f_post + sampling_conversion_factor)

            # assign data inside trigger event window from df to analog
            analog[f_pre_a:f_post_a, 2:] = df.iloc[f_pre_d:f_post_d, :]

            # assign data filename to the same window
            # Note: this overrides the column that has where the trigger
            # event starts.
            analog[f_pre_a:f_post_a, 1] = k

            # add one to our counter variable, necessary for analog_1 and
            # trig_indices to pull the correct event window.
            cnt += 1
    return analog


def merge_aligned_data(analog, df_to_merge, channels, filename):
    ''' 
    This is where we merge our data. 
    '''
    analog_columns = ['t', 'data_file'] + channels
    df_merged = pd.DataFrame(analog, columns=analog_columns)
    for key, df in df_to_merge.items():
        df_merged = pd.merge(df_merged, df, on='t', how='outer')
    df_merged.to_csv(filename)
    return df_merged