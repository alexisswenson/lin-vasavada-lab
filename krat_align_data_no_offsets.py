# -*- coding: utf-8 -*-
"""
#######################################################
Aligning Kangaroo rat data
#######################################################

Lab affiliations for this project:
Dr. David Lin | Washington State University
Dr. Craig McGowan | University of Southern California

Required files:
    - Kinematics from DeepLabCut
    - BORIS labels
    - Trigger Offsets from scraping_text_from_videos.mlx
    - Features optional, if you don't have you can calulate by calling:
                
features = {file: kpm.create_feature_array(
    df, features_a, features_d) for file, df in kinematics.items()}            
    

This code is meant to align our lab's analog data with the corresponding
video data. This code assumes that short segments of videos are concatenated 
together into a larger video, and they have labels at the bottom of "frame" and
"trigger offset". This trigger offset allows us to match the video-derived
data with the analog data. The analog data contains a "trigger" column. When
the value of trigger goes above ~5 for the first time, this corresponds to when
the trigger offset on the video is 0.

########################################################
The challenges that this code was meant to overcome are:
######################################################## 

    - Different sampling rates. 
Analog will (should) be a higher sampling rate than the video sampling rate.

    - Different data collection windows
Analog has a larger data collection window, so we need to narrow this down by 
identifying the center of the window we need for an event, found by the first
instance where the "trigger" column > 5, and then pulling data from before and
after. The window is found from the trigger offsets file by finding the first
negative number (start of window) in the video sequence and the last positive 
number (end of window) in the video sequence. This code is in align_data() fnc.

    - Slight robustness to some of the idiosyncrasies I found in the data.
Not robust to new quirks in the data. The main ones I can note were the issue
of having multiple data files appended on top of each other in the analog .txt
files, and slight typos in the header columns. I used read_csv to read the .txt
files, but I realized somewhat later (working on a different project) that 
using something like np.loadtxt() would have been better. So if you run into
issues with loading the data, I reccomend making this change. It'll probably
require some decent restructuring of the data. This is probably close to
what you would want to load in all the data inside "raw_sensor_data" folder.
I tried this, and this breaks if there are multiple headers, so you'd need to 
fix that. 

raw_data = {f:np.loadtxt(os.path.join("raw_sensor_data",f),     
                  delimiter='\t', 
                  skiprows=6,     
                  dtype='object') for f in os.listdir("raw_sensor_data")}
data is tab delimited (maybe you don't need that though), header is 6 rows,
and you want an object datatype because you have a mix of datatypes (we have 
numeric and strings)

#######################################################
The merged data will have the following columns:
    t - time of the video, in analog sampling frequency
    data_file - the analog filename that the data is pulled from
    t_real - the time column from the analog file, to check where the data is
             coming from
    next few will be columns from analog data, self-explanatory
    
    (next columns are video-derived data which will be sparse, due to being 
     sampled less frequently)
    
    kinematic marker data from deeplabcut
    BORIS labels
    features (angles, distances, velocities)
        check the krat_pred_model_functions to see how these are calculated
    
#######################################################
If you want to understand how this works better, I *highly* reccomend using
a debugger to understand the different functions when it works. DD139 is a good
example of an animal that runs through this correctly. You can put breakpoints
at parts of a function to see what it is inputted, and what is happening.


All code was written by Nicholas Ozanich, working in Dr. Lin's lab at WSU.
For questions about the code, contact ozanichnicholas@gmail.com

#######################################################
Before you run, make sure your IDE is in the directory that contains all the
data folders and the krat_..._functions.py files because these are a local
import and require you to be in the same directory. I put most the functions
in those because they're several hundred lines of code and its cleaner to just
separate them.
#######################################################
"""
import krat_align_data_functions as kad
# you don't need the kpm functions, unless you want to calculate features here
#import krat_pred_model_functions as kpm
import numpy as np
import pandas as pd
import os


# These should be all you *should* have to change.
###############################################################################
animal_n = "DD118A"  # The string you input should match the folder in

# no_trigger_offsets_align_data
analog_sampling_rate = 0.00025      #DD has 300?
video_sampling_rate = 0.0025  # SCDD10 has 400 Hz video sampling rate
sampling_conversion_factor = int(video_sampling_rate / analog_sampling_rate)

BORIS = False
FEATURES = False
###############################################################################

# Some paths
# path to the animal folder
animal_path = os.path.join("no_trigger_offsets_align_data", animal_n)
animal_kinematics_path = os.path.join(
    animal_path, "kinematics")      # path to kinematics
# path to raw sensor .txt files
animal_sensor_path = os.path.join(animal_path, "raw_sensor_data")

# Load the kinematics (dict)
kinematics = {f.split('DLC')[0]: kad.dlc_to_xy(
    os.path.join(animal_kinematics_path, f)) for f in os.listdir(animal_kinematics_path)}
# Example for sorting by date and time
#from datetime import datetime

#kinematic_labels = [key for key in kinematics.keys()]

# kinematic_labels.sort(key=lambda date: datetime.strptime(date[:16], "%Y-%m-%d_%H-%M")) # year, month, d, hour, minute. check doc for strptime to see upper/lower differences

#new_k = {k: kinematics[k] for k in kinematic_labels}


# These are the video-derived datas you want to merge with analog.
# Delete entries you don't want.
df_to_merge = {"Kinematics": kinematics}

# If we have labels
if BORIS:
    BORIS_labels_path = os.path.join(animal_path, "BORIS")
    
    # Load BORIS labels from CSV files
    labels = {}
    for f in os.listdir(BORIS_labels_path):
        if f.endswith('.csv'):
            BORIS_df=pd.read_csv(os.path.join(BORIS_labels_path, f))
            BORIS_df.rename(columns={'Time': 't'}, inplace=True)
            labels[f.split('No')[0][:-1]] = BORIS_df
            
            # If 'Time' column exists in labels, rename it to 't'
            # if 'labels' in df_to_merge and 'Time' in df_to_merge['labels']:
            # df_to_merge['labels'].rename(columns={'Time': 't'}, inplace=True)
            # df_to_merge["Labels"] = labels


    # Correct the BORIS labels - see wiki for detailed explanation
    # You may want to comment this out, depending on state of labels
   # for key, value in labels.items():
    #    label = kad.correct_BORIS_labels(value, ['Stance'])  # .iloc[:-1,:]
     #   if kinematics[key].shape[0] != label.shape[0]:
      #      labels[key] = label.iloc[:-1, :]
       # else:
        #    labels[key] = label


# If we have features
if FEATURES:
    # Load feature arrays (dict)
    # (if you don't have, you need to run create them with the kpm functions)
    features = {f[:-13]: pd.read_csv(os.path.join('features', f)).drop('Unnamed: 0', axis=1).dropna()
                for f in os.listdir('features')}
    df_to_merge["Features"] = features


print(f"You are analyzing the kangaroo rat with tag:{animal_n}")
print(
    f"You provided a video sampling rate of {video_sampling_rate} (or {1/video_sampling_rate} Hz)")
print(
    f"You provided an analog sampling rate of {analog_sampling_rate} (or {1/analog_sampling_rate} Hz)")
print(
    f"The conversion factor between these two is {sampling_conversion_factor} analog samples for every video sample")
print("It is crucial to this script that these have common sampling rates. (the conversion factor should be an integer) \n")


# don't need these, just to clean up variable explorer
#del label, key, value
# %% Analog data from .txt files ----> Python dictionary data structure
# Read in the analog data (dict)
'''
Read in the raw data using pd.read_csv() which is robust to multiple headers in
the data, which we fix in our code using split_data_files(value). All the
functions are susceptible to change in experimental design, and were meant
to combat the idiosyncrasies in the data. But there will be new idiosyncrasies.

A note, lines that are getting the Range are unnecessary, unless you want the 
10V from the data.

Going forward, the raw data should probably be corrected before doing this, but
I coded it in instead. Its difficult to explain, but this is a place where I
reccomend using the debugger and breakpoints to see what is happening at each
step.
'''
raw_data = {f: pd.read_csv(os.path.join(animal_sensor_path, f))
            for f in os.listdir(animal_sensor_path)}

datas = {}
# Iterate over data from analog .txt files
for key, value in raw_data.items():
    # key: analog filename | value: data inside that file
    # split_data_file will split the data if there are multiple headers in it
    dfs = kad.split_data_files(value)
    for i, df in enumerate(dfs):
        # i: index (0:length(dfs)) | df: data corresponding to that header
        print(f"{key} | i={i} | length = {len(df)}")
        # Scrape header entries. Some are unused.
        Interval, ExcelDateTime, TimeFormat, DateFormat, ChannelTitle = kad.scrape_metadata(
            df)

        # Fix the channel title (very sensitive to change)
        ChannelTitle = kad.fix_ChannelTitle(ChannelTitle)

        # Getting the range of voltage values into a list
        Range = df.iloc[4].str.split().to_list()[0]

        # Fix Range
        Range = kad.fix_Range(Range)

        columns = ChannelTitle[1:]
        columns.insert(0, 't_real')

        # get rid of header
        df = df.iloc[5:, :]
        df = pd.DataFrame(df.iloc[:, 0].str.split().to_list(), columns=columns)

        # just getting the .txt filename without .txt, and i corresponds to
        # if we had multiple headers
        fmtstring = f"{key[:-4]}_{i}"
        # assign data to dictionary datas
        datas[fmtstring] = df

del Interval, ExcelDateTime, TimeFormat, DateFormat, ChannelTitle, Range,
columns, dfs, i, fmtstring
print("Done importing data and getting it into our data structures!")
# %% One video
''' 
Names of variables:
analog_sampling_rate - float, sampling rate of analog data in s
video_sampling_rate - float, sampling rate of video data in s
sampling_conversion_factor - int, # of analog samples for every video sample

df_to_merge - dict of Type of data: dataframe of data to merge, is flexible
              to different combinations of merged data
df_TOs - list of trigger offset filenames we are using
TOs - dict of trigger offset data
TO_0 - array of the frames when the trigger offsets == 0
trigged - dict of sensor data filesnames : list of indices where trigger offset
          first goes above 5.
analog_keys - list of analog data files
channels - all the channels from the sensor data
pre_trig - all indices where we have a positive value followed by a negative
           value in our trigger offsets. (UNCLEAN)
post_trig - all indices where we have a negative value followed by a positive
           value in our trigger offsets. (UNCLEAN)
pre_trigs - cleaned array of the start of trigger events (frame,offset)
post_trigs - cleaned array of the end of trigger events (frame,offset)
trig_indices - combined pre_trigs, and post_trigs

analog - array that is changed several times. These are the iterations it goes
through:
    1. Zeros array created with video-time column in first col (setup for align
                                                                ment)
    2. Adding 1's to the indices where the trigger offset (from video) == 0
                                                        (add_TO_to_start)
    3. Adding analog data to this array by centering the extraction on when
       the trigger (in raw data) first goes above 5, and pulling data that
       corresponds to the start of the event and end of the event from 
       trig_indices. We then assign this to the same range, but centered on
       when the trigger offset == 0. (align_data). Now ready to merge.
df_merged - dataframe of merged data, with the analog being the first part, and 
            followed by the data from df_to_merge. Flexible to varying 
            different datas as long as they are the same length.
'''


# Get the indices where the analog data is triggered
trigged = kad.get_triggered_times(datas,
                                  threshold=4,
                                  animal_n=animal_n)

#forDD118B
#trigged["DD118B_tr04_0"] = trigged["DD118B_tr04_0"][1:]

#forDD118A
trigged["DD118A_trial08_0"] = [trigged["DD118A_trial08_0"][0]]
trigged["DD118A_trial10_0"] = [trigged["DD118A_trial10_0"][0], trigged["DD118A_trial10_0"][3]]

# %%
# Get the analog filenames
analog_keys = list(trigged.keys())

# Columns from the data
channels = list(datas[analog_keys[0]].columns)

# Loop over events?


def setup_for_alignment_noTO(kinematics,
                             video_key,
                             video_sampling_rate,
                             analog_sampling_rate,
                             channels
                             ):
    # Create a time column to add to video-derived data

    # -1 to have 0 start
    video_time = kinematics[video_key].index.values

    video_time = video_time * video_sampling_rate  # multiply frames by sampling rate

    # round the floats for consistency
    video_time = [round(x, 5) for x in video_time]

    # Iterate over df you want to merge
    for key, df in kinematics.items():
        kinematics[video_key]['t'] = video_time

    # Get sampling conversion rate
    sampling_conv = int(video_sampling_rate/analog_sampling_rate)

    # Get the necessary length of an array to hold analog data that corresponds
    # to the video data

    video_len = len(kinematics[video_key])
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
    return analog, kinematics


def align_data(trigged, analog, datas, analog_sampling_rate, sampling_conversion_factor):

    pass


def merge_aligned_data_noTO(analog, df_to_merge, channels, filename):
    ''' 
    This is where we merge our data. 
    '''
    analog_columns = ['t', 'data_file'] + channels
    df_merged = pd.DataFrame(analog, columns=analog_columns)
    for key, df in df_to_merge.items():
        df_merged = pd.merge(df_merged, df, on='t', how='outer')
    df_merged.to_csv(filename)
    return df_merged


# List so we can iterate over these when we iterate over kinematics
# key, trigger activated idx for key, list of trigger activateds in trigged.items() for trigger activated in list
trigged_l = [(key, v) for key, value in trigged.items() for v in value]
''' 
trigged_l was hard to iteratively unpack when kept in the data structure it was
so this changes the data structure that is easier to unpack iteratively. 
change is just to have vid name, analog idx, rather than vid name, [analog idxs]
'''
print("Length of kinematics:", len(kinematics))
print("Length of trigged_l:", len(trigged_l))

analogs = []
# we iterate over the videos
#for i, kv in enumerate(kinematics.items()):                                           this was Nick's OG code
    
for i, (kv, trigged_item) in enumerate(zip(kinematics.items(), trigged_l)):
    """
    i is our iterator
    kv is a tuple with our (key,value), we unpack it immedietly though into k,v

    This loop will do the alignment video-by-video, compared to when we have
    trigger offsets where we want to combine everything together into one array
    these will all be separate arrays.

    """
    # This key should be the same for kinematics and BORIS
    # Kinematics key
    k = kv[0]

    # Kinematics data
    v = kv[1]

    # trigger data key
    d_k = trigged_l[i][0]

    # trigger data value
    d_v = trigged_l[i][1]

    # analog data to pull from
    df = datas[d_k]
    
    df = df.rename(columns={'t_real': 't'})
    
    print(i, k)
    print(d_k, d_v)

    analog, kinematics = setup_for_alignment_noTO(kinematics,
                                                  k,
                                                  video_sampling_rate,
                                                  analog_sampling_rate,
                                                  channels)
    
    # add in BORIS label to df_to_merge here
   # df_to_merge = {"kinematics": kinematics[k], "labels":labels[k]}
    df_to_merge = {"kinematics": kinematics[k]}
    

    # Need to add points where we will center on analog
    # add_TO_start_to_a

    v_len = (len(v)//2) * video_sampling_rate
    print(f"Video Length: {v_len}s")
    # # of analog frames to go back and forwarad
    f_pp = int(v_len / analog_sampling_rate)

    f_pre = d_v - f_pp  # number of analog data samples to go back
    # number of analog data samples to go forward
    f_post = d_v + f_pp + sampling_conversion_factor
    # Assign our window of analog data that matches our videos
    analog[:, 2:] = df.iloc[f_pre:f_post, :]
    analog[:, 1] = d_k      # assign video name to a column
    filename = f"no_trigger_offsets_align_data\DD118A\merged_data\{d_k}{i}.csv"
    
    if 't' in df.columns:
        print("DataFrame 'df' contains a column named 't'")
    else:
        print("DataFrame 'df' does not contain a column named 't'")
    
    # Merge aligned data
    df_merged = kad.merge_aligned_data(analog, df_to_merge, channels, filename)  # Properly indented
