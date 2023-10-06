import warnings
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import import_functions.MA_simple_import_functions_2 as ma_import
import import_functions.MA_data_format_functions as ma_process
import import_functions.MA_filter_functions as ma_filter
import import_functions.MA_process_ledpanels_topic as ma_ledpanels
import import_functions.MA_get_angvel as ma_get_angvel
import os
from fnmatch import fnmatch
import pickle

#data directory info

exp = 'MA_031021'

print(exp)
dataDir = '/Users/fponce/Documents/magno_arena_opto/MA_031021/data'

pattern_ma_data = "*.hdf5"

#processed data to be saved:
preprocessed_dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/'
pik_filename = preprocessed_dataDir+exp + '_preprocessed.p'
###############################################################################

# experiment info

time_each_trial = [180,180,900]

type_each_trial = ['dark','moving','dark']

expected_fps = 30

# some experiments have a all_ai_voltages topic, some don't
# 0: no all_ai_voltages topic saved
# 1: all_ai_voltages exists and needs to be analyzed
controller_voltages_rostopic = 0
###############################################################################
###############################################################################

#import data
datapaths = []
for path, subdirs, files in os.walk(dataDir):
    for name in files:
        if fnmatch(name, pattern_ma_data):
            datapaths.append(os.path.join(path, name))

print(str(len(datapaths))+' files')

###############################################################################
number_trials = len(time_each_trial)
experiment_time = np.sum(np.asarray(time_each_trial))
dark_trials = [i for i,x in enumerate(type_each_trial) if x=='dark']
###############################################################################
# get data from hdf5 to a dictionary with all the ros topics
# each key is a list of arrays, list len is number of files/flies
all_rostopics_rawdata = ma_import.get_topics_from_hdf5_to_dict(datapaths)
#print(list(all_rostopics_data.keys()))

# read lists from the dictionary
all_trial_index = all_rostopics_rawdata['all_trial_index']
all_ros_ts = all_rostopics_rawdata['all_ros_ts']
all_magnotether_tstamps = all_rostopics_rawdata['all_magnotether_ros_tstamps']
all_magnotether_angle = all_rostopics_rawdata['all_magnotether_angle']

all_ledpanels_ros_tstamps = all_rostopics_rawdata['all_ledpanels_ros_tstamps']
all_ledpanels_command = all_rostopics_rawdata['all_ledpanels_command']
all_ledpanels_1 = all_rostopics_rawdata['all_ledpanels_1']

###############################################################################
###############################################################################
# MA arena data topic analysis
# get start and end times/frames of each trial in all the different nodes that publish at different rates

# create list of lists with start and end trial times
# list lenght = number of files/flies

# get the ros time stamps that correspond  to start and end times of trials
all_start_times_ma, all_end_times_ma = ma_process.get_all_trial_start_n_end_times(number_trials,
                                                                                        all_ros_ts,
                                                                                        all_trial_index)

#get the closest times that correspond to start and end times in all_magnotether_tstamps topic
all_start_times_m, all_end_times_m = ma_process.get_closest_start_n_end_times (all_start_times_ma,
                                                                                all_end_times_ma,
                                                                                all_magnotether_tstamps)

all_start_frames_m, all_end_frames_m = ma_process.get_closest_start_n_end_frames (all_start_times_ma,
                                                                                    all_end_times_ma,
                                                                                    all_magnotether_tstamps)

print(str(len(all_start_frames_m[0]))+' start times/frames')
#print(all_start_frames_m[0])

###############################################################################
###############################################################################
# magnotether angle data topic analysis
# use the start/end times to get the magnotether angles from the start to end of experiment

all_magnotether_tstamps_exp = ma_process.get_start_end_exp_data(all_magnotether_tstamps,
                                                    all_start_frames_m,
                                                    all_end_frames_m)


all_magnotether_angles_exp = ma_process.get_start_end_exp_data(all_magnotether_angle,
                                                    all_start_frames_m,
                                                    all_end_frames_m)

###############################################################################
# interpolation and filtering of magnotether angles
# create regularly spaced time vector for data interpolation
t = [0, experiment_time]
reg_t = np.linspace(t[0], t[-1],(experiment_time*expected_fps)+1)

# create a regularly spaced time vector for each fly
# not used for interpolation
all_reg_ts = []
for i in range(len(all_magnotether_tstamps)):
    t = [all_magnotether_tstamps_exp[i][0], all_magnotether_tstamps_exp[i][-1]]
    reg_ts = np.linspace(t[0], t[-1],(experiment_time*expected_fps)+1)
    all_reg_ts.append(reg_t)

#get ellapsed time vectors for each fly
all_reg_ts_e = []
for i in range(len(all_magnotether_tstamps)):
    te = all_reg_ts[i] - all_reg_ts[i][0]
    all_reg_ts_e.append(te)

#interpolate magnotether angles
#get elapsed time for interpolation
all_magnotether_te_exp = ma_process.get_all_elapsed_time(all_magnotether_tstamps_exp)

all_magnotether_angles_exp_uw_interp = ma_process.get_all_magnotether_interp_angles (all_magnotether_angles_exp,
                                                                                    all_magnotether_te_exp,
                                                                                    reg_t)

#filter magnotether angles
all_angles = ma_filter.get_all_magnotether_filt_angles(all_magnotether_angles_exp_uw_interp)

#get filtered angular velocities
all_angvelos = ma_get_angvel.get_all_angvel(reg_t, all_magnotether_angles_exp_uw_interp)
###############################################################################

all_start_times_m_e = []
all_end_times_m_e = []
for i in range(len(all_start_times_m)):
    start_times_m_e = all_start_times_m[i] - all_start_times_m[i][0]
    end_times_m_e = all_end_times_m[i] - all_start_times_m[i][0]
    all_start_times_m_e.append(start_times_m_e)
    all_end_times_m_e.append(end_times_m_e)

all_start_times, all_end_times = ma_process.get_start_n_end_times_m (all_start_times_m_e, all_end_times_m_e, reg_t)
all_start_frames, all_end_frames = ma_process.get_start_n_end_frames_m (all_start_times_m_e, all_end_times_m_e, reg_t)

###############################################################################
###############################################################################
# ledpanels data topic analysis

#get elapsed time of panels
all_ledpanels_te = ma_process.get_all_elapsed_time(all_ledpanels_ros_tstamps)

all_gains = ma_ledpanels.get_ledpanels_gains (all_ledpanels_command, all_ledpanels_1, dark_trials)

###############################################################################
###############################################################################
#save processed data to dictionary

all_processed_rawdata = {}
all_processed_rawdata["experiment_name"] = exp
all_processed_rawdata["time_each_trial"] = time_each_trial
all_processed_rawdata["type_each_trial"] = type_each_trial
all_processed_rawdata["all_trial_start_times"] = all_start_times
all_processed_rawdata["all_trial_end_times"] = all_end_times
all_processed_rawdata["all_trial_start_frames"] = all_start_frames
all_processed_rawdata["all_trial_end_frames"] = all_end_frames
all_processed_rawdata["all_gains"] = all_gains
all_processed_rawdata["all_angles"] = all_angles
all_processed_rawdata["all_velos"] = all_angvelos
all_processed_rawdata["reg_t"] = reg_t
all_processed_rawdata["datapaths"] = datapaths

#save dictionary to a pickel
pickle.dump((all_processed_rawdata), open(pik_filename, "wb" ) )


print(str(exp)+' pre processed data saved to pickel')
