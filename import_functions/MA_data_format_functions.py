#import packages
import numpy as np
from scipy.interpolate import interp1d
import warnings

##########################################################################################

def get_all_trial_start_n_end_times(number_trials, all_elapsed_time, all_trial_index):
#this gets the start and end times (using all_elapsed_time and all_trial_index)
#when the trial changes in the virtual desert node
#generates a list of arrays, list len is number of files
#datapaths is a list of paths to each file in dataset
#number_trials is the number of trials in experiment

    all_start_times = []
    all_end_times = []
    for i in range(len(all_elapsed_time)):
        start_times_trials = []
        end_times_trials = []
        for j in range(number_trials):
            start_time = all_elapsed_time[i][np.where(all_trial_index[i]==j)][0]
            end_time = all_elapsed_time[i][np.where(all_trial_index[i]==j)][-1]
            start_times_trials.append(start_time)
            end_times_trials.append(end_time)
        all_start_times.append(start_times_trials)
        all_end_times.append(end_times_trials)

    return all_start_times, all_end_times


def get_all_trial_start_n_end_frames(number_trials, all_trial_index):
#this gets the start and end frames (using all_trial_index)
#when the trial changes in the virtual desert node
#generates a list of arrays, list len is number of files
#datapaths is a list of paths to each file in dataset
#number_trials is the number of trials in experiment

    all_start_frames = []
    all_end_frames = []
    for i in range(len(all_trial_index)):
        start_frames_trials = []
        end_frames_trials = []
        for j in range(number_trials):
            start_frame = [np.where(all_trial_index[i]==j)][0][0][0]
            end_frame = [np.where(all_trial_index[i]==j)][0][0][-1]
            start_frames_trials.append(start_frame)
            end_frames_trials.append(end_frame)
        all_start_frames.append(start_frames_trials)
        all_end_frames.append(end_frames_trials)

    return all_start_frames, all_end_frames

##########################################################################################

# using the start and end times list, find the closest corresponding times in a different time topic

def get_closest_start_n_end_times (all_start_times, all_end_times, tt):
#get the closest times that correspond
#to start and end times (that come from vdesert node) in the reg_t
    all_start_times_m = []
    all_end_times_m = []
    for i in range(len(all_start_times)):
        start_times_trials_m = []
        end_times_trials_m = []
        for j in range(len(all_start_times[0])):
            start_times_m = find_nearest(tt[i], all_start_times[i][j])
            end_times_m = find_nearest(tt[i], all_end_times[i][j])
            start_times_trials_m.append(start_times_m)
            end_times_trials_m.append(end_times_m)
        all_start_times_m.append(start_times_trials_m)
        all_end_times_m.append(end_times_trials_m)
    return all_start_times_m, all_end_times_m

def get_closest_start_n_end_frames (all_start_times, all_end_times, tt):
#get the closest times that correspond
#to start and end times (that come from vdesert node) in the reg_t
    all_start_frames_m = []
    all_end_frames_m = []
    for i in range(len(all_start_times)):
        start_frames_trials_m = []
        end_frames_trials_m = []
        for j in range(len(all_start_times[0])):
            start_frames_m = find_nearest_idx(tt[i], all_start_times[i][j])
            end_frames_m = find_nearest_idx(tt[i], all_end_times[i][j])
            start_frames_trials_m.append(start_frames_m)
            end_frames_trials_m.append(end_frames_m)
        all_start_frames_m.append(start_frames_trials_m)
        all_end_frames_m.append(end_frames_trials_m)
    return all_start_frames_m, all_end_frames_m

##########################################################################################

def get_all_elapsed_time(my_list):
#for list of lists
#getting elapsed time of time stamps to use
    all_t_ellapsed = []
    for i in range(len(my_list)):
        t_ellapsed = my_list[i] - my_list[i][0]
        all_t_ellapsed.append(t_ellapsed)
    return all_t_ellapsed

def get_elapsed_time_2lists(my_list_1, my_list_2):
#for list of arrays
#gets elapsed time of each array in list and first element
#of each array in my_list_2
    all_t_ellapsed = []
    for i in range(len(my_list_1)):
        t_ellapsed = my_list_1[i] - my_list_2[i][0]
        all_t_ellapsed.append(t_ellapsed)
    return all_t_ellapsed

##########################################################################################

def find_nearest(array, value):
#for list of lists
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

##########################################################################################

def get_start_n_end_times_m (all_start_times, all_end_times, reg_t):
#get the closest times that correspond
#to start and end times (that come from vdesert node) in the reg_t
    all_start_times_m = []
    all_end_times_m = []
    for i in range(len(all_start_times)):
        start_times_trials_m = []
        end_times_trials_m = []
        for j in range(len(all_start_times[0])):
            start_times_m = find_nearest(reg_t, all_start_times[i][j])
            end_times_m = find_nearest(reg_t, all_end_times[i][j])
            start_times_trials_m.append(start_times_m)
            end_times_trials_m.append(end_times_m)
        all_start_times_m.append(start_times_trials_m)
        all_end_times_m.append(end_times_trials_m)
    return all_start_times_m, all_end_times_m

def get_start_n_end_frames_m (all_start_times, all_end_times, reg_t):
#get the closest times that correspond
#to start and end times (that come from vdesert node) in the reg_t
    all_start_frames_m = []
    all_end_frames_m = []
    for i in range(len(all_start_times)):
        start_frames_trials_m = []
        end_frames_trials_m = []
        for j in range(len(all_start_times[0])):
            start_frames_m = find_nearest_idx(reg_t, all_start_times[i][j])
            end_frames_m = find_nearest_idx(reg_t, all_end_times[i][j])
            start_frames_trials_m.append(start_frames_m)
            end_frames_trials_m.append(end_frames_m)
        all_start_frames_m.append(start_frames_trials_m)
        all_end_frames_m.append(end_frames_trials_m)
    return all_start_frames_m, all_end_frames_m

##########################################################################################

# use the start/end times to get the magnotether angles from the start to end of experiment

def get_start_end_exp_data(all_data_list, all_start_frames, all_end_frames):

    all_data_list_exp = []
    for i in range(len(all_data_list)):
        # get the angles and timestamps from start to end of experiment
        data_list_exp = all_data_list[i][all_start_frames[i][0]:all_end_frames[i][-1]]
        all_data_list_exp.append(data_list_exp)

    return all_data_list_exp

##########################################################################################

def get_first_n_last_minute_trial (all_start_frames_m, all_end_frames_m, number_frames_per_sec):
#get start frames for last minute and end frames first min
    all_start_frames_m_lm = []
    all_end_frames_m_fm = []
    for i in range(len(all_start_frames_m)):
        start_frames_m_lm_trials = []
        end_frames_m_fm_trials = []
        for j in range(len(all_start_frames_m[0])):
            start_frames_m_lm = all_end_frames_m[i][j] - number_frames_per_sec*60 #getting start frames last min
            end_frames_m_fm = all_start_frames_m[i][j] + number_frames_per_sec*60 #getting end frames first min
            start_frames_m_lm_trials.append(start_frames_m_lm)
            end_frames_m_fm_trials.append(end_frames_m_fm)
        all_start_frames_m_lm.append(start_frames_m_lm_trials)
        all_end_frames_m_fm.append(end_frames_m_fm_trials)
    return all_start_frames_m_lm, all_end_frames_m_fm

##########################################################################################

##########################################################################################
#interpolation of magnotether angles

def get_all_magnotether_interp_angles (all_magnotether_angle, all_t_ellapsed, reg_t):
#all_magnotether_angle: list of magnotether angles
#all_t_ellapsed: list of elapsed times of each file
#reg_t: vector of evenly spaced time

    all_magnotether_interp_angles = []
    for i in range(len(all_magnotether_angle)):
        mysecs = all_t_ellapsed[i]
        myangles = all_magnotether_angle[i]
        myangles_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(myangles)))
        f_a = interp1d(mysecs, myangles_unwrapped, bounds_error=False)
        reg_a = f_a(reg_t)
        all_magnotether_interp_angles.append(reg_a)
    return all_magnotether_interp_angles

##########################################################################################
