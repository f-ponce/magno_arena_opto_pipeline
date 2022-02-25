#import packages
import numpy as np
import copy
import pandas
import os
import imp
import pickle
from scipy.interpolate import interp1d
import warnings
import time
import matplotlib.pyplot as plt
import inspect
import types
import math as mat
import seaborn as sns
import matplotlib.gridspec as gridspec
from matplotlib.colors import LogNorm
import datetime
from matplotlib import animation
from IPython.display import HTML
from tempfile import NamedTemporaryFile
import matplotlib.backends.backend_pdf
from scipy import interpolate
from itertools import groupby
import re
import csv
import copy
from fnmatch import fnmatch
import json
import h5py
import copy
import math
pi = math.pi

##########################################################################################

#import functions
#used to import data to python notebooks

def get_topics_from_hdf5(datapaths):
#gets data from hdf5 file to arrays and generates a list per topic
#each list len is the number of files in dataset

    all_params_ts = []
    all_data_params = []

    all_ros_ts = []
    all_ts = []
    all_elapsed_time = []
    all_trial_index = []
    all_trial_elapsed_time = []
    all_angle_for_autostep = []
    all_init_angle = []
    all_autostep_running = []
    all_autostep_started = []
    all_autostep_stopped = []
    all_flow_running = []
    all_flow_started = []
    all_flow_stopped = []
    all_panels_running = []
    all_panels_started = []
    all_panels_stopped = []

    all_magnotether_angle = []
    all_magnotether_ros_tstamps = []
    all_magnotether_tstamps = []

    all_motion_ros_tstamps = []
    all_motion_tstamps = []
    all_motion_setpoint = []
    all_motion_position = []

    all_ledpanels_ros_tstamps = []
    all_ledpanels_command = []
    all_ledpanels_1 = []
    all_ledpanels_2 = []
    all_ledpanels_3 = []
    all_ledpanels_4 = []
    all_ledpanels_5 = []
    all_ledpanels_6 = []

    all_alicat_ros_tstamps = []
    all_alicat_devices = []

    all_sun_ros_tstamps = []
    all_sun_red = []
    all_sun_green = []
    all_sun_blue = []
    all_sun_message = []
    all_sun_led_number = []

    for i in range(len(datapaths)):
        f = h5py.File(datapaths[i], "r")

        #parameters topic
    #     params_ts = np.asarray(f['data_params_ros_tstamps'])
    #     data_params = np.asarray(f['data_params'])

        #virtual_desert topic
        ros_ts = np.asarray(f['ros_tstamps'])
        ts = np.asarray(f['tstamps'])
        elapsed_time = np.asarray(f['elapsed_time'])
        trial_index = np.asarray(f['current_trial_index'])
        trial_elapsed_time = np.asarray(f['trial_e_time'])
        angle_for_autostep = np.asarray(f['angle'])
        init_angle = np.asarray(f['init_angle'])

        #actions
        autostep_running = np.asarray(f['autostep_action_running'])
        autostep_started = np.asarray(f['autostep_action_started'])
        autostep_stopped = np.asarray(f['autostep_action_stopped'])

        flow_running = np.asarray(f['flow_action_running'])
        flow_started = np.asarray(f['flow_action_started'])
        flow_stopped = np.asarray(f['flow_action_stopped'])

        panels_running = np.asarray(f['panels_action_running'])
        panels_started = np.asarray(f['panels_action_started'])
        panels_stopped = np.asarray(f['panels_action_stopped'])

        #magnotether_angle topic
        magnotether_angle = np.asarray(f['magnotether_angle'])
        magnotether_ros_tstamps = np.asarray(f['magnotether_ros_tstamps'])
        magnotether_tstamps = np.asarray(f['magnotether_tstamps'])

        #motion_data topic
        motion_ros_tstamps = np.asarray(f['motion_data_ros_tstamps'])
        motion_tstamps = np.asarray(f['motion_data_tstamps'])
        motion_setpoint = np.asarray(f['motion_data_setpoint'])
        motion_position = np.asarray(f['motion_data_position'])

        #ledpanels topic
        ledpanels_ros_tstamps = np.asarray(f['ledpanels_ros_tstamps'])
        ledpanels_command = np.asarray(f['ledpanels_panels_command'])
        ledpanels_1 = np.asarray(f['ledpanels_panels_arg1'])
        ledpanels_2 = np.asarray(f['ledpanels_panels_arg2'])
        ledpanels_3 = np.asarray(f['ledpanels_panels_arg3'])
        ledpanels_4 = np.asarray(f['ledpanels_panels_arg4'])
        ledpanels_5 = np.asarray(f['ledpanels_panels_arg5'])
        ledpanels_6 = np.asarray(f['ledpanels_panels_arg6'])

        #alicat topic
        alicat_ros_tstamps = np.asarray(f['alicat_ros_tstamps'])
        alicat_devices = np.asarray(f['alicat_devices'])

        #sun topic
        sun_ros_tstamps = np.asarray(f['sun_ros_tstamps'])
        sun_red = np.asarray(f['sun_red'])
        sun_green = np.asarray(f['sun_green'])
        sun_blue = np.asarray(f['sun_blue'])
        sun_message = np.asarray(f['sun_message'])
        sun_led_number = np.asarray(f['sun_led_number'])

    #     all_params_ts.append(params_ts)
    #     all_data_params.append(data_params)

        all_ros_ts.append(ros_ts)
        all_ts.append(ts)
        all_elapsed_time.append(elapsed_time)
        all_trial_index.append(trial_index)
        all_trial_elapsed_time.append(trial_elapsed_time)
        all_angle_for_autostep.append(angle_for_autostep)
        all_init_angle.append(init_angle)
        all_autostep_running.append(autostep_running)
        all_autostep_started.append(autostep_started)
        all_autostep_stopped.append(autostep_stopped)
        all_flow_running.append(flow_running)
        all_flow_started.append(flow_started)
        all_flow_stopped.append(flow_stopped)
        all_panels_running.append(panels_running)
        all_panels_started.append(panels_started)
        all_panels_stopped.append(panels_stopped)
        all_magnotether_angle.append(magnotether_angle)
        all_magnotether_ros_tstamps.append(magnotether_ros_tstamps)
        all_magnotether_tstamps.append(magnotether_tstamps)
        all_motion_ros_tstamps.append(motion_ros_tstamps)
        all_motion_tstamps.append(motion_tstamps)
        all_motion_setpoint.append(motion_setpoint)
        all_motion_position.append(motion_position)
        all_ledpanels_1.append(ledpanels_1)
        all_ledpanels_2.append(ledpanels_2)
        all_ledpanels_3.append(ledpanels_3)
        all_ledpanels_4.append(ledpanels_4)
        all_ledpanels_5.append(ledpanels_5)
        all_ledpanels_6.append(ledpanels_6)
        all_ledpanels_command.append(ledpanels_command)
        all_ledpanels_ros_tstamps.append(ledpanels_ros_tstamps)
        all_alicat_ros_tstamps.append(alicat_ros_tstamps)
        all_alicat_devices.append(alicat_devices)
        all_sun_ros_tstamps.append(sun_ros_tstamps)
        all_sun_red.append(sun_red)
        all_sun_green.append(sun_green)
        all_sun_blue.append(sun_blue)
        all_sun_led_number.append(sun_led_number)

        return  all_ts,\
all_elapsed_time,\
all_trial_index,\
all_trial_elapsed_time,\
all_angle_for_autostep,\
all_init_angle,\
all_autostep_running,\
all_autostep_started,\
all_autostep_stopped,\
all_flow_running,\
all_flow_started,\
all_flow_stopped,\
all_panels_running,\
all_panels_started,\
all_panels_stopped,\
all_magnotether_angle,\
all_magnotether_ros_tstamps,\
all_magnotether_tstamps,\
all_motion_ros_tstamps,\
all_motion_tstamps,\
all_motion_setpoint,\
all_motion_position,\
all_ledpanels_1,\
all_ledpanels_2,\
all_ledpanels_3,\
all_ledpanels_4,\
all_ledpanels_5,\
all_ledpanels_6,\
all_ledpanels_command,\
all_ledpanels_ros_tstamps,\
all_alicat_ros_tstamps,\
all_alicat_devices,\
all_sun_ros_tstamps,\
all_sun_red,\
all_sun_green,\
all_sun_blue,\
all_sun_led_number

##########################################################################################

def get_all_trial_start_n_end_times(datapaths, number_trials, all_elapsed_time, all_trial_index):
#this gets the start and end times (using all_elapsed_time and all_trial_index)
#when the trial changes in the virtual desert node
#generates a list of arrays, list len is number of files
#datapaths is a list of paths to each file in dataset
#number_trials is the number of trials in experiment

    all_start_times = []
    all_end_times = []
    for i in range(len(datapaths)):
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

##########################################################################################

def get_all_trial_start_n_end_frames(datapaths, number_trials, all_trial_index):
#this gets the start and end frames (using all_trial_index)
#when the trial changes in the virtual desert node
#generates a list of arrays, list len is number of files
#datapaths is a list of paths to each file in dataset
#number_trials is the number of trials in experiment

    all_start_frames = []
    all_end_frames = []
    for i in range(len(datapaths)):
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

def get_elapsed_time(my_list):
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
#interpolation of magnotether angles

def get_all_magnotether_interp_angles (all_magnotether_angle, all_t_ellapsed, reg_t):
#all_magnotether_angle: list of magnotether angles
#all_t_ellapsed: list of elapsed times of each file
#reg_t: vector of evenly spaced time

    all_magnotether_interp_angles = []
    for i in range(len(all_magnotether_angle)):
        mysecs_np = all_t_ellapsed[i]
        myangles_np = all_magnotether_angle[i]
        f_a = interp1d(mysecs_np, myangles_np, bounds_error=False)
        reg_a = f_a(reg_t)
        all_magnotether_interp_angles.append(reg_a)
    return all_magnotether_interp_angles

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

def get_idx_panels_commands(ledpanels_command):
#gets the indeces of the set pattern id command and
# the gain command
    idx_pat_command = [i for i, x in enumerate(ledpanels_command)
               if x == b'set_pattern_id']
    idx_gain_command = [i for i, x in enumerate(ledpanels_command)
               if x == b'send_gain_bias']
    idx_stop_command = [i for i, x in enumerate(ledpanels_command)
               if x == b'stop']
    return idx_pat_command, idx_gain_command, idx_stop_command

##########################################################################################


import matplotlib.pyplot as plt
import matplotlib
import tempfile
import shutil
import os

###################################################################################################
# Adjust Spines (Dickinson style, thanks to Andrew Straw)
###################################################################################################

# NOTE: smart_bounds is disabled (commented out) in this function. It only works in matplotlib v >1.
# to fix this issue, try manually setting your tick marks (see example below)
def adjust_spines(ax,spines, spine_locations={}, smart_bounds=True, xticks=None, yticks=None, linewidth=1, spineColor='black'): # ivo: spineColor
    if type(spines) is not list:
        spines = [spines]

    # get ticks
    if xticks is None:
        xticks = ax.get_xticks()
    if yticks is None:
        yticks = ax.get_yticks()

#    spine_locations_dict = {'top': 10, 'right': 10, 'left': 10, 'bottom': 10}
    spine_locations_dict = {'top': 6, 'right': 6, 'left': 6, 'bottom': 6}
    for key in spine_locations.keys():
        spine_locations_dict[key] = spine_locations[key]

    if 'none' in spines:
        for loc, spine in ax.spines.iteritems():
            spine.set_color('none') # don't draw spine
        ax.yaxis.set_ticks([])
        ax.xaxis.set_ticks([])
        return

    for loc, spine in ax.spines.items():
        if loc in spines:
            spine.set_position(('outward',spine_locations_dict[loc])) # outward by x points
            spine.set_linewidth(linewidth)
            spine.set_color(spineColor) #ivo
            ax.tick_params(colors=spineColor) #ivo
            ax.tick_params(length=linewidth*4) #ivo
            ax.tick_params(pad=linewidth*4) #ivo
            ax.tick_params(direction="in") #ivo
        else:
            spine.set_color('none') # don't draw spine

    # smart bounds, if possible
    if int(matplotlib.__version__[0]) > 0 and smart_bounds:
        for loc, spine in ax.spines.items():
            if loc in ['left', 'right']:
                ticks = yticks
            if loc in ['top', 'bottom']:
                ticks = xticks
            spine.set_bounds(ticks[0], ticks[-1])

    # turn off ticks where there is no spine
    if 'left' in spines:
        ax.yaxis.set_ticks_position('left')
    elif 'right' in spines:
        ax.yaxis.set_ticks_position('right')
    else:
        # no yaxis ticks
        ax.yaxis.set_ticks([])

    if 'bottom' in spines:
        ax.xaxis.set_ticks_position('bottom')
    if 'top' in spines:
        ax.xaxis.set_ticks_position('top')
    else:
        # no xaxis ticks
        ax.xaxis.set_ticks([])

    if 'left' in spines or 'right' in spines:
        ax.set_yticks(yticks)
    if 'top' in spines or 'bottom' in spines:
        ax.set_xticks(xticks)

    for line in ax.get_xticklines() + ax.get_yticklines():
        #line.set_markersize(6)
        line.set_markeredgewidth(linewidth)


def kill_spines(ax):
    return adjust_spines(ax,'none',
                  spine_locations={},
                  smart_bounds=True,
                  xticks=None,
                  yticks=None,
                  linewidth=1)

def kill_labels(ax):
    #ax = ax['axis']
    for tl in ax.get_xticklabels() + ax.get_yticklabels():
            tl.set_visible(False)

def kill_all_spines(layout):
    [kill_spines(ax) for ax in layout.axes.values()]

def kill_all_labels(layout):
    [kill_labels(ax) for ax in layout.axes.values()]

def set_fontsize(fig,fontsize):
    """
    For each text object of a figure fig, set the font size to fontsize
    """
    def match(artist):
        return artist.__module__ == "matplotlib.text"

    for textobj in fig.findobj(match=match):
        textobj.set_fontsize(fontsize)

def set_fontfamily(fig,fontfamily): #ivo
    """
    For each text object of a figure fig, set the font size to fontsize
    """
    def match(artist):
        return artist.__module__ == "matplotlib.text"

    for textobj in fig.findobj(match=match):
        textobj.set_family(fontfamily)

def fix_mpl_svg(file_path, pattern='miterlimit:100000;', subst='miterlimit:1;'):
    """used to fix problematic outputs from the matplotlib svg
    generator, for example matplotlib creates exceptionaly large meterlimis"""
    fh, abs_path = tempfile.mkstemp()
    with open(abs_path,'w') as new_file:
        with open(file_path) as old_file:
            for line in old_file:
                new_file.write(line.replace(pattern, subst))
    os.close(fh)

    os.remove(file_path)

    shutil.move(abs_path, file_path)
    return
