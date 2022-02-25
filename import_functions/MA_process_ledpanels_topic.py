#import packages
import numpy as np
import warnings
import copy
import import_functions.MA_data_format_functions as ma_process

################################################################################
## get ledpanels info basic operations

#gets the indeces of the set pattern id command

def get_ledpanels_command_idxs(all_ledpanels_command):
    all_idx_pat_command = []
    all_idx_gain_command = []
    all_idx_stop_command = []
    for i in range(len(all_ledpanels_command)):
        idx_pat_command, idx_gain_command, idx_stop_command = get_idx_panels_commands(all_ledpanels_command[i])
        all_idx_pat_command.append(idx_pat_command)
        all_idx_gain_command.append(idx_gain_command)
        all_idx_stop_command.append(idx_stop_command)

    return all_idx_pat_command, all_idx_gain_command, all_idx_stop_command

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

################################################################################

#get times when ledpanels node sent a gain command

# def get_ledpanels_gain_te(all_idx_gain_command, all_ledpanels_te):
#     all_ledpanels_gain_te  = []
#     for i in range(len(all_idx_gain_command)):
#         ledpanels_gain_te = []
#         for j in range(len(all_idx_gain_command[i])):
#             t = all_ledpanels_te][all_idx_gain_command[i][j]]
#             ledpanels_gain_te.append(t)
#         all_ledpanels_gain_te.append(ledpanels_gain_te)
#     return all_ledpanels_gain_te

################################################################################

# #get arg1 for pannel_id command and gains
# #get the arg1 for the rows of pattern_id command
# all_pattern_id = []
# for i in range(len(all_idx_pat_command)):
#     pattern_id_file = []
#     for j in (all_idx_pat_command[i]):
#         arg1 = all_ledpanels_1[i][j]
#         pattern_id_file.append(arg1)
#     all_pattern_id.append(pattern_id_file)

################################################################################

def find_ledpanels_gains_value(all_ledpanels_1, all_idx_gain_command):
    #get the arg1 for the rows of gain command
    all_gains_values = []
    for i in range(len(all_idx_gain_command)):
        gains_file = []
        for j in (all_idx_gain_command[i]):
            arg1 = all_ledpanels_1[i][j]
            gains_file.append(arg1)
        all_gains_values.append(gains_file)
    return all_gains_values

################################################################################
# temporary fix for dark trials gain

def insert_dark_trial_gain(all_gains_values, dark_trials):
    #dark trials is list of dark trials
    all_gains = copy.deepcopy(all_gains_values)
    for i in range(len(all_gains_values)):
        for j in range(len(dark_trials)):
            all_gains[i].insert(dark_trials[j],0)
    return all_gains

################################################################################

def get_ledpanels_gains (all_ledpanels_command, all_ledpanels_1, dark_trials):
    all_idx_pat_command, all_idx_gain_command, all_idx_stop_command = get_ledpanels_command_idxs(all_ledpanels_command)
    all_gains_value = find_ledpanels_gains_value(all_ledpanels_1, all_idx_gain_command)
    all_gains = insert_dark_trial_gain(all_gains_value, dark_trials)
    return all_gains
