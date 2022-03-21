import warnings
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import basic_analysis_functions.MA_basic_analysis_functions as ma_analysis
import basic_analysis_functions.MA_plot_functions_basic_analysis as ma_plot
import import_functions.MA_filter_functions as ma_filter
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_run_sac_ID_ivo as ma_runsacc
###############################################################################

dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/'
filename = "MA_081621_10s_processed.p"
pickelDir = dataDir+filename

pickle_file = open(pickelDir, "rb")
all_processed_rawdata= []
all_processed_rawdata = ma_analysis.open_pickel(pickle_file)

g_t_av_f = -1.82 #gain to angular velocity factor

#are the trials randomized?
trials_rand = 1

# for key, value in all_processed_rawdata.items():
#     print (key)
###############################################################################
# read lists from the dictionary
all_velos = all_processed_rawdata["all_velos"]
all_angles = all_processed_rawdata["all_angles"]
reg_t = all_processed_rawdata["reg_t"]
all_start_times = all_processed_rawdata["all_trial_start_times"]
all_end_times = all_processed_rawdata["all_trial_end_times"]
all_start_frames = all_processed_rawdata["all_trial_start_frames"]
all_end_frames = all_processed_rawdata["all_trial_end_frames"]
all_gains = all_processed_rawdata["all_gains"]
type_each_trial = all_processed_rawdata["type_each_trial"]
datapaths = all_processed_rawdata["datapaths"]

###############################################################################
all_pattern_velos = []
for i in range(len(all_gains)):
    pattern_velos = np.array(all_gains[i])*g_t_av_f
    pattern_velos = np.where(pattern_velos==-0, 0, pattern_velos)
    all_pattern_velos.append(list(pattern_velos))

###############################################################################
# #plot raw data
# for i in range(len(datapaths)):
#
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     t = reg_t
#     av = all_velos[i]
#     st = all_start_times[i]
#     et = all_end_times[i]
#     g = all_gains[i]
#     tt = type_each_trial
#
#     ma_plot.plot_exp_angvelo(t, av, st, et, g, g_t_av_f, tt, fly)
#     plt.show()
###############################################################################
#sort data because trials were randomized

#find 'moving' trials
bias_trials = [i for i,x in enumerate(type_each_trial) if x=='moving']

#find moving trials with (0, +,- gains)
all_posgain_bias_trials = []
all_neggain_bias_trials = []
all_static_bias_trials = []

for i in range(len(all_gains)):
    mt = np.array(bias_trials)
    gains_test_trials = np.array(all_gains[i])[mt]

    posgain_bias_trials = mt[np.where(gains_test_trials==100)[0]]
    neggain_bias_trials = mt[np.where(gains_test_trials==-100)[0]]
    static_bias_trials = mt[np.where(gains_test_trials==0)[0]]

    all_posgain_bias_trials.append(posgain_bias_trials.tolist())
    all_neggain_bias_trials.append(neggain_bias_trials.tolist())
    all_static_bias_trials.append(static_bias_trials.tolist())

###############################################################################

###############################################################################
#find angular velocity medians in each trial
#
# t = reg_t
# av = median_ang_vel
# st = np.nanmean(all_start_times, axis=0)
# et = np.nanmean(all_end_times, axis=0)
# g = all_gains[0]
# tt = type_each_trial
#
# ma_plot.plot_angvelo(t, av, st, et, g, g_t_av_f , tt, 'mean')
###############################################################################

all_trial_avs = []
for i in range(len(all_velos)):
    trial_avs = []
    for j in range(len(all_start_frames[i])):

        s = all_start_frames[i][j]
        e = all_end_frames[i][j]

        v = all_velos[i][s:e]
        # get 3 trial angular velocities
        tav_1, tav_2, tav_3 = ma_analysis.get_3_med_ang_vels(v)
        trial_avs.append([tav_1, tav_2, tav_3])
    all_trial_avs.append(trial_avs)

###############################################################################
#saccade analysis

all_nosaccs_velos_2, all_saccs_id, all_saccs_l_idx, all_saccs_r_idx, all_overfilt_velos = ma_runsacc.get_all_nosacc_velos_loop(all_velos, reg_t)

all_nosaccs_velos = ma_runsacc.get_all_interp_after_saccs(all_nosaccs_velos_2, reg_t)
###############################################################################

for i in range(len(datapaths)):

    v = all_velos[i]
    v1 = all_nosaccs_velos[i]
    ov = all_overfilt_velos[i]
    test_trials = all_neggain_bias_trials[i]
    start_frames = all_start_frames[i]
    end_frames = all_end_frames[i]
    pattern_velos = all_pattern_velos[i]
    datapath = datapaths[i]

    #ma_plot.plot_ma081621_test_trials(v, v1, ov, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath, 0, i)


###############################################################################
#plot all trials and median
# test_trials_list = [all_posgain_bias_trials, all_neggain_bias_trials, all_static_bias_trials]
#
# for i in range(len(test_trials_list)):
#     test_trials = test_trials_list[i]
#
#     mm = ma_analysis.get_test_trial_median_ma081621(all_velos, test_trials, all_start_frames, all_end_frames)
#     ma_plot.plot_all_ma081621_test_trials(all_nosaccs_velos, reg_t, test_trials, all_start_frames, all_end_frames, all_pattern_velos, datapaths, mm, 1, i)
###############################################################################

#plot all trials and median
# test_trials_list = [all_posgain_bias_trials, all_neggain_bias_trials, all_static_bias_trials]
#
# for i in range(len(test_trials_list)):
#     test_trials = test_trials_list[i]
#
#     mm = ma_analysis.get_test_trial_median(all_velos, test_trials, all_start_frames, all_end_frames)
#     ma_plot.plot_all_ma081621_test_trials(all_nosaccs_velos, reg_t, test_trials, all_start_frames, all_end_frames, all_pattern_velos, datapaths, mm, 1, i)
###############################################################################
# all_rf = []
# for i in range(len(all_start_frames)):
#     trial_rf = []
#     for j in range(len(all_start_frames[i])):
#         if j>0:
#             s = all_start_frames[i][j] - 0
#         else:
#             s = all_start_frames[i][j]
#
#         e = all_end_frames[i][j]
#         ang_velo_trial = all_velos[i][s:e]
#         pat_ang_velo = all_pattern_velos[i][j]
#
#         rf = ma_analysis.get_trial_rise_frame(ang_velo_trial, pat_ang_velo)
#
#         trial_rf.append(rf+s)
#     all_rf.append(trial_rf)

###############################################################################
