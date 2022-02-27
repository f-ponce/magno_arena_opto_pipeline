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
# #plot test trials
# plt.figure()
# for i in range(len(datapaths)):
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     for j in range(len(all_posgain_bias_trials[i])):
#
#         trial = all_posgain_bias_trials[i][j]
#         s = all_start_frames[i][trial - 2]
#         e = all_end_frames[i][trial + 1]
#
#         v = all_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         plt.plot(t,v, c='k', linewidth=1)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
# plt.show()

###############################################################################
#find angular velocity medians in each trial

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

# for i in [0]:#range(len(datapaths)):
#
#     v = all_velos[i]
#     test_trials = all_posgain_bias_trials[i]
#     start_frames = all_start_frames[i]
#     end_frames = all_end_frames[i]
#     pattern_velos = all_pattern_velos[i]
#     datapath = datapaths[i]
#     ma_plot.plot_ma081621_test_trials(v, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath)

###############################################################################
#saccade analysis

dt = reg_t[1] - reg_t[0]
fps = 1/dt

#overfilter angular velocity data
highcut = 0.2
all_overfilt_velos_1 = ma_runsacc.get_all_overfiltered_velos(all_velos, dt, highcut)

#remove saccade
threshold = 3.5
all_nosaccs_velos_1, all_saccs_idx_1,  all_saccs_l_idx_1, all_saccs_r_idx_1, all_saccs_thresh_1 = ma_runsacc.get_all_nosacc_velos(all_velos,
                                                                                                    all_overfilt_velos_1, threshold)

###############################################################################
#
# for i in [0]:#range(len(datapaths)):
#
#     v = all_nosaccs_velos_1[i]
#     ov = all_overfilt_velos_1[i]
#     test_trials = all_posgain_bias_trials[i]
#     start_frames = all_start_frames[i]
#     end_frames = all_end_frames[i]
#     pattern_velos = all_pattern_velos[i]
#     datapath = datapaths[i]
#
#     ma_plot.plot_ma081621_test_trials(v, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath, ov)

###############################################################################
from scipy.interpolate import interp1d

all_nosaccs_velos_1_interp = []
for i in range(len(all_nosaccs_velos_1)):
    v1 = all_nosaccs_velos_1[i]
    v1_nan_idx = ~np.isnan(v1)

    v = np.array(v1[v1_nan_idx])
    t = np.array(reg_t[v1_nan_idx])

    f_a = interp1d(t, v, bounds_error=False)
    reg_a = f_a(reg_t)

    all_nosaccs_velos_1_interp.append(reg_a)

highcut = 0.2
all_overfilt_velos = ma_runsacc.get_all_overfiltered_velos(all_nosaccs_velos_1_interp, dt, highcut)

#remove saccade
threshold = 2.
all_nosaccs_velos, all_saccs_idx,  all_saccs_l_idx, all_saccs_r_idx, all_saccs_thresh = ma_runsacc.get_all_nosacc_velos(all_velos,
                                                                                                    all_overfilt_velos, threshold)

###############################################################################

for i in [0]:#range(len(datapaths)):

    v = all_nosaccs_velos[i]
    ov = all_overfilt_velos[i]
    test_trials = all_posgain_bias_trials[i]
    start_frames = all_start_frames[i]
    end_frames = all_end_frames[i]
    pattern_velos = all_pattern_velos[i]
    datapath = datapaths[i]
    overfilt_velos= all_overfilt_velos[i]
    ma_plot.plot_ma081621_test_trials(v, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath, ov)

# for i in [0]:#range(len(datapaths)):
#
#     v = all_velos[i]
#     test_trials = all_posgain_bias_trials[i]
#     start_frames = all_start_frames[i]
#     end_frames = all_end_frames[i]
#     pattern_velos = all_pattern_velos[i]
#     datapath = datapaths[i]
#     ma_plot.plot_ma081621_test_trials(v, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath)

#plotting removed saccades

# for i in [0]:#range(len(datapaths)):
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     for j in range(len(all_posgain_bias_trials[i])):
#         plt.figure(figsize=(14,7))
#
#         trial = all_static_bias_trials[i][j]
#         s = all_start_frames[i][trial - 2]
#         e = all_end_frames[i][trial + 1]
#
#         v = all_velos[i][s:e]
#         v_ns = all_nosaccs_velos[i][s:e]
#         v_overfilt = all_overfilt_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))*(1/fps)
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         fly_setpoint = all_trial_avs[i][trial- 2][2]
#
#         #plt.plot(v - v_overfilt, c='magenta')
#         #plt.plot(t,v, c='k', linewidth=1)
#         plt.plot(t,v, c='k', linewidth=1, zorder = 3)
#         plt.plot(t,v_ns, c='b', linewidth=1, zorder = 3)
#
#         sdt = all_saccs_thresh[i]
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         # plt.axhline(y=sdt, color='gray', linestyle='-', linewidth= 0.5)
#         # plt.axhline(y=-sdt, color='gray', linestyle='-', linewidth= 0.5)
# plt.show()

###############################################################################
# #plot filt - overfilt
# for i in [0]:#range(len(datapaths)):
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     for j in range(len(all_posgain_bias_trials[i])):
#         plt.figure(figsize=(14,7))
#
#         trial = all_neggain_bias_trials[i][j]
#         s = all_start_frames[i][trial - 2]
#         e = all_end_frames[i][trial + 1]
#
#         v = all_velos[i][s:e]
#         vf = all_overfilt_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))*(1/30.)
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         fly_setpoint = all_trial_avs[i][trial- 2][2]
#
#         plt.plot(t,v-vf, c='gray', linewidth=1)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         # plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
# plt.show()

# plt.plot(all_overfilt_velos[0])
# plt.show()
################################################################################
#plotting test trials angular vel over overfiltered ang vel

# for i in range(len(datapaths)):
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     for j in range(len(all_posgain_bias_trials[i])):
#         plt.figure(figsize=(14,7))
#
#         trial = all_neggain_bias_trials[i][j]
#         s = all_start_frames[i][trial - 2]
#         e = all_end_frames[i][trial + 1]
#
#         v = all_velos[i][s:e]
#         vf = all_overfilt_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))*(1/30.)
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         fly_setpoint = all_trial_avs[i][trial- 2][2]
#
#         plt.plot(t,v, c='gray', linewidth=1)
#         plt.plot(t,vf, c='k', linewidth=1)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
# plt.show()

###############################################################################
#plot filt - overfilt
# for i in range(len(datapaths)):
#     fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
#     print(fly)
#     for j in range(len(all_posgain_bias_trials[i])):
#         plt.figure(figsize=(14,7))
#
#         trial = all_neggain_bias_trials[i][j]
#         s = all_start_frames[i][trial - 2]
#         e = all_end_frames[i][trial + 1]
#
#         v = all_velos[i][s:e]
#         vf = all_overfilt_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))*(1/30.)
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         fly_setpoint = all_trial_avs[i][trial- 2][2]
#
#         plt.plot(t,v-vf, c='gray', linewidth=1)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         # plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
# plt.show()


# plt.figure()
# for i in range(len(all_posgain_bias_trials)):
#     c = ['k', 'g']
#     for j in range(len(all_posgain_bias_trials[i])):
#         trial = all_neggain_bias_trials[i][j]
#         pre = trial - 2
#         post = trial + 1
#
#         xt1 = random.uniform(0.0,0.4)
#         xt2 = random.uniform(0.6,1.0)
#
#         plt.scatter([xt1,xt1+2,xt1+4],all_trial_avs[i][pre], c=c[0], alpha = 0.5)
#         plt.scatter([xt2,xt2+2,xt2+4],all_trial_avs[i][post], c=c[1], alpha = 0.5)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         plt.ylim([-400,400])
#
#
# plt.figure()
# for i in range(len(all_posgain_bias_trials)):
#     c = ['k', 'g']
#     for j in range(len(all_posgain_bias_trials[i])):
#         trial = all_posgain_bias_trials[i][j]
#         pre = trial - 2
#         post = trial + 1
#
#         xt1 = random.uniform(0.0,0.4)
#         xt2 = random.uniform(0.6,1.0)
#
#         plt.scatter([xt1,xt1+2,xt1+4],all_trial_avs[i][pre], c=c[0], alpha = 0.5)
#         plt.scatter([xt2,xt2+2,xt2+4],all_trial_avs[i][post], c=c[1], alpha = 0.5)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         plt.ylim([-400,400])
# plt.show()



###############################################################################

###############################################################################

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
#     ma_plot.plot_angvelo(t, av, st, et, g, g_t_av_f, tt, fly)
#     plt.show()
###############################################################################

# median_ang_vel = ma_analysis.get_median_angvelos(all_velos)
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

# plt.figure()
#fly_color = ['k', 'r', 'b', 'g', 'y', 'm', 'c']
# for i in range(len(all_start_frames)):
#     for j in range(len(all_start_frames[i])):
#         s = all_start_frames[i][j]
#         e = all_end_frames[i][j]
#         ang_velo_trial = all_velos[i][s:e]
#         rf = all_rf[i][j]
#         plt.scatter(j, reg_t[rf]-reg_t[s], c = fly_color[i])
# plt.show()

###############################################################################
#plotting rt-- needs to be made a function at some point
#for i in range(len(all_start_frames)):
# for j in range(len(all_start_frames[i])):
#     for i in range(len(all_start_frames)):
#         s = all_start_frames[i][j]
#         e = all_start_frames[i][j] + 90
#         v = all_velos[i][s:e]
#         t = np.linspace(0,  len(v), len(v))
#         pat_ang_velo = all_pattern_velos[i][j]
#         rf = all_rf[i][j] - s
#         if pat_ang_velo<0:
#             plt.figure()
#             plt.plot(t,v, c = 'k')
#             plt.plot((t[0],t[-1]), (pat_ang_velo, pat_ang_velo), 'lime', linewidth= 1)
#             plt.scatter(t[rf], v[rf], c = 'r', s = 5, zorder=20)
#
# plt.show()
###############################################################################
