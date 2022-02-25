import warnings
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
import basic_analysis_functions.MA_basic_analysis_functions as ma_analysis
import basic_analysis_functions.MA_plot_functions_basic_analysis as ma_plot
import import_functions.MA_filter_functions as ma_filter
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
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

#plotting test trials
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
#         t = np.linspace(0,  len(v), len(v))
#         g = all_gains[i][j]
#         tt = type_each_trial[j]
#
#         fly_setpoint = all_trial_avs[i][trial- 2][2]
#
#         plt.plot(t,v, c='k', linewidth=1)
#         plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#         plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
#         plt.show()

###############################################################################
#saccade analysis

dt = reg_t[1] - reg_t[0]
fps = 1/dt

#overfiltering angular velocity
all_overfilt_velos = []
for i in range(len(all_velos)):
    overfilt_velo = ma_filter.overfiltVec(all_velos[i], dt)
    all_overfilt_velos.append(overfilt_velo)

#run saccade ider
all_nosaccs_velos = []
all_saccs_lr_idx = []
all_saccs_l_idx = []
all_saccs_r_idx = []
all_saccs_thresh = []
for i in range(len(all_velos)):
    #sacc detection threshold
    sdt = np.nanstd(all_velos[i] - all_overfilt_velos[i])*3
    headingVeloThresholdDeg = sdt
    #detrend angular velocity trace
    detrended_angvelo = all_velos[i] - all_overfilt_velos[i]
    #run saccade ider function on detrended angular velocity
    velo_detrended_nosaccs, saccs_lr, saccs_l, saccs_r = ma_sacc.findSacs(detrended_angvelo, all_velos[i], headingVeloThresholdDeg)
    #find and remove saccade from actual angular velocity signal
    velo_nosaccs = np.copy(all_velos[i])
    velo_nosaccs[np.isnan(velo_detrended_nosaccs)] = np.nan

    all_nosaccs_velos.append(velo_nosaccs)
    all_saccs_lr_idx.append(saccs_lr)
    all_saccs_l_idx.append(saccs_l)
    all_saccs_r_idx.append(saccs_r)
    all_saccs_thresh.append(headingVeloThresholdDeg)

for i in [0]:#range(len(datapaths)):
    fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
    print(fly)
    for j in range(len(all_posgain_bias_trials[i])):
        plt.figure(figsize=(14,7))

        trial = all_static_bias_trials[i][j]
        s = all_start_frames[i][trial - 2]
        e = all_end_frames[i][trial + 1]

        v = all_velos[i][s:e]
        v_ns = all_nosaccs_velos[i][s:e]
        v_overfilt = all_overfilt_velos[i][s:e]
        t = np.linspace(0,  len(v), len(v))*(1/fps)
        g = all_gains[i][j]
        tt = type_each_trial[j]

        fly_setpoint = all_trial_avs[i][trial- 2][2]

        #plt.plot(v - v_overfilt, c='magenta')
        #plt.plot(t,v, c='k', linewidth=1)
        plt.plot(t,v, c='k', linewidth=1, zorder = 3)
        plt.plot(t,v_ns, c='b', linewidth=1, zorder = 3)

        sdt = all_saccs_thresh[i]
        plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
        # plt.axhline(y=sdt, color='gray', linestyle='-', linewidth= 0.5)
        # plt.axhline(y=-sdt, color='gray', linestyle='-', linewidth= 0.5)
plt.show()

#plot saccadeless traces

    # plt.figure(figsize=(14,7))
    #
    # plt.plot(heading_noSacs, c='k', zorder =10000)
    # all_velos_sacs=np.copy(all_velos[i])*np.nan
    #
    # all_velos_sacs[~np.isnan(saccs)]=all_velos[i][~np.isnan(saccs)]
    # plt.plot(all_velos_sacs,  linestyle='None',  marker='.', color='b', markersize=4,    zorder=10000000)
    #
    #
    # plt.plot(all_velos[i], c='gray')
    #
    # plt.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
    # plt.axhline(y=sdt, color='gray', linestyle='-', linewidth= 0.5)
    # plt.axhline(y=-sdt, color='gray', linestyle='-', linewidth= 0.5)
    #
    # plt.plot(all_velos[i] - all_overfilt_velos[i], c='magenta')
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
