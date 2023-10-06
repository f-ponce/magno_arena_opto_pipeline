import warnings
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
from collections import Counter
import copy
from scipy.interpolate import interp1d

import basic_analysis_functions.MA_basic_analysis_functions as ma_analysis
import plot_functions.MA_plot_functions_basic_analysis as ma_plot
import plot_functions.MA_plot_functions_MA_081621 as ma_exp_plot
import plot_functions.plot_magno_funcs as ma_exp_plot2
import import_functions.MA_filter_functions as ma_filter
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_run_sac_ID_ivo as ma_runsacc
###############################################################################

exp_s = 100

#exp = 'MA_081621_'+str(exp_s)+'s'
exp = 'MA_032122_'+str(exp_s)+'s'

filename = exp+'_preprocessed.p'

dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/'
pickelDir = dataDir+filename

pickle_file = open(pickelDir, "rb")
all_processed_rawdata = []
all_processed_rawdata = ma_analysis.open_pickel(pickle_file)

#are the trials randomized?
trials_rand = 1
fps = 30
# for key, value in all_processed_rawdata.items():
#     print (key)
plot_raw_data = 0
#processed data to be saved:
save_data = 0
processed_dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_processed_data/'
pik_filename = processed_dataDir+exp + '_processed.p'

###############################################################################
# read lists from the dictionary
exp = all_processed_rawdata["experiment_name"]
type_each_trial = all_processed_rawdata["type_each_trial"]
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
g_t_av_f = -1.82 #gain to angular velocity factor

all_pattern_velos = []
for i in range(len(all_gains)):
    pattern_velos = np.array(all_gains[i])*g_t_av_f
    pattern_velos = np.where(pattern_velos==-0, 0, pattern_velos)
    all_pattern_velos.append(list(pattern_velos))

###############################################################################
#plot raw data
if plot_raw_data == 1:

    for i in range(len(datapaths)):

        fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])
        print(fly)
        t = reg_t
        av = all_velos[i]
        st = all_start_times[i]
        et = all_end_times[i]
        g = all_gains[i]
        tt = type_each_trial

        ma_plot.plot_exp_angvelo(t, av, st, et, g, g_t_av_f, tt, fly)
        plt.show()
else:
    pass
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
#get a single start/end frame number for all flies in exp, for plotting etc
exp_start_frames = []
exp_end_frames = []
for j in range(len(all_start_times[0])):
    s_trial_frames = []
    e_trial_frames = []
    for i in range(len(all_start_times)):
        s_trial_frames.append(all_start_frames[i][j])
        e_trial_frames.append(all_end_frames[i][j])
    s_freq = Counter(s_trial_frames)
    e_freq = Counter(e_trial_frames)
    s_frame = s_freq.most_common()
    e_frame = e_freq.most_common()
    exp_start_frames.append(s_frame[0][0])
    exp_end_frames.append(e_frame[0][0])
# print(exp_start_frames)
# print(exp_end_frames)
###############################################################################
#find and remove saccades

all_nosaccs_velos_2, all_saccs_id, all_saccs_l_idx, all_saccs_r_idx, all_overfilt_velos = ma_runsacc.get_all_nosacc_velos_loop(all_velos, reg_t)
all_nosaccs_velos = ma_runsacc.get_all_interp_after_saccs(all_nosaccs_velos_2, reg_t)

###############################################################################
#find angular velocity medians in each trial

all_trial_avs = []
for i in range(len(all_nosaccs_velos)):
    trial_avs = []
    for j in range(len(exp_start_frames)):

        s = exp_start_frames[j]
        e = exp_end_frames[j]

        v = all_nosaccs_velos[i][s:e]
        # get 3 trial angular velocities
        tav_1, tav_2, tav_3 = ma_analysis.get_3_med_ang_vels(v)
        trial_avs.append([tav_1, tav_2, tav_3])
    all_trial_avs.append(trial_avs)
###############################################################################
#plot single trial, single fly

for i in [0]:#range(len(datapaths)):

    v = all_velos[i]
    v1 = all_nosaccs_velos_2[i]
    ov = all_overfilt_velos[i]
    test_trials = all_static_bias_trials[i]
    start_frames = all_start_frames[i]
    end_frames = all_end_frames[i]
    pattern_velos = all_pattern_velos[i]
    datapath = datapaths[i]
    sl = all_saccs_l_idx[i]
    sr = all_saccs_r_idx[i]

    # trial_sl_n, trial_sr_n, sacBins, binW = ma_runsacc.sacHist_per_trial(reg_t, test_trials, start_frames, end_frames, sl, sr)
    #
    # plt.figure()
    # plt.bar (sacBins, trial_sl_n[0], binW, color = 'blue', ec = 'blue', zorder=10 )
    # plt.bar (sacBins, trial_sl_n[1], binW, color = 'gray', ec = 'gray', zorder=10 )
    # plt.show()
    #
    # plt.figure()
    # plt.bar (sacBins, trial_sr_n[0], binW, color = 'red', ec = 'red', zorder=10 )
    # plt.bar (sacBins, trial_sr_n[1], binW, color = 'magenta', ec = 'magenta', zorder=10 )
    # plt.show()

    # ss = np.abs((reg_t)[~np.isnan(sl)])
    # n, bins = np.histogram(ss, 120, density=False)
    # binW = bins[2] - bins[1]
    # hlfBinWidth = binW/2.
    # sacBins = bins[:-1]+hlfBinWidth
    # plt.bar (sacBins, n, binW, color = 'blue', ec = 'gray', zorder=10 )
    # plt.show()
    # print(sacBins)
    # print(n)

    ##################
    # ss = np.abs((reg_t)[~np.isnan(sl)])
    # n, bins, patches = plt.hist(ss, 120, density=False)
    #
    # print(np.sum((n/np.sum(n))))
    # if i==0:
    #     n_array = (n/np.sum(n))[...,np.newaxis]
    # else:
    #     n_array = np.concatenate((n_array, (n/np.sum(n))[...,np.newaxis]),axis=1)
    # print(n_array.shape)
    # if i==len(datapaths)-1:
    #     binW=bins[2]-bins[1]
    #     hlfBinWidth = binW/2.
    #     sacHistN = n_array.mean(axis=1)
    #     sacBins = bins[:-1]+hlfBinWidth
    #     binTime=(sacBins/sacBins[-1])*(reg_t[-1]/60.)
    #     binTimeW=binTime[2]-binTime[1]
    #     print(reg_t/60.)
    #     plt.figure()
    #     #plt.bar (sacBins, sacHistN, binW, color = 'red', ec = "none", zorder=10 )
    #     plt.bar (binTime, sacHistN, binTimeW, color = 'red', ec = "none", zorder=10 )
    #     plt.show()

    # ma_exp_plot.plot_ma081621_test_trials(v, v1, ov, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath, 0, i)
    # ma_exp_plot.plot_ma081621_test_trials_saccs(v, v1, ov, reg_t, test_trials, start_frames, end_frames, sl, sr,pattern_velos, datapath, 0, i)


###############################################################################

#plot all trials, all flies and median

test_trials_list = [all_posgain_bias_trials, all_neggain_bias_trials, all_static_bias_trials]

for i in range(len(test_trials_list)):
    test_trials = test_trials_list[i]

    mm, ci = ma_analysis.get_test_trial_median_ma081621(all_nosaccs_velos, test_trials, exp_start_frames, exp_end_frames)

    # ma_exp_plot.plot_all_ma081621_test_trials(all_nosaccs_velos, reg_t, test_trials, exp_start_frames, exp_end_frames, all_pattern_velos,
    #                                        datapaths, mm, 0, i)
    # ma_exp_plot.plot_all_ma081621_test_trials_perfly(all_nosaccs_velos, reg_t, test_trials, exp_start_frames, exp_end_frames, all_pattern_velos,
    #                                        datapaths, mm, 0, i)
###############################################################################

# #plot all test trials with 95 ci filled in
test_trials_list = [all_posgain_bias_trials, all_neggain_bias_trials]#, all_static_bias_trials]
#test_trials_list = [all_static_bias_trials]

test_trials_mmeans = []
test_trials_ci_95s = []
for i in range(len(test_trials_list)):
    test_trials = test_trials_list[i]

    mm_sd, ci_95_sd = ma_analysis.get_test_trial_median_ma081621_2(all_nosaccs_velos, test_trials, exp_start_frames, exp_end_frames)
    actual_time = np.linspace(0, (120+30+exp_s), (120+30+exp_s)*30)

    if len(mm_sd)>=len(actual_time):
        mm_sdd = mm_sd[0:len(actual_time)]
        ci_95_sdd0 = ci_95_sd[0][0:len(actual_time)]
        ci_95_sdd1 = ci_95_sd[1][0:len(actual_time)]
    else:
        mm_sdd = mm_sd
        ci_95_sdd0 = ci_95_sd[0]
        ci_95_sdd1 = ci_95_sd[1]

    f_a = interp1d(actual_time[0:len(mm_sdd)], mm_sdd, bounds_error=False)
    reg_a = f_a(actual_time)
    f_aci0 = interp1d(actual_time[0:len(mm_sdd)], ci_95_sdd0, bounds_error=False)
    reg_aci0 = f_aci0(actual_time)
    f_aci1 = interp1d(actual_time[0:len(mm_sdd)], ci_95_sdd1, bounds_error=False)
    reg_aci1 = f_aci1(actual_time)
    test_trials_mmeans.append(reg_a)
    test_trials_ci_95s.append([reg_aci0,reg_aci1])


ma_exp_plot2.plot_medang_dur_sd(test_trials_mmeans, test_trials_ci_95s, exp_s)
#ma_exp_plot2.plot_medang_dur_sd_control(test_trials_mmeans, test_trials_ci_95s, exp_s)


# ma_exp_plot.plot_all_ma081621_test_trials_sd(all_nosaccs_velos, reg_t, test_trials_list, exp_start_frames, exp_end_frames, all_pattern_velos,
#                                              datapaths, test_trials_mmeans, test_trials_ci_95s, 0, i)



# #normalized test_trials plot_decay_traces
# test_n_trials_mmeans = []
# test_n_trials_ci_95s = []
# for i in range(len(test_trials_list)):
#     test_trials = test_trials_list[i]
#
#     mm_sd, ci_95_sd = ma_analysis.get_test_trial_n_median_ma081621(all_nosaccs_velos, test_trials, exp_start_frames, exp_end_frames)
#     test_trials_ci_95s.append(ci_95_sd)
#     test_trials_mmeans.append(mm_sd)

###############################################################################
#pool data from both directions and get mean and ci

test_trials_list_p = []
for i in range(len(all_posgain_bias_trials)):
    posneg_bias_trials = [all_posgain_bias_trials[i][0], all_posgain_bias_trials[i][1],
                          all_neggain_bias_trials[i][0], all_neggain_bias_trials[i][1]]
    test_trials_list_p.append(posneg_bias_trials)

all_nosaccs_velos_absv = copy.deepcopy(all_nosaccs_velos)
for i in range(len(all_nosaccs_velos)):
    for j in range(len(all_posgain_bias_trials[i])):
        trial = all_posgain_bias_trials[i][j]
        s = exp_start_frames[trial]
        e = exp_end_frames[trial+1]
        all_nosaccs_velos_absv[i][s:e] = all_nosaccs_velos_absv[i][s:e]*-1

mm_p, ci_95_p = ma_analysis.get_test_trial_median_ma081621(all_nosaccs_velos_absv, test_trials_list_p, exp_start_frames, exp_end_frames)

# #plot all traces
# ma_exp_plot.plot_all_ma081621_test_trials(all_nosaccs_velos_absv, reg_t, test_trials_list_p, exp_start_frames, exp_end_frames, all_pattern_velos,
#                                            datapaths, mm_p, 0, i)
#
# #plot median, 95 ci
# ma_exp_plot.plot_all_ma081621_test_trials_sd(all_nosaccs_velos_absv, reg_t, [test_trials_list_p], exp_start_frames, exp_end_frames, all_pattern_velos,
#                                              datapaths, [mm_p], [ci_95_p], 0, i)

###############################################################################
#get and plot decay times

decayt_posgain_30 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_posgain_bias_trials, exp_start_frames, exp_end_frames, 30, fps)
decayt_posgain_50 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_posgain_bias_trials, exp_start_frames, exp_end_frames, 50, fps)

decayt_neggain_30 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_neggain_bias_trials, exp_start_frames, exp_end_frames, 30, fps)
decayt_neggain_50 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_neggain_bias_trials, exp_start_frames, exp_end_frames, 50, fps)

decayt_static_30 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_static_bias_trials, exp_start_frames, exp_end_frames, 30, fps)
decayt_static_50 = ma_analysis.get_all_decay_times(all_nosaccs_velos_absv, all_static_bias_trials, exp_start_frames, exp_end_frames, 50, fps)

test_trials_decayt_30 = [decayt_posgain_30, decayt_neggain_30, decayt_static_30]
test_trials_decayt_50 = [decayt_posgain_50, decayt_neggain_50, decayt_static_50]
###############################################################################
if save_data==1:
    #save processed data to dictionary
    all_processed_data = {}
    all_processed_data["experiment_name"] = exp
    all_processed_data["type_each_trial"] = type_each_trial
    all_processed_data["datapaths"] = datapaths
    all_processed_data["all_angles"] = all_angles
    all_processed_data["all_velos"] = all_velos
    all_processed_data["exp_start_frames"] = exp_start_frames
    all_processed_data["exp_end_frames"] = exp_end_frames
    all_processed_data["all_nosaccs_velos"] = all_nosaccs_velos
    all_processed_data["reg_t"] = reg_t
    all_processed_data["test_trials_list"] = test_trials_list
    all_processed_data["test_trials_mmeans"] = test_trials_mmeans
    all_processed_data["test_trials_ci_95s"] = test_trials_ci_95s
    all_processed_data["test_trials_decayt_30"] = test_trials_decayt_30
    all_processed_data["test_trials_decayt_50"] = test_trials_decayt_50
    all_processed_data["test_trials_pool_mmean"] = mm_p
    all_processed_data["test_trials_pool_ci_95"] = ci_95_p

    #save dictionary to a pickel
    pickle.dump((all_processed_data), open(pik_filename, "wb" ) )


    print(str(exp)+' processed data saved to pickel')

else:
    pass

###############################################################################
###############################################################################
###############################################################################
#dark periods analysis
# type_each_trial
# all_posgain_bias_trials
# all_nosaccs_velos

#zeroing_trials = [1, 6, 11, 16, 21, 26]

# plt.figure()
# ccs = []
# for i in range(len(all_nosaccs_velos)):
#     #for j in range(len(zeroing_trials)):
#     for j in range(len(all_posgain_bias_trials[i])):
#         try:
#             #trial = zeroing_trials[j]
#             trial = all_posgain_bias_trials[i][j]
#             s1 = exp_start_frames[trial+1]
#             e1 = exp_end_frames[trial+1]
#             s2 = exp_start_frames[trial+3]
#             e2 = exp_end_frames[trial+3]
#
#             #plt.plot(reg_t[0:len(all_nosaccs_velos[i][s1:e1])], (all_nosaccs_velos[i][s1:e1]), c = 'blue')
#             plt.plot(reg_t[0:len(all_nosaccs_velos[i][s2:e2])], (all_nosaccs_velos[i][s2:e2]), c = 'gray')
#
#         except:
#             pass
#
# plt.show()




# plt.figure()
# mm0=np.median(ccs,axis=0)
# plt.plot(mm0, c = 'k')
#
# upperHAF = np.nanmean(ccs, axis=0) + np.std(ccs, axis=0)
# lowerHAF = np.nanmean(ccs, axis=0) - np.std(ccs, axis=0)
# plt.fill_between(np.arange(1,len(upperHAF)+1,1), lowerHAF, upperHAF, color='k', zorder=1,
#                  edgecolor='none', alpha=.2)
#
# plt.show()

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
