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

exp = 'MA_031021'

filename = exp+'_preprocessed.p'

dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/'
pickelDir = dataDir+filename

pickle_file = open(pickelDir, "rb")
all_processed_rawdata = []
all_processed_rawdata = ma_analysis.open_pickel(pickle_file)

fps = 30

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
# ###############################################################################

# ###############################################################################
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
print(exp_start_frames)
print(exp_end_frames)
# ###############################################################################
#find and remove saccades

all_nosaccs_velos_2, all_saccs_id, all_saccs_l_idx, all_saccs_r_idx, all_overfilt_velos = ma_runsacc.get_all_nosacc_velos_loop(all_velos, reg_t)

all_flipped_vels = []
for i in range(len(all_nosaccs_velos_2)):
    if all_gains[i][1]==-99:
        flipped_vel = all_nosaccs_velos_2[i]*-1
    else:
        flipped_vel = all_nosaccs_velos_2[i]
    all_flipped_vels.append(flipped_vel)

all_nosaccs_velos = ma_runsacc.get_all_interp_after_saccs(all_flipped_vels, reg_t)


med_angvel = np.nanmean(np.asarray(all_nosaccs_velos), axis=0)
ci_95s = ma_analysis.get_95_confidence_intervals(all_nosaccs_velos, iterations=100)

ma_exp_plot2.plot_medang_longdur(reg_t, med_angvel, ci_95s)
#



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
