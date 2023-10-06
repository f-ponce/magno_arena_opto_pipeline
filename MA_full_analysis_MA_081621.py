import warnings
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
from collections import Counter
from numpy import linspace
from matplotlib import cm
import seaborn as sns
import itertools

import basic_analysis_functions.MA_basic_analysis_functions as ma_analysis
import plot_functions.MA_plot_functions_basic_analysis as ma_plot
import plot_functions.MA_plot_functions_MA_081621 as ma_exp_plot
import plot_functions.plot_magno_funcs as ma_exp_plot2
import import_functions.MA_filter_functions as ma_filter
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_run_sac_ID_ivo as ma_runsacc
import basic_analysis_functions.MA_sac_ID_functions as ma_saccf
###############################################################################

timepoints = [3, 10, 30, 60, 90, 180]
timepoints = [40, 60, 80, 100]
timepoints = [3, 10, 30, 40, 60, 80, 90, 100, 180]

#colors = ['purple', 'darkslateblue', 'cadetblue', 'red', 'darkturquoise', 'darkturquoise', 'seagreen', 'seagreen', 'yellow']
colors = ['purple', 'darkslateblue', 'cadetblue', 'blue', 'red', 'red', 'magenta', 'green','k','yellow']
colors = ['purple', 'darkslateblue', 'cadetblue', 'blue', 'red', 'magenta', 'green','k','orange']
colors = ['#fdca26', '#fb9f3a', '#ed7953', '#d8576b', '#bd3786', '#9c179e', '#7201a8', '#46039f', '#0d0887']

experiments_names = ['MA_081621', 'MA_032122']
#experiments_names = ['MA_032122']
#experiments_names = ['MA_081621']
processed_dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_processed_data/'

# all_processed_data = []
# for i in range(len(timepoints)):
#     #filename = 'MA_081621_'+str(timepoints[i])+'s_processed.p'
#     filename = 'MA_032122_'+str(timepoints[i])+'s_processed.p'
#     pickelDir = processed_dataDir+filename
#     pickle_file = open(pickelDir, "rb")
#     processed_data = ma_analysis.open_pickel(pickle_file)
#     all_processed_data.append(processed_data)

all_processed_data = []
for i in range(len(timepoints)):
    for j in range(len(experiments_names)):
        filename_str = experiments_names[j]
        try:
            filename = filename_str+'_'+str(timepoints[i])+'s_processed.p'
            pickelDir = processed_dataDir+filename
            pickle_file = open(pickelDir, "rb")
            processed_data = ma_analysis.open_pickel(pickle_file)
            all_processed_data.append(processed_data)
        except:
            pass

        #filename = 'MA_081621_'+str(timepoints[i])+'s_processed.p'
        #filename = 'MA_032122_'+str(timepoints[i])+'s_processed.p'
        # pickelDir = processed_dataDir+filename
        # pickle_file = open(pickelDir, "rb")
        # processed_data = ma_analysis.open_pickel(pickle_file)
        # all_processed_data.append(processed_data)

###############################################################################

#plot all the decay traces

ma_exp_plot2.plot_exp_decay_traces(all_processed_data, colors)

#ma_exp_plot.plot_exp_decay_2traces(all_processed_data)
#ma_exp_plot.plot_exp_decay_per_tp(all_processed_data)


###############################################################################

#defining color map

# plt.figure()
# for i in range(len(all_processed_data)):
#     test1 = all_processed_data[i]["test_trials_decayt_50"][0]*-1
#     test2 = all_processed_data[i]["test_trials_decayt_50"][1]
#     randnums1 = np.random.uniform(i, i-0.3, size=(len(test1)))
#     randnums2 = np.random.uniform(i, i-0.3, size=(len(test2)))
#     plt.scatter(randnums1, test1, color = colors[i], alpha = 0.5)
#     plt.scatter(randnums2, test2, color = colors[i], alpha = 0.5)
# plt.show()
all_mtests = []
for i in range(len(all_processed_data)):
    test1 = all_processed_data[i]["test_trials_decayt_50"][0]
    test2 = all_processed_data[i]["test_trials_decayt_50"][1]
    mtests = []
    for j in range(len(test1)):
        mt = np.median([test1[j], test2[j]])
        mtests.append(mt)
    all_mtests.append(mtests)

###############################################################################
ts = [3, 10, 30, 40, 60, 80, 90, 100, 180]

ma_exp_plot2.plot_exp_decay_times(all_processed_data, colors, all_mtests, ts)


###############################################################################
#ts = [3, 10, 30, 40, 60, 80, 90, 100, 180]
# sns.set_style("ticks")
# fig, axs = plt.subplots(figsize=(7,4), facecolor='w', edgecolor='k')
#
# for i in range(len(all_processed_data)):
#
#     test1 = all_processed_data[i]["test_trials_decayt_50"][0]
#     test2 = all_processed_data[i]["test_trials_decayt_50"][1]
#     tests = all_mtests[i]
#
#     mtest1 = np.median(test1)
#     mtest2 = np.median(test2)
#     mtest = np.median([mtest1, mtest2])
#     tt = [test1, test2]
#     ftt =  list(itertools.chain(*tt))
#     # cc1 = ma_analysis.(ftt)
#     # cc2 = ma_analysis.(test2)
#     # cc11= (cc1[0]+cc2[0])/2
#     # cc22= (cc1[1]+cc2[1])/2
#
#
#
#     randnums = np.random.uniform(ts[i]-1, ts[i]+1, size=(len(tests)))
#     randnums1 = np.random.uniform(i, i-0.5, size=(len(test1)))
#     randnums2 = np.random.uniform(i, i-0.5, size=(len(test2)))
#
#     plt.scatter((randnums), tests, color = colors[i], alpha=0.6, s = 45)
#     # plt.scatter(randnums1, test1, color = colors[i], alpha=0.5)
#     # plt.scatter(randnums2, test2, color = colors[i], alpha=0.5)
#     # plt.scatter(i, mtest1, color = 'r')
#     # plt.scatter(i+0.2, mtest2, color = 'b')
#     plt.scatter((ts[i]), mtest, color = 'k', s = 45,marker='s')
#     #plt.errorbar((ts[i]), mtest, yerr = cc1,fmt='o',ecolor = colors[i],color='black')
#
#     #################
#     #plot line at 0
#     #axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#     axs.axhline(y=182.0, color='green', linestyle='--', linewidth= 1.5)
#     #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
#
#     #plot limits
#     axs.set_ylim([-0, 250])
#     axs.set_xlim([0, 200])
#
#     #plot labels
#     axs.set_xlabel('time (s)',fontsize=20)
#     axs.set_ylabel('angular velocity (deg/s)',fontsize=20)
#
#     #plot format
#     sns.despine()
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#     axs.yaxis.set_tick_params(labelsize=20)
#     axs.xaxis.set_tick_params(labelsize=20)
#
# fig.tight_layout()
# plt.show()

###############################################################################
# ts = [3, 10, 30, 40, 60, 63, 80, 90, 100, 180]
# fig, axs = plt.subplots(figsize=(7,3.5), facecolor='w', edgecolor='k')
# sns.set_style("ticks")
# for i in range(len(all_processed_data)):
#
#     #test1 = all_processed_data[i]["test_trials_decayt_50"][0]
#     test1 = all_processed_data[i]["test_trials_decayt_50"][2]
#
#     mtest1 = np.median(test1)
#     mtest = np.median([mtest1])
#
#     randnums = np.random.uniform(ts[i], ts[i]-0.3, size=(len(test1)))
#     randnums1 = np.random.uniform(i, i-0.3, size=(len(test1)))
#     randnums2 = np.random.uniform(i, i-0.3, size=(len(test2)))
#
#     plt.scatter((randnums), test1, color = colors[i], alpha=0.5)
#     # plt.scatter(randnums1, test1, color = colors[i], alpha=0.5)
#     # plt.scatter(randnums2, test2, color = colors[i], alpha=0.5)
#     # plt.scatter(i, mtest1, color = 'r')
#     # plt.scatter(i+0.2, mtest2, color = 'b')
#     plt.scatter((ts[i]), mtest, color = 'k')
#
#     #################
#     #plot line at 0
#     axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#     axs.axhline(y=182.0, color='green', linestyle='-', linewidth= 0.9)
#     #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
#
#     #plot limits
#     axs.set_ylim([-300, 300])
#     #axs.set_xlim([t[0], t[-1]])
#     #plot labels
#     axs.set_xlabel('time (s)',fontsize=15)
#     axs.set_ylabel('angular velocity (deg/s)',fontsize=15)
#
#     #plot format
#     sns.despine()
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#     axs.yaxis.set_tick_params(labelsize=15)
#     axs.xaxis.set_tick_params(labelsize=15)
#
# fig.tight_layout()
# plt.show()

# ###############################################################################
# fig, axs = plt.subplots(figsize=(7,3.5), facecolor='w', edgecolor='k')
# sns.set_style("ticks")
# for i in range(len(all_processed_data)):
#
#     test1 = all_processed_data[i]["test_trials_decayt_50"][0]
#     test2 = all_processed_data[i]["test_trials_decayt_50"][1]
#     tests = np.median()
#
#     mtest1 = np.median(test1)
#     mtest2 = np.median(test2)
#     mtest = np.median([mtest1, mtest2])
#     randnums1 = np.random.uniform(i, i-0.3, size=(len(test1)))
#     randnums2 = np.random.uniform(i, i-0.3, size=(len(test2)))
#
#     plt.scatter(randnums1, test1, color = colors[i], alpha=0.5)
#     plt.scatter(randnums2, test2, color = colors[i], alpha=0.5)
#     # plt.scatter(i, mtest1, color = 'r')
#     # plt.scatter(i+0.2, mtest2, color = 'b')
#     plt.scatter(i, mtest, color = 'k')
#
#     #################
#     #plot line at 0
#     axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
#     #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)
#
#     #plot limits
#     axs.set_ylim([-100, 450])
#     #axs.set_xlim([t[0], t[-1]])
#     #plot labels
#     axs.set_xlabel('time (s)',fontsize=15)
#     axs.set_ylabel('angular velocity (deg/s)',fontsize=15)
#
#     #plot format
#     sns.despine()
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#
#     axs.tick_params(direction='in', length=8, width=2)
#     sns.despine(offset=10, trim=False);
#     axs.spines['left'].set_linewidth(2)
#     axs.spines['bottom'].set_linewidth(2)
#     axs.yaxis.set_tick_params(labelsize=15)
#     axs.xaxis.set_tick_params(labelsize=15)
# fig.tight_layout()
# plt.show()
