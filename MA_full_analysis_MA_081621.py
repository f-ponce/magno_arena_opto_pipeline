import warnings
import numpy as np
import pickle
import random
import matplotlib.pyplot as plt
from collections import Counter
from numpy import linspace
from matplotlib import cm
import seaborn as sns

import basic_analysis_functions.MA_basic_analysis_functions as ma_analysis
import plot_functions.MA_plot_functions_basic_analysis as ma_plot
import plot_functions.MA_plot_functions_MA_081621 as ma_exp_plot
import import_functions.MA_filter_functions as ma_filter
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_run_sac_ID_ivo as ma_runsacc
import basic_analysis_functions.MA_sac_ID_functions as ma_saccf
###############################################################################

timepoints = [3, 10, 30, 60, 90, 180]
colors = ['purple', 'darkslateblue', 'cadetblue', 'darkturquoise', 'seagreen', 'yellow']

processed_dataDir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_processed_data/'

all_processed_data = []
for i in range(len(timepoints)):
    filename = 'MA_081621_'+str(timepoints[i])+'s_processed.p'
    pickelDir = processed_dataDir+filename
    pickle_file = open(pickelDir, "rb")
    processed_data = ma_analysis.open_pickel(pickle_file)
    all_processed_data.append(processed_data)

###############################################################################

#plot all the decay traces
#ma_exp_plot.plot_exp_decay_traces(all_processed_data, colors)
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
#     plt.scatter(randnums1, test1, color = colors[i])
#     plt.scatter(randnums2, test2, color = colors[i])
# plt.show()

fig, axs = plt.subplots(figsize=(7,3.5), facecolor='w', edgecolor='k')
sns.set_style("ticks")
for i in range(len(all_processed_data)):

    test1 = all_processed_data[i]["test_trials_decayt_30"][0]
    test2 = all_processed_data[i]["test_trials_decayt_30"][1]
    mtest1 = np.median(test1)
    mtest2 = np.median(test2)
    mtest = np.median([mtest1, mtest2])
    randnums1 = np.random.uniform(i, i-0.3, size=(len(test1)))
    randnums2 = np.random.uniform(i, i-0.3, size=(len(test2)))

    plt.scatter(randnums1, test1, color = colors[i])
    plt.scatter(randnums2, test2, color = colors[i])
    # plt.scatter(i, mtest1, color = 'r')
    # plt.scatter(i+0.2, mtest2, color = 'b')
    plt.scatter(i, mtest, color = 'k')

    #################
    #plot line at 0
    axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
    #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

    #plot limits
    axs.set_ylim([-100, 450])
    #axs.set_xlim([t[0], t[-1]])
    #plot labels
    axs.set_xlabel('time (s)',fontsize=15)
    axs.set_ylabel('angular velocity (deg/s)',fontsize=15)

    #plot format
    sns.despine()
    axs.tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs.spines['left'].set_linewidth(2)
    axs.spines['bottom'].set_linewidth(2)

    axs.tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs.spines['left'].set_linewidth(2)
    axs.spines['bottom'].set_linewidth(2)
    axs.yaxis.set_tick_params(labelsize=15)
    axs.xaxis.set_tick_params(labelsize=15)
fig.tight_layout()
plt.show()
