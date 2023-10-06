import warnings
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy
from matplotlib import cm
from numpy import linspace
import basic_analysis_functions.MA_sac_ID_functions as ma_saccf
###############################################################################
###############################################################################
# MA_081621
###############################################################################

def plot_ma081621_test_trials(velos, velos1, ovelos, reg_t, test_trials, start_frames, end_frames, pattern_velos, datapath, savemyfig, flynumber):
    savefigdir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/MA_081621_10s_figs/'
    exp = 'MA_081621_10s'

    for i in range(len(test_trials)):
        fig, axs = plt.subplots(figsize=(7,3.6), facecolor='w', edgecolor='k')
        sns.set_style("ticks")

        trial = test_trials[i]

        #get start/end frames of trials to plot
        s = start_frames[trial - 2]
        e = end_frames[trial + 1]

        s_bt = start_frames[trial - 1]
        e_bt = end_frames[trial - 1]

        s_t = start_frames[trial]
        e_t = end_frames[trial]

        s_at = start_frames[trial + 1]
        e_at = end_frames[trial + 1]

        #get the velo and time trace
        v = velos[s:e]
        v1 = velos1[s:e]
        ov = ovelos[s:e]
        t = reg_t[s:e] - reg_t[s]

        pv_t = pattern_velos[trial]
        pv_bt = pattern_velos[trial-1]

        #get the trial start/end times
        t_t = reg_t[s_t:e_t] - reg_t[s]
        t_bt = reg_t[s_bt:e_bt] - reg_t[s]
        t_at = reg_t[s_at:e_at] - reg_t[s]

        fly = 'fly_'+str(flynumber)+'_'+str(datapath[-24:-5])

        #plot traces
        # axs.plot(t,ov+250, c='g', linewidth=1)
        axs.plot(t,v, c='k', linewidth=1.5, zorder=0)
        #axs.plot(t,v1, c='r', linewidth=1, zorder=2)

        #plot formatting

        #plot background color
        axs.axvspan(t_bt[0], t_bt[-1], facecolor='darkseagreen', alpha=0.2)
        axs.axvspan(t_t[0], t_t[-1], facecolor='darkseagreen', alpha=0.4)

        #plot pattern angular velocities
        axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

        #plot line at 0
        axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
        #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

        #plot limits
        axs.set_ylim([-250, 250])
        axs.set_xlim([t[0], t[-1]])


        #plot labels
        axs.set_xlabel('time (s)',fontsize=20)
        axs.set_ylabel('angular velocity (deg/s)',fontsize=20)
        axs.set_title(str(fly),fontsize=20)

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
        axs.yaxis.set_tick_params(labelsize=20)
        axs.xaxis.set_tick_params(labelsize=20)


        if savemyfig == 1:
            fig.savefig(savefigdir+exp+'_trialdata_'+str(fly)+'_'+str(i)+'.pdf')  #bbox_inches='tight'
        else:
            pass
    #fig.tight_layout()
    plt.show()

###############################################################################

def plot_ma081621_test_trials_saccs(velos, velos1, ovelos, reg_t, test_trials, start_frames, end_frames, saccs_l_idx, saccs_r_idx, pattern_velos, datapath, savemyfig, flynumber):
    # savefigdir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/MA_081621_10s_figs/'
    # exp = 'MA_081621_10s'

    for i in range(len(test_trials)):
        fig, (ax1, ax2) = plt.subplots(2, figsize=(14,7), facecolor='w', edgecolor='k', gridspec_kw={'height_ratios': [1, 7]}, sharex=True)

        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        sns.set_style("ticks")

        trial = test_trials[i]

        #get start/end frames of trials to plot
        s = start_frames[trial - 2]
        e = end_frames[trial + 1]

        s_bt = start_frames[trial - 1]
        e_bt = end_frames[trial - 1]

        s_t = start_frames[trial]
        e_t = end_frames[trial]

        s_at = start_frames[trial + 1]
        e_at = end_frames[trial + 1]

        #get the velo and time trace
        v = velos[s:e]
        v1 = velos1[s:e]
        ov = ovelos[s:e]
        t = reg_t[s:e] - reg_t[s]

        pv_t = pattern_velos[trial]
        pv_bt = pattern_velos[trial-1]

        #get the trial start/end times
        t_t = reg_t[s_t:e_t] - reg_t[s]
        t_bt = reg_t[s_bt:e_bt] - reg_t[s]
        t_at = reg_t[s_at:e_at] - reg_t[s]

        #subplot 2 (bottom)
        fly = 'fly_'+str(flynumber)+'_'+str(datapath[-24:-5])

        #plot traces
        ax2.plot(t,ov+250, c='g', linewidth=1)
        ax2.plot(t,v, c='k', linewidth=1.5, zorder=0)
        ax2.plot(t,v1, c='r', linewidth=1, zorder=2)

        #subplot 1 (top)
        sl = np.argwhere(~np.isnan(saccs_l_idx[s:e]))
        sr = np.argwhere(~np.isnan(saccs_r_idx[s:e]))

        for k in range(len(sr)):
            ax1.eventplot(reg_t[sr[k]], lineoffsets=10, color = 'blue', linelengths=10, linewidths=2, alpha = 0.7)
        for k in range(len(sl)):
            ax1.eventplot(reg_t[sl[k]], lineoffsets=1, color = 'red', linelengths=10, linewidths=2, alpha = 0.7)

        #plot formatting (subplot2)
        #plot background color
        ax2.axvspan(t_bt[0], t_bt[-1], facecolor='darkseagreen', alpha=0.2, ymin=-1, ymax=0.9)
        ax2.axvspan(t_t[0], t_t[-1], facecolor='darkseagreen', alpha=0.4, ymin=-1, ymax=0.9)

        #plot pattern angular velocities
        ax2.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
        ax2.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
        ax2.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

        #plot line at 0
        ax2.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
        #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

        # #plot limits
        # axs.set_ylim([-500, 500])
        # axs.set_xlim([t[0], t[-1]])

        #plot labels
        ax2.set_xlabel('time (s)',fontsize=20)
        ax2.set_ylabel('angular velocity (deg/s)',fontsize=20)

        #plot format
        sns.despine()
        ax2.tick_params(direction='in', length=8, width=2)
        sns.despine(offset=10, trim=False);
        ax2.spines['left'].set_linewidth(2)
        ax2.spines['bottom'].set_linewidth(2)

        ax2.tick_params(direction='in', length=8, width=2)
        sns.despine(offset=10, trim=False);
        ax2.spines['left'].set_linewidth(2)
        ax2.spines['bottom'].set_linewidth(2)
        ax2.yaxis.set_tick_params(labelsize=20)
        ax2.xaxis.set_tick_params(labelsize=20)

        #plot format (subplot2)
        ax1.set_title(str(fly),fontsize=20)
        ax1.axis('off')
        # ax1.yaxis.set_tick_params(left=False)
        # ax1.xaxis.set_tick_params(bottom=False)
        # ax1.set(xlabel=None)
        # ax1.set(ylabel=None)
        # ax1.set(yticklabels=[])
        # #ax1.set(xticklabels=[])

        if savemyfig == 1:
            fig.savefig(savefigdir+exp+'_trialdata_'+str(fly)+'_'+str(i)+'.pdf')  #bbox_inches='tight'
        else:
            pass
    #fig.tight_layout()
    #plt.show()
###############################################################################

def plot_all_ma081621_test_trials(all_velos, reg_t, all_test_trials, exp_start_frames, exp_end_frames, all_pattern_velos, datapaths, trial_median, savemyfig, trialc):

    # savefigdir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/MA_081621_10s_figs/'
    # exp = 'MA_081621_10s'

    fig, axs = plt.subplots(figsize=(7,3.5), facecolor='w', edgecolor='k')
    sns.set_style("ticks")
    for i in range(len(all_velos)):
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial - 2]
            e = exp_end_frames[trial + 1]

            s_bt = exp_start_frames[trial - 1]
            e_bt = exp_end_frames[trial - 1]

            s_t = exp_start_frames[trial]
            e_t = exp_end_frames[trial]

            s_at = exp_start_frames[trial + 1]
            e_at = exp_end_frames[trial + 1]

            #get the velo and time trace
            v = all_velos[i][s:e]
            #ov = all_ovelfilt_velos[i][s:e]
            t = reg_t[s:e] - reg_t[s]

            #trial median
            mm = trial_median
            tmm = reg_t[0:len(mm)]

            pv_t = all_pattern_velos[i][trial]
            pv_bt = all_pattern_velos[i][trial-1]

            #get the trial start/end times
            t_t = reg_t[s_t:e_t] - reg_t[s]
            t_bt = reg_t[s_bt:e_bt] - reg_t[s]
            t_at = reg_t[s_at:e_at] - reg_t[s]

            fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])

            #plot traces
            axs.plot(t,v, c='gray', linewidth=1)
            #axs.plot(t,ov, c='b', linewidth=1)
            axs.plot(tmm,mm, c='k', linewidth=1, zorder = 5)

    #plot formatting

    #plot pattern angular velocities
    axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
    axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
    axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

    #plot background color
    axs.axvspan(t_bt[0], t_bt[-1], facecolor='darkseagreen', alpha=0.2)
    axs.axvspan(t_t[0], t_t[-1], facecolor='darkseagreen', alpha=0.4)

    #plot pattern angular velocities
    axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
    axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
    axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

    #plot line at 0
    axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
    #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

    #plot limits
    axs.set_ylim([-500, 500])
    axs.set_xlim([t[0], t[-1]])


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

    if savemyfig == 1:
        fig.savefig(savefigdir+exp+'_trialdata_'+str(trialc)+'.pdf')  #bbox_inches='tight'
    else:
        pass

    fig.tight_layout()
    plt.show()

###############################################################################

def plot_all_ma081621_test_trials_perfly(all_velos, reg_t, all_test_trials, exp_start_frames, exp_end_frames, all_pattern_velos, datapaths, trial_median, savemyfig, trialc):

    # savefigdir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/MA_081621_10s_figs/'
    # exp = 'MA_081621_10s'

    for i in range(len(all_velos)):
        fig, axs = plt.subplots(figsize=(7,3.5), facecolor='w', edgecolor='k')
        sns.set_style("ticks")
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial - 2]
            e = exp_end_frames[trial + 1]

            s_bt = exp_start_frames[trial - 1]
            e_bt = exp_end_frames[trial - 1]

            s_t = exp_start_frames[trial]
            e_t = exp_end_frames[trial]

            s_at = exp_start_frames[trial + 1]
            e_at = exp_end_frames[trial + 1]

            #get the velo and time trace
            v = all_velos[i][s:e]
            #ov = all_ovelfilt_velos[i][s:e]
            t = reg_t[s:e] - reg_t[s]

            #trial median
            mm = trial_median
            tmm = reg_t[0:len(mm)]

            pv_t = all_pattern_velos[i][trial]
            pv_bt = all_pattern_velos[i][trial-1]

            #get the trial start/end times
            t_t = reg_t[s_t:e_t] - reg_t[s]
            t_bt = reg_t[s_bt:e_bt] - reg_t[s]
            t_at = reg_t[s_at:e_at] - reg_t[s]

            fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])

            #plot traces
            axs.plot(t,v, c='gray', linewidth=1)
            #axs.plot(t,ov, c='b', linewidth=1)
            axs.plot(tmm,mm, c='k', linewidth=1, zorder = 5)

        #plot formatting

        #plot pattern angular velocities
        axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

        #plot background color
        axs.axvspan(t_bt[0], t_bt[-1], facecolor='darkseagreen', alpha=0.2)
        axs.axvspan(t_t[0], t_t[-1], facecolor='darkseagreen', alpha=0.4)

        #plot pattern angular velocities
        axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
        axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

        #plot line at 0
        axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
        #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

        #plot limits
        axs.set_ylim([-500, 500])
        axs.set_xlim([t[0], t[-1]])


        #plot labels
        axs.set_xlabel('time (s)',fontsize=15)
        axs.set_ylabel('angular velocity (deg/s)',fontsize=15)
        axs.set_title(str(fly),fontsize=20)

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

        if savemyfig == 1:
            fig.savefig(savefigdir+exp+'_trialdata_'+str(trialc)+'.pdf')  #bbox_inches='tight'
        else:
            pass

        fig.tight_layout()
        plt.show()

###############################################################################

def plot_all_ma081621_test_trials_sd(all_velos, reg_t, test_trials_list, exp_start_frames, exp_end_frames, all_pattern_velos, datapaths, trial_median, ci_95, savemyfig, trialc):

    # savefigdir = '/Users/fponce/Documents/magno_arena_opto/analysis_MA_2022/pickled_data/MA_081621_180s_figs/'
    # exp = 'MA_081621_180s'
    #10-7.36
    fig, axs = plt.subplots(figsize=(7.36,3.5), facecolor='w', edgecolor='k')
    sns.set_style("ticks")
    for i in range(len(all_velos)):
        for k in range(len(test_trials_list)):
            for j in range(len(test_trials_list[k][i])):

                trial = test_trials_list[k][i][j]

                #get start/end frames of trials to plot
                s = exp_start_frames[trial - 2]
                e = exp_end_frames[trial + 1]

                s_bt = exp_start_frames[trial - 1]
                e_bt = exp_end_frames[trial - 1]

                s_t = exp_start_frames[trial]
                e_t = exp_end_frames[trial]

                s_at = exp_start_frames[trial + 1]
                e_at = exp_end_frames[trial + 1]

                #get the velo and time trace
                v = all_velos[i][s:e]
                #ov = all_ovelfilt_velos[i][s:e]
                t = reg_t[s:e] - reg_t[s]

                #trial median
                mm = trial_median[k]
                tmm = reg_t[0:len(mm)]

                wn = 0.05
                # mean_y = mm
                b, a = scipy.signal.butter(3, wn)
                # #mean_y = scipy.signal.filtfilt(b,a,mean_y)
                # third = int(0.25*len(mean_y))
                # baseline_shift = np.mean(mean_y[0:1800])#[0:third] )
                # mean_y -= baseline_shift
                norm_factor = 1#np.max(mean_y)

                lo_filtered = scipy.signal.filtfilt(b,a,ci_95[k][0])
                hi_filtered = scipy.signal.filtfilt(b,a,ci_95[k][1])

                lowerHAF = (lo_filtered-0.)/norm_factor#(lo_filtered-baseline_shift)/norm_factor
                upperHAF = (hi_filtered-0.)/norm_factor#(hi_filtered-baseline_shift)/norm_factor

                #get the pattern velo
                pv_t = all_pattern_velos[i][trial]
                pv_bt = all_pattern_velos[i][trial-1]

                #get the trial start/end times
                t_t = reg_t[s_t:e_t] - reg_t[s]
                t_bt = reg_t[s_bt:e_bt] - reg_t[s]
                t_at = reg_t[s_at:e_at] - reg_t[s]

                fly = 'fly_'+str(i)+'_'+str(datapaths[i][-24:-5])

                #plot traces
                #axs.plot(t,v, c='gray', linewidth=1)
                #axs.plot(t,ov, c='b', linewidth=1)
                axs.plot(tmm,mm, c='k', linewidth=1, zorder = 5)
                axs.fill_between(tmm, lowerHAF, upperHAF, facecolor='gray', alpha=0.2, edgecolor='none')

                if i ==0:
                    #plot pattern angular velocities
                    axs.plot((t_bt[0], t_bt[-1]), (pv_bt, pv_bt), 'lime', linewidth= 2)
                    axs.plot((t_t[0], t_t[-1]), (pv_t, pv_t), 'lime', linewidth= 2)
                    axs.plot((t_t[0], t_t[0]), (pv_bt, pv_t), 'lime', linewidth= 2)

    #plot formatting

    #plot background color
    axs.axvspan(t_bt[0], t_bt[-1], facecolor='darkseagreen', alpha=0.2)
    axs.axvspan(t_t[0], t_t[-1], facecolor='darkseagreen', alpha=0.4)

    #plot line at 0
    axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
    #plt.axhline(y=fly_setpoint, color='blue', linestyle='--', linewidth= 0.9)

    #plot limits
    axs.set_ylim([-250, 250])
    axs.set_xlim([t[0], t[-1]])


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

    if savemyfig == 1:
        fig.savefig(savefigdir+exp+'_trialdata_'+str(trialc)+'.pdf')  #bbox_inches='tight'
    else:
        pass

    fig.tight_layout()
    plt.show()

###############################################################################
###############################################################################
###############################################################################

def plot_exp_decay_traces(all_processed_data, colors):

    # #defining color map
    # start = 0.0
    # stop = 1.0
    # number_of_lines = len(all_processed_data)
    # cm_subsection = linspace(start, stop, number_of_lines)
    # colors = [ cm.viridis(x) for x in cm_subsection ]

    fig, axs = plt.subplots(1, 2, figsize=(6,5), facecolor='w', edgecolor='k', gridspec_kw={'width_ratios': [50, 1]})
    fig.subplots_adjust(wspace=0.02)
    sns.set_style("ticks")

    for i in range(len(all_processed_data)):
        exp_mm = all_processed_data[i]["test_trials_pool_mmean"][-1840:-1]
        t = all_processed_data[i]["reg_t"][0:len(exp_mm)]

        #plot 95 ci
        wn = 0.05
        mean_y = exp_mm
        b, a = scipy.signal.butter(3, wn)
        mean_y = scipy.signal.filtfilt(b,a,mean_y)
        third = int(0.25*len(mean_y))
        baseline_shift = np.mean( mean_y[0:180])#[0:third] )
        mean_y -= baseline_shift
        norm_factor = 1#np.max(mean_y)

        b, a = scipy.signal.butter(3, wn)

        ci_951 = all_processed_data[i]["test_trials_pool_ci_95"][0][-1840:-1]
        ci_952 = all_processed_data[i]["test_trials_pool_ci_95"][1][-1840:-1]
        lo_filtered = scipy.signal.filtfilt(b,a,ci_951)
        hi_filtered = scipy.signal.filtfilt(b,a,ci_952)

        lowerHAF = (lo_filtered-0.)/norm_factor#(lo_filtered-baseline_shift)/norm_factor
        upperHAF = (hi_filtered-0.)/norm_factor#(hi_filtered-baseline_shift)/norm_factor

        #plot traces
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm, 17*1), c = colors[i])
        #axs[0].fill_between(t, lowerHAF, upperHAF, facecolor='gray', alpha=0.3, edgecolor='none')

    #plot pattern angular pattern_angular_velocity
    axs[0].plot((t[0], t[40]), (180, 180), 'k', linewidth= 2)

    #plot colorbar
    # gradient = np.linspace(0, 1, 256)
    # gradient = np.vstack((gradient)*-1)
    #axs.imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
    #axs[1].imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))

    #plot formatting

    #plot limits
    axs[0].set_ylim([0, 250])
    axs[0].set_xlim([t[0], t[-1]])

    #plot labels
    axs[0].set_xlabel('time (s)',fontsize=20)
    axs[0].set_ylabel('angular velocity (deg/s)',fontsize=20)

    #plot format
    sns.despine()
    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)

    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].yaxis.set_tick_params(labelsize=20)
    axs[0].xaxis.set_tick_params(labelsize=20)

    axs[1].set_axis_off()

    fig.tight_layout()
    plt.show()
###############################################################################

def plot_exp_decay_2traces(all_processed_data):

    #defining color map
    start = 0.0
    stop = 1.0
    number_of_lines = len(all_processed_data)
    cm_subsection = linspace(start, stop, number_of_lines)
    colors = [ cm.viridis(x) for x in cm_subsection ]

    fig, axs = plt.subplots(1, 2, figsize=(6,5), facecolor='w', edgecolor='k', gridspec_kw={'width_ratios': [10000, 1]})
    fig.subplots_adjust(wspace=0.02)
    sns.set_style("ticks")

    for i in range(len(all_processed_data)):
        exp_mm1 = all_processed_data[i]["test_trials_mmeans"][0][-1840:-1]
        exp_mm2 = all_processed_data[i]["test_trials_mmeans"][1][-1840:-1]

        t = all_processed_data[i]["reg_t"][0:len(exp_mm1)]
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm1, 17*1)*-1, c = colors[i])
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm2, 17*1), c = colors[i])

    #plot pattern angular pattern_angular_velocity
    axs[0].plot((t[0], t[40]), (180, 180), 'k', linewidth= 2)

    #plot colorbar
    # gradient = np.linspace(0, 1, 256)
    # gradient = np.vstack((gradient)*-1)

    #axs.imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
    #axs[1].imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
    #plot formatting

    #plot limits
    axs[0].set_ylim([0, 250])
    axs[0].set_xlim([t[0], t[-1]])

    #plot labels
    axs[0].set_xlabel('time (s)',fontsize=20)
    axs[0].set_ylabel('angular velocity (deg/s)',fontsize=20)

    #plot format
    sns.despine()
    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)

    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].yaxis.set_tick_params(labelsize=20)
    axs[0].xaxis.set_tick_params(labelsize=20)

    # axs[1].set_axis_off()

    fig.tight_layout()
    plt.show()

###############################################################################

def plot_exp_decay_per_tp(all_processed_data):

    #defining color map
    start = 0.0
    stop = 1.0
    number_of_lines = len(all_processed_data)
    cm_subsection = linspace(start, stop, number_of_lines)
    colors = [ cm.viridis(x) for x in cm_subsection ]

    for i in range(len(all_processed_data)):

        fig, axs = plt.subplots(1, 2, figsize=(6,5), facecolor='w', edgecolor='k', gridspec_kw={'width_ratios': [1000, 1]})
        fig.subplots_adjust(wspace=0.02)
        sns.set_style("ticks")

        exp_mm1 = all_processed_data[i]["test_trials_mmeans"][0][-1840:-1]
        exp_mm2 = all_processed_data[i]["test_trials_mmeans"][1][-1850:-1]
        exp_mm = all_processed_data[i]["test_trials_pool_mmean"][-1840:-1]

        exp_ci951 = all_processed_data[i]["test_trials_ci_95s"][0][-1840:-1]
        exp_ci952 = all_processed_data[i]["test_trials_ci_95s"][1][-1840:-1]

        t = all_processed_data[i]["reg_t"][0:len(exp_mm1)]
        t2 = all_processed_data[i]["reg_t"][0:len(exp_mm2)]
        #plot traces
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm1, 17*1)*-1, c = 'r')
        axs[0].plot(t2,ma_saccf.overfilter_ang_vel(exp_mm2, 17*1), c = 'b')
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm, 17*1), c = 'k')

        #plot pattern angular pattern_angular_velocity
        axs[0].plot((t[0], t[40]), (180, 180), 'k', linewidth= 2)

        #plot colorbar
        gradient = np.linspace(0, 1, 256)
        gradient = np.vstack((gradient)*-1)

        #axs.imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
        #axs[1].imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
        #plot formatting

        #plot limits
        axs[0].set_ylim([0, 250])
        axs[0].set_xlim([t[0], t[-1]])

        #plot labels
        exp = all_processed_data[i]['experiment_name']
        axs[0].set_title(exp,fontsize=20)
        axs[0].set_xlabel('time (s)',fontsize=20)
        axs[0].set_ylabel('angular velocity (deg/s)',fontsize=20)

        #plot format
        sns.despine()
        axs[0].tick_params(direction='in', length=8, width=2)
        sns.despine(offset=10, trim=False);
        axs[0].spines['left'].set_linewidth(2)
        axs[0].spines['bottom'].set_linewidth(2)

        axs[0].tick_params(direction='in', length=8, width=2)
        sns.despine(offset=10, trim=False);
        axs[0].spines['left'].set_linewidth(2)
        axs[0].spines['bottom'].set_linewidth(2)
        axs[0].yaxis.set_tick_params(labelsize=20)
        axs[0].xaxis.set_tick_params(labelsize=20)

        axs[1].set_axis_off()

        fig.tight_layout()
        plt.show()
###############################################################################

###############################################################################

def plot_decay_traces_tp(all_processed_data, colors):

    # #defining color map
    # start = 0.0
    # stop = 1.0
    # number_of_lines = len(all_processed_data)
    # cm_subsection = linspace(start, stop, number_of_lines)
    # colors = [ cm.viridis(x) for x in cm_subsection ]

    fig, axs = plt.subplots(1, 2, figsize=(6,5), facecolor='w', edgecolor='k', gridspec_kw={'width_ratios': [50, 1]})
    fig.subplots_adjust(wspace=0.02)
    sns.set_style("ticks")

    for i in range(len(all_processed_data)):
        exp_mm = all_processed_data[i]["test_trials_pool_mmean"][-1840:-1]
        t = all_processed_data[i]["reg_t"][0:len(exp_mm)]

        #plot traces
        axs[0].plot(t,ma_saccf.overfilter_ang_vel(exp_mm, 17*1), c = colors[i])
        axs[0].fill_between(t, lowerHAF, upperHAF, facecolor='gray', alpha=0.3, edgecolor='none')

    #plot pattern angular pattern_angular_velocity
    axs[0].plot((t[0], t[40]), (180, 180), 'k', linewidth= 2)

    #plot colorbar
    # gradient = np.linspace(0, 1, 256)
    # gradient = np.vstack((gradient)*-1)
    #axs.imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))
    #axs[1].imshow(gradient, aspect='auto', cmap=plt.get_cmap('viridis'))

    #plot formatting

    #plot limits
    axs[0].set_ylim([0, 250])
    axs[0].set_xlim([t[0], t[-1]])

    #plot labels
    axs[0].set_xlabel('time (s)',fontsize=20)
    axs[0].set_ylabel('angular velocity (deg/s)',fontsize=20)

    #plot format
    sns.despine()
    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)

    axs[0].tick_params(direction='in', length=8, width=2)
    sns.despine(offset=10, trim=False);
    axs[0].spines['left'].set_linewidth(2)
    axs[0].spines['bottom'].set_linewidth(2)
    axs[0].yaxis.set_tick_params(labelsize=20)
    axs[0].xaxis.set_tick_params(labelsize=20)

    axs[1].set_axis_off()

    fig.tight_layout()
    plt.show()
