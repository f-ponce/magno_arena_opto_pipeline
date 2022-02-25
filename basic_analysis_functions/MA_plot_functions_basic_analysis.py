import warnings
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
###############################################################################

def plot_ledpanels_angvelos (trial_start_frames, trial_end_frames, gains, g_t_av_f):
    sns.set_style("ticks")

    for i in range(len(trial_start_frames)):
        ang_vel_pattern = gains[i]*g_t_av_f #gain to ang vel conversion

        s = trial_start_frames[i]
        e = trial_end_frames[i]

        plt.axvspan(s, e, facecolor='darkseagreen', alpha=0.2)
        plt.plot((s,e), (ang_vel_pattern, ang_vel_pattern), 'lime', linewidth= 1)

def color_trial_background(trial_start_times, trial_end_times, type_trials):
    sns.set_style("ticks")

    background_color = {}
    background_color["dark"] = ['gray', 0.2]
    background_color["moving"] = ['darkseagreen', 0.5]
    background_color["static"] = ['darkseagreen', 0.1]

    for i in range(len(trial_start_times)):

        s = trial_start_times[i]
        e = trial_end_times[i]

        bc = background_color[type_trials[i]][0]
        bc_a = background_color[type_trials[i]][1]
        plt.axvspan(s, e, facecolor=bc, alpha=bc_a)

def plot_exp_angvelo(reg_t, angvelo, trial_start_times, trial_end_times, gains, g_t_av_f, type_trials, title):
    sns.set_style("ticks")

    fig, axs = plt.subplots(figsize=(14,7), facecolor='w', edgecolor='k')
    axs.plot(reg_t, angvelo, color='k', alpha = 1)

    #plot ledpanels angular velocity
    plot_ledpanels_angvelos (trial_start_times, trial_end_times, gains, g_t_av_f)

    #color background of trials
    color_trial_background(trial_start_times, trial_end_times, type_trials)

    #horizontal x=0 line
    #axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)

    #plot limits
    plt.ylim([-500, 500])
    plt.xlim([0, reg_t[-1]])

    #plot label
    axs.set_xlabel('time (s)',fontsize=20)
    axs.set_ylabel('angular velocity (deg/s)',fontsize=20)
    axs.set_title(str(title),fontsize=20)

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

    plt.show()

def plot_meanangvelo(reg_t, angvelo, mean_trial_start_times, mean_trial_end_times, gains, g_t_av_f, type_trials, title):
    sns.set_style("ticks")

    fig, axs = plt.subplots(figsize=(14,7), facecolor='w', edgecolor='k')
    axs.plot(reg_t, angvelo, color='k', alpha = 1)

    #plot ledpanels angular velocity
    plot_ledpanels_angvelos (trial_start_times, trial_end_times, gains, g_t_av_f)

    #color background of trials
    color_trial_background(trial_start_times, trial_end_times, type_trials)

    #horizontal x=0 line
    #axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)

    #plot limits
    plt.ylim([-500, 500])
    plt.xlim([0, reg_t[-1]])

    #plot label
    axs.set_xlabel('time (s)',fontsize=20)
    axs.set_ylabel('angular velocity (deg/s)',fontsize=20)
    axs.set_title(str(title),fontsize=20)

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

    plt.show()
