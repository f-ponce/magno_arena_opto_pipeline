import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib
import itertools

import basic_analysis_functions.MA_sac_ID_functions as ma_saccf
###############################################################################
#global params
sns.set_theme(font='Arial', style='ticks')
plt.rcParams["font.family"] = 'Arial'
fs = 8
cm = 1/2.54
lw = 1.0
tml = 3 #tick mark len
tlw = 0.5 #trace line width
tlw_decay = 0.7
tpd = 2 #tick padding
lbpd = 2 #axis labelpad
p_lw = 1.0 #pattern linewidth
patcolor = '#31a354'
save_dir = '/Users/fponce/Desktop/MI_figs'

###############################################################################

def plot_medang_dur_sd(test_trials_mmeans, test_trials_ci_95s, exp_s):

    save_im_name = 'ang_vel_test_trials_'+str(exp_s)+'.svg'

    fig = plt.figure(figsize=((7)*cm, (3)*cm))
    gs = gridspec.GridSpec(nrows=1, ncols=1, left=0.25, right=0.94, top=0.9, bottom=0.3)
    axs = fig.add_subplot(gs[0, :])

    test_time = np.linspace(0, (120+30+exp_s), (120+30+exp_s)*30)

    for tt in range(len(test_trials_mmeans)):

        axs.plot(test_time, test_trials_mmeans[tt], color='k', zorder=2, linewidth=tlw)
        axs.plot((0, test_time[-1]), (0, 0), 'gray', linewidth=0.5, zorder=1)
        #ci
        axs.fill_between(test_time, test_trials_ci_95s[tt][0], test_trials_ci_95s[tt][1],  color= 'gray')

    set_pattern_vels(axs, exp_s)
    set_axes_params_angvel(axs, exp_s)
    set_axes_colors(axs, exp_s)

    #fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600)
    plt.show()

###############################################################################

def plot_medang_dur_sd_control(test_trials_mmeans, test_trials_ci_95s, exp_s):

    save_im_name = 'ang_vel_control_trials_'+str(exp_s)+'.svg'

    fig = plt.figure(figsize=((7)*cm, (3)*cm))
    gs = gridspec.GridSpec(nrows=1, ncols=1, left=0.2, right=0.94, top=0.9, bottom=0.3)
    axs = fig.add_subplot(gs[0, :])

    test_time = np.linspace(0, (120+30+exp_s), (120+30+exp_s)*30)

    for tt in range(len(test_trials_mmeans)):

        axs.plot(test_time, test_trials_mmeans[tt], color='k', zorder=2, linewidth=tlw)
        axs.plot((0, test_time[-1]), (0, 0), 'gray', linewidth=0.5, zorder=1)
        #ci
        axs.fill_between(test_time, test_trials_ci_95s[tt][0], test_trials_ci_95s[tt][1],  color= 'gray')

    set_pattern_vels2(axs, exp_s)
    set_axes_params_angvel(axs, exp_s)
    set_axes_colors(axs, exp_s)

    #fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600, transparent=True)
    plt.show()

###############################################################################

# def plot_exp_decay_traces(all_processed_data, colors):
#
#     save_im_name = 'decay_traces'+'.svg'
#
#     fig = plt.figure(figsize=((7)*cm, (3)*cm))
#     mywidths = [1,0.07]
#     gs = gridspec.GridSpec(nrows=1, ncols=2, left=0.2, right=0.85, top=0.9, bottom=0.25, wspace=0.25, width_ratios=mywidths)
#     ax0 = fig.add_subplot(gs[0, 0])
#     ax1 = fig.add_subplot(gs[0, 1])
#
#     for i in range(len(all_processed_data)):
#         exp_mm = all_processed_data[i]["test_trials_pool_mmean"][-1840:-1]
#         t = all_processed_data[i]["reg_t"][0:len(exp_mm)]
#
#         #plot traces
#         ax0.plot(t,ma_saccf.overfilter_ang_vel(exp_mm, 17*1), c = colors[i])
#
#     #plot pattern angular pattern_angular_velocity
#     ax0.plot((t[0], t[50]), (182.5, 182.5), color=patcolor, linewidth= 2)
#
#     #plot colorbar
#     gradient = np.linspace(0, 1, 256)
#     gradient = np.vstack((gradient)*-1)
#     ax1.imshow(gradient[1:-2], aspect='auto', cmap=plt.get_cmap('plasma'))
#
#     #plot formatting
#     set_axes_params_decay_angvel(ax0)
#     set_axes_params_colorbar(ax1)
#     make_box(ax0)
#
#     fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600, transparent=True)
#     plt.show()

def plot_exp_decay_traces(all_processed_data, colors):

    save_im_name = 'decay_traces'+'.svg'

    fig = plt.figure(figsize=((5.5)*cm, (4.5)*cm))
    gs = gridspec.GridSpec(nrows=1, ncols=1, left=0.2, right=0.94, top=0.9, bottom=0.2)
    ax0 = fig.add_subplot(gs[0, 0])

    for i in range(len(all_processed_data)):
        exp_mm = all_processed_data[i]["test_trials_pool_mmean"][-1840:-1]
        t = all_processed_data[i]["reg_t"][0:len(exp_mm)]

        #plot traces
        ax0.plot(t,ma_saccf.overfilter_ang_vel(exp_mm, 17*1), c = colors[i], linewidth=tlw_decay)

    #plot pattern angular pattern_angular_velocity
    ax0.plot((t[0], t[50]), (182.5, 182.5), color=patcolor, linewidth= 2)

    # #plot colorbar
    # gradient = np.linspace(0, 1, 256)
    # gradient = np.vstack((gradient)*-1)
    # ax1.imshow(gradient[1:-2], aspect='auto', cmap=plt.get_cmap('plasma'))

    #plot formatting
    set_axes_params_decay_angvel(ax0)
    # set_axes_params_colorbar(ax1)
    make_box(ax0)

    fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600, transparent=True)
    plt.show()

###############################################################################

def plot_exp_decay_times(all_processed_data, colors, all_mtests, ts):

    save_im_name = 'decay_times'+'.svg'

    fig = plt.figure(figsize=((7.5)*cm, (4.5)*cm))
    gs = gridspec.GridSpec(nrows=1, ncols=1, left=0.16, right=0.98, top=0.9, bottom=0.2)
    ax0 = fig.add_subplot(gs[0, 0])

    for i in range(len(all_processed_data)):

        test1 = all_processed_data[i]["test_trials_decayt_50"][0]
        test2 = all_processed_data[i]["test_trials_decayt_50"][1]
        tests = all_mtests[i]

        mtest1 = np.median(test1)
        mtest2 = np.median(test2)
        mtest = np.median([mtest1, mtest2])
        tt = [test1, test2]
        ftt =  list(itertools.chain(*tt))


        randnums = np.random.uniform(ts[i]-1, ts[i]+1, size=(len(tests)))
        randnums1 = np.random.uniform(i, i-0.5, size=(len(test1)))
        randnums2 = np.random.uniform(i, i-0.5, size=(len(test2)))

        plt.scatter((randnums), tests, color = colors[i], alpha=0.6, s =15)
        plt.scatter((ts[i]), mtest, color = 'k', s=15, marker='_')

        #################
        #plot line at 0
        #axs.axhline(y=0.0, color='gray', linestyle='-', linewidth= 0.9)
        ax0.plot((3,180),(182.5,182.5), color=patcolor, linewidth=p_lw, zorder=0)
        set_axes_params_decay_times(ax0)

    fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600, transparent=True)
    plt.show()

###############################################################################

def plot_medang_longdur(t, med_angvel, ci_95s):

    save_im_name = 'ang_vel_long_duration2'+'.svg'

    fig = plt.figure(figsize=((13)*cm, (3)*cm))
    gs = gridspec.GridSpec(nrows=1, ncols=1, left=0.14, right=0.96, top=0.8, bottom=0.3)
    axs = fig.add_subplot(gs[0, :])

    #axs.plot(t, med_angvel, color='k', zorder=3, linewidth=tlw)
    #axs.plot((0, t[-1]), (0, 0), 'gray', linewidth=0.5, zorder=1)
    #ci
    axs.fill_between(t, ci_95s[0], ci_95s[1],  color= 'gray')

    set_pattern_vels3(axs)
    set_axes_params_angvel2(axs, t[-1])
    set_axes_colors2(axs)

    fig.savefig(save_dir+'/'+save_im_name, format='svg', dpi=600, transparent=True)
    plt.show()
###############################################################################
###############################################################################
###############################################################################
#formatting functions

def set_axes_params_angvel(ax, exp_s, fs=fs, lw=lw):
    sns.despine(ax=ax, offset=5, trim=False)
    ax.tick_params(axis='y', direction='in', length=tml, width=lw, pad=tpd)
    ax.tick_params(axis='x', direction='in', length=tml, width=lw, pad=tpd)
    txlim = 60+30+60+exp_s
    ax.spines['left'].set_bounds(-182.5, 182.5)
    ax.spines['bottom'].set_bounds(0, txlim)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.set_xlim([0,txlim])
    ax.set_xticks([0,60,90,(90+exp_s+60)])
    ax.set_xticklabels(['0','60','90',str(90+exp_s+60)])
    ax.set_ylim([-300,300])
    ax.set_yticks([-182.5,0,182.5])
    ax.set_yticklabels(['-182.5','0','182.5'])
    ax.set_xlabel('Time (s)', fontsize=fs, labelpad=lbpd)
    ax.set_ylabel('Ang Vel (degs/s)', fontsize=fs, labelpad=lbpd, loc='center')

def set_axes_params_angvel2(ax, txlim, fs=fs, lw=lw):
    sns.despine(ax=ax, offset=5, trim=False)
    ax.tick_params(axis='y', direction='in', length=tml, width=lw, pad=tpd)
    ax.tick_params(axis='x', direction='in', length=tml, width=lw, pad=tpd)
    ax.spines['left'].set_bounds(-182.5, 182.5)
    ax.spines['bottom'].set_bounds(0, txlim)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.set_xlim([0,1260])
    ax.set_xticks([0,180,360,1260])
    ax.set_xticklabels(['0','180','360','1260'])
    ax.set_ylim([-300,300])
    ax.set_yticks([-182.5,0,182.5])
    ax.set_yticklabels(['-182.5','0','182.5'])
    ax.set_xlabel('Time (s)', fontsize=fs, labelpad=lbpd)
    ax.set_ylabel('Ang Vel\n(degs/s)', fontsize=fs, labelpad=lbpd, loc='center')

def set_pattern_vels(ax, exp_s):
    ax.plot((60,90),(0,0), color=patcolor, linewidth=p_lw, zorder=1)
    ax.plot((90,90),(0,182.5), color=patcolor, linewidth=p_lw, zorder=1)
    ax.plot((90,90+exp_s),(182.5,182.5), color=patcolor, linewidth=p_lw, zorder=1)
    ax.plot((90,90),(-182.5,0), color=patcolor, linewidth=p_lw, zorder=1)
    ax.plot((90,90+exp_s),(-182.5,-182.5), color=patcolor, linewidth=p_lw, zorder=1)

def set_pattern_vels2(ax, exp_s):
    ax.plot((60,90),(0,0), color=patcolor, linewidth=lw, zorder=1)
    ax.plot((90,90+exp_s),(0,0), color=patcolor, linewidth=p_lw, zorder=1)

def set_pattern_vels3(ax):
    ax.plot((180,360),(-182.5,-182.5), color=patcolor, linewidth=lw, zorder=1)

def set_axes_colors(ax, exp_s):
    ax.axvspan((60), (90), facecolor='seagreen', alpha=0.1, edgecolor='None')
    ax.axvspan((90), (90+exp_s), facecolor='seagreen', alpha=0.2, edgecolor='None')

def set_axes_colors2(ax):
    ax.axvspan((180), (360), facecolor='seagreen', alpha=0.2, edgecolor='None')

def set_axes_params_decay_angvel(ax, fs=fs, lw=lw):
    sns.despine(ax=ax, offset=5, trim=False)
    ax.tick_params(axis='y', direction='in', length=tml, width=lw, pad=tpd)
    ax.tick_params(axis='x', direction='in', length=tml, width=lw, pad=tpd)
    ax.spines['left'].set_bounds(0, 200)
    ax.spines['bottom'].set_bounds(2, 62)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.set_xlim([0,63])
    ax.set_xticks([2,62])
    ax.set_xticklabels(['0','60'])
    ax.axvspan((0), (2), facecolor='seagreen', alpha=0.1, edgecolor='None')
    ax.set_ylim([0, 200])
    ax.set_yticks([0,100,200])
    ax.set_yticklabels(['0','100','200'])
    ax.set_xlabel('Time (s)', fontsize=fs, labelpad=lbpd)
    ax.set_ylabel('Ang Vel (degs/s)', fontsize=fs, labelpad=lbpd, loc='center')

def make_box(ax):
    ax.plot((52,62),(1,1), color='k', linewidth=0.5, zorder=3)
    ax.plot((52,62),(200,200), color='k', linewidth=0.5, zorder=3)
    ax.plot((52,52),(1,200), color='k', linewidth=0.5, zorder=3)
    ax.plot((62,62),(1,200), color='k', linewidth=0.5, zorder=3)

def set_axes_params_colorbar(ax, fs=fs, lw=lw):
    sns.despine(ax=ax, offset=5, trim=False, right=False, left=False)
    ax.tick_params(axis='y', direction='in', length=tml, width=lw, pad=tpd)
    ax.tick_params(axis='x', direction='in', length=0, width=0)
    ax.yaxis.tick_right()
    ax.yaxis.set_ticks_position('right')
    ax.spines['right'].set_bounds(0, 182.5)
    #ax.spines['bottom'].set_bounds(4, 64)
    ax.spines['right'].set_linewidth(lw)
    ax.spines['left'].set_linewidth(0)
    ax.spines['right'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(0)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.set_xticklabels([])
    ax.set_ylim([0, 182.5])
    ax.set_yticks([0,182.5])
    ax.set_yticklabels(['3','180'])
    ax.yaxis.set_label_position("right")
    ax.set_ylabel('Bias Duration (s)', fontsize=fs, labelpad=-8)

def set_axes_params_decay_times(ax, fs=fs, lw=lw):
    sns.despine(ax=ax, offset=5, trim=False)
    ax.tick_params(axis='y', direction='in', length=tml, width=lw, pad=tpd)
    ax.tick_params(axis='x', direction='in', length=tml, width=lw, pad=tpd)
    ax.spines['left'].set_bounds(0, 250)
    ax.spines['bottom'].set_bounds(3, 180)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)
    ax.yaxis.set_tick_params(labelsize=fs)
    ax.xaxis.set_tick_params(labelsize=fs)
    ax.set_xlim([0,190])
    ax.set_xticks([3,10,30,60,87,100,180])
    ax.set_xticklabels(['3','10','30','60','90','100','180'])
    ax.set_ylim([-4, 250])
    ax.set_yticks([0,100,250])
    ax.set_yticklabels(['0','100','250'])
    ax.set_xlabel('Bias Duration (s)', fontsize=fs, labelpad=lbpd, loc='center')
    ax.set_ylabel('Ang Vel (degs/s)', fontsize=fs, labelpad=lbpd, loc='center')

def remove_inner_ticklabels(fig):
    for ax in fig.axes:
        try:
            ax.label_outer()
        except:
            pass
