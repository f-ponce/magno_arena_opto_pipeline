import numpy as np
import copy
from scipy.interpolate import interp1d
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_sac_ID_functions as ma_saccf

#run saccade ider on angular velocity traces

def get_all_overfiltered_velos(all_velos, window):
    #overfiltering angular velocity
    all_overfilt_velos = []
    for i in range(len(all_velos)):
        #overfilt_velo = ma_saccf.overfiltVec(all_velos[i], dt, highcut)
        overfilt_velo = ma_saccf.overfilter_ang_vel(all_velos[i], window)
        all_overfilt_velos.append(overfilt_velo)

    return all_overfilt_velos

# def get_all_overfiltered_velos_test(all_velos, window, highcut):
#     f = 30.
#     dt = 1/30.
#     #overfiltering angular velocity
#     all_overfilt_velos = []
#     for i in range(len(all_velos)):
#         overfilt_velo = ma_saccf.overfilter_ang_vel(all_velos[i], window)
#         overfilt_velo = ma_saccf.overfiltVec(overfilt_velo, f, highcut)
#         all_overfilt_velos.append(overfilt_velo)
#
#     return all_overfilt_velos

def get_all_nosacc_velos(all_velos, all_overfilt_velos, threshold):

    #run saccade ider on every fly
    all_nosaccs_velos = []
    all_saccs_lr_idx = []
    all_saccs_l_idx = []
    all_saccs_r_idx = []
    all_saccs_thresh = []
    for i in range(len(all_velos)):
        #sacc detection threshold
        VeloThresholdDeg = np.nanstd(all_velos[i] - all_overfilt_velos[i])*threshold

        #detrend angular velocity trace
        detrended_angvelo = all_velos[i] - all_overfilt_velos[i]
        #run saccade ider function on detrended angular velocity
        velo_detrended_nosaccs, saccs_lr, saccs_l, saccs_r = ma_sacc.findSacs(detrended_angvelo, all_velos[i], VeloThresholdDeg)
        #find and remove saccade from actual angular velocity signal
        velo_nosaccs = np.copy(all_velos[i])
        velo_nosaccs[np.isnan(velo_detrended_nosaccs)] = np.nan

        all_nosaccs_velos.append(velo_nosaccs)
        all_saccs_lr_idx.append(saccs_lr)
        all_saccs_l_idx.append(saccs_l)
        all_saccs_r_idx.append(saccs_r)
        all_saccs_thresh.append(VeloThresholdDeg)

    return all_nosaccs_velos, all_saccs_lr_idx, all_saccs_l_idx, all_saccs_r_idx, all_saccs_thresh

def get_all_nosacc_velos_loop(all_velos, reg_t):

    dt = reg_t[1] - reg_t[0]
    fps = 1/dt

    #overfilter angular velocity data
    highcut = 0.2
    window = 17*3
    all_overfilt_velos_1 = get_all_overfiltered_velos(all_velos, window)

    threshold = 3.5
    all_nosaccs_velos_1, all_saccs_idx_1,  all_saccs_l_idx_1, all_saccs_r_idx_1, all_saccs_thresh_1 = get_all_nosacc_velos(all_velos,
                                                                                                        all_overfilt_velos_1, threshold)

    ############################################################################
    #interpolating velo to run again through saccade ider
    all_nosaccs_velos_1_interp = []
    for i in range(len(all_nosaccs_velos_1)):
        v1 = all_nosaccs_velos_1[i]
        v1_nan_idx = ~np.isnan(v1)

        v = np.array(v1[v1_nan_idx])
        t = np.array(reg_t[v1_nan_idx])

        f_a = interp1d(t, v, bounds_error=False)
        reg_a = f_a(reg_t)

        all_nosaccs_velos_1_interp.append(reg_a)
    ############################################################################
    highcut = 0.2
    window = 17*2+1
    all_overfilt_velos = get_all_overfiltered_velos(all_nosaccs_velos_1_interp, window)

    #remove saccade
    threshold = 3.
    all_nosaccs_velos_2, all_saccs_idx_2,  all_saccs_l_idx_2, all_saccs_r_idx_2, all_saccs_thresh_2 = get_all_nosacc_velos(all_nosaccs_velos_1,
                                                                                                        all_overfilt_velos, threshold)

    ############################################################################
    #add indices found in 2 saccade ider runs
    all_saccs_idx = []
    all_saccs_l_idx = []
    all_saccs_r_idx = []
    for i in range(len(all_saccs_idx_2)):
        saccs_id = all_nosaccs_velos_2[i] + all_nosaccs_velos_1[i]
        #saccs_l_idx = all_saccs_l_idx_2[i] + all_saccs_l_idx_1[i]
        saccs_l_idx = np.array([ a1 if not np.isnan(a1) else b1 for a1,b1 in zip(all_saccs_l_idx_2[i], all_saccs_l_idx_1[i])])
        #saccs_r_idx = all_saccs_l_idx_2[i] + all_saccs_l_idx_1[i]
        saccs_r_idx = np.array([ a1 if not np.isnan(a1) else b1 for a1,b1 in zip(all_saccs_r_idx_2[i], all_saccs_r_idx_1[i])])
        all_saccs_idx.append(saccs_id)
        all_saccs_l_idx.append(saccs_l_idx)
        all_saccs_r_idx.append(saccs_r_idx)

    #remove saccades from original velocity

    all_nosaccs_velos = copy.deepcopy(all_velos)
    for i in range(len(all_nosaccs_velos)):
        idxs = np.isnan(all_saccs_idx[i])
        all_nosaccs_velos[i][idxs] = np.nan

    return all_nosaccs_velos, all_saccs_idx, all_saccs_l_idx, all_saccs_r_idx, all_overfilt_velos

def get_all_interp_after_saccs(all_nosaccs_velos, reg_t):
    #interpolating velo to run again through saccade ider
    all_nosaccs_velos_interp = []
    for i in range(len(all_nosaccs_velos)):
        v1 = all_nosaccs_velos[i]
        v1_nan_idx = ~np.isnan(v1)

        v = np.array(v1[v1_nan_idx])
        t = np.array(reg_t[v1_nan_idx])

        f_a = interp1d(t, v, bounds_error=False)
        reg_a = f_a(reg_t)

        all_nosaccs_velos_interp.append(reg_a)

    return all_nosaccs_velos_interp
