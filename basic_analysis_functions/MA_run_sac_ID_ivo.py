import numpy as np
import basic_analysis_functions.MA_sac_ID_ivo as ma_sacc
import basic_analysis_functions.MA_sac_ID_functions as ma_saccf

#run saccade ider on angular velocity traces

def get_all_overfiltered_velos(all_velos, dt, highcut):
    #overfiltering angular velocity
    all_overfilt_velos = []
    for i in range(len(all_velos)):
        overfilt_velo = ma_saccf.overfiltVec(all_velos[i], dt, 0.5)
        all_overfilt_velos.append(overfilt_velo)

    return all_overfilt_velos

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
