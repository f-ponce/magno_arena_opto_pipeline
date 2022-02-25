#run saccade ider on angular velocity traces


#overfiltering angular velocity
all_overfilt_velos = []
for i in range(len(all_velos)):
    fps = reg_t[1] - reg_t[1]
    overfilt_velo = ma_filter.overfiltVec(all_velos[i], fps)
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

    all_nosaccs_velos.append(velo_nosaccs1)
    all_saccs_lr_idx.append(saccs_lr)
    all_saccs_l_idx.append(saccs_l)
    all_saccs_r_idx.append(saccs_r)
    all_saccs_thresh.append(headingVeloThresholdDeg)
