import pickle
import warnings
import numpy as np
###############################################################################

def open_pickel(pickle_file):
    while True:
        try:
            all_processed_rawdata = (pickle.load(pickle_file))
        except EOFError:
            break
    pickle_file.close()
    return all_processed_rawdata

###############################################################################

def get_median_angvelos(all_velos):
    angvelo_median =(np.nanmean(all_velos, axis=0))
    return angvelo_median

###############################################################################

def get_trial_rise_frame(angvelo, pattern_angular_velocity):
    tc = 0.63
    if pattern_angular_velocity>0:
        rise_frame = np.where(angvelo>=pattern_angular_velocity*tc)[0][0]
    elif pattern_angular_velocity<0:
        rise_frame = np.where(angvelo<=pattern_angular_velocity*tc)[0][0]
    else:
        rise_frame = np.where(np.abs(angvelo)<=np.abs(100*tc))[0][0]

    return rise_frame

###############################################################################

def get_3_med_ang_vels (trial_ang_velo):
    vi = int(len(trial_ang_velo)/3)
    median_v1 = np.nanmedian(trial_ang_velo[0:vi])
    median_v2 = np.nanmedian(trial_ang_velo[vi:2*vi])
    median_v3 = np.nanmedian(trial_ang_velo[2*vi:3*vi])
    return median_v1, median_v2, median_v3
