#import packages
import numpy as np
from scipy.interpolate import interp1d
import warnings
import scipy.signal as sig
import scipy
from scipy.signal import butter, filtfilt

###############################################################################

def get_all_magnotether_filt_angles(all_magnotether_angles_uw):
#all_magnotether_angle: list of magnotether angles

    all_magnotether_angles_filt = []
    for i in range(len(all_magnotether_angles_uw)):
        magnotether_angles_filt = filter_angles_savgol(all_magnotether_angles_uw[i])
        all_magnotether_angles_filt.append(magnotether_angles_filt)
    return all_magnotether_angles_filt

def filter_angles_savgol(angles_deg_uw):
    #unwrap angles
    #angles_unwrapped = np.rad2deg(np.unwrap(np.deg2rad(angles_deg)))
    angles_unwrapped = angles_deg_uw
    #Savitzky–Golay filter
    window = 17
    polyorder = 3
    filtered_angles = sig.savgol_filter(angles_unwrapped , window, polyorder, 0)

    #wrap angles
    filtered_angles_rad = np.deg2rad(filtered_angles)
    filtered_angles_wr_rad = np.arctan2(np.sin(filtered_angles_rad), np.cos(filtered_angles_rad))
    filtered_angles_wr_deg = np.rad2deg(filtered_angles_wr_rad)

    return filtered_angles_wr_deg
