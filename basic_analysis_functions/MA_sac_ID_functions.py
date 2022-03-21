#import packages
import numpy as np
from scipy.interpolate import interp1d
import warnings
import scipy.signal as sig
import scipy
from scipy.signal import butter, filtfilt
import import_functions.MA_get_angvel as ma_get_angvel

###############################################################################

def butter_lowpass(th_order=4, highcut=12., fs=30.):
    nyq             = 0.5 * fs
    high            = highcut / nyq
    return butter(th_order, high, btype='low')

def overfiltVec(inArr, kFf, highcut):
    d, c            = butter_lowpass(th_order=4, highcut=highcut , fs=30.)
    return filtfilt(d, c, inArr)

#get overfiltered ang velo
def overfilter_ang_vel(ang_vel, w):
    #Savitzky–Golay filter
    window = w
    polyorder = 1
    overfiltered_angvel = sig.savgol_filter(ang_vel , window, polyorder, 0)

    return overfiltered_angvel
