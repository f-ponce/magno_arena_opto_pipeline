#import packages
import numpy as np
from scipy.interpolate import interp1d
import warnings
import scipy.signal as sig
import scipy
from scipy.signal import butter, filtfilt

###############################################################################

def butter_lowpass(th_order=4, highcut=12., fs=30.):
    nyq             = 0.5 * fs
    high            = highcut / nyq
    return butter(th_order, high, btype='low')

def overfiltVec(inArr, kFf, highcut):
    d, c            = butter_lowpass(th_order=4, highcut=0.5 , fs=30.)
    return filtfilt(d, c, inArr)
