import sys
import numpy as np
import scipy as sp
import scipy.signal as sig
import matplotlib.pyplot as plt

def get_angvel(t,p,window=17,polyorder=3):
    p = p + 180.0
    p_2x = (2*p) % 360.0
    p_2x_unwrap = np.unwrap(p_2x,discont=180.0)
    p_2x_filt = sig.savgol_filter(p_2x_unwrap, window, polyorder, 0)
    dp_filt = sig.savgol_filter(p_2x_unwrap, window, polyorder, 1)/2.0
    info = {'p_2x': p_2x_unwrap, 'p_2x_filt': p_2x_filt}
    return dp_filt, info

#get filtered ang velo from all files in the list
def get_all_angvel(reg_t, all_magnotether_angles_exp_uw_interp):
    dt = reg_t[1] - reg_t[0]
    all_angvels = []
    for i in range(len(all_magnotether_angles_exp_uw_interp)):
        #wrap angles
        angles_rad = np.deg2rad(all_magnotether_angles_exp_uw_interp[i])
        angles_wr_rad = np.arctan2(np.sin(angles_rad), np.cos(angles_rad))
        angles_wr_deg = np.rad2deg(angles_wr_rad)
        #filter angular velocity of each file/fly
        w = 17
        po = 3
        ang_vel, info = get_angvel(reg_t, angles_wr_deg,window=w,polyorder=po)
        all_angvels.append(ang_vel/dt)

    return all_angvels


# datafile = sys.argv[1]
# data = np.loadtxt(datafile)
#
# t = data[:,0]
# t = t - t[0]
# p = data[:,1]
#
# dt = np.diff(t).mean()
# print(dt)
#
# dp_filt, info = get_angvel(t,p)
# dp_filt =dp_filt/dt
#
#
# fig, ax = plt.subplots(2,1,sharex=True)
# ax[0].plot(t,info['p_2x'],'b')
# ax[0].plot(t,info['p_2x_filt'],'r')
# ax[0].set_ylabel('2x position (p)')
# ax[0].grid(True)
#
# ax[1].plot(t,dp_filt)
# ax[1].set_ylabel('angular velocity (p)')
# ax[1].grid(True)
# plt.show()
