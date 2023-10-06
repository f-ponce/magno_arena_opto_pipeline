import pickle
import warnings
import numpy as np
import scipy
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

def flatten_list(_2d_list):
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list

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

###############################################################################

def get_test_trial_median_ma081621(all_velos, all_test_trials, exp_start_frames, exp_end_frames):

    all_fly_trials = []
    for i in range(len(all_velos)):
        fly_trials = []
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial - 2]
            e = exp_end_frames[trial + 1]

            #get the velo
            v = all_velos[i][s:e]

            fly_trials.append(v)

        all_fly_trials.append(fly_trials)

    all_len_fly_trials = []
    for i in range(len(all_fly_trials)):
        for j in range(len(all_fly_trials[i])):
            lt = len(all_fly_trials[i][j])
            all_len_fly_trials.append(lt)

    smallest_trial_fly = min(all_len_fly_trials)
    #smallest_trial_fly = 4571

    all_fly_trials_n = []
    all_fly_trials_n_ci = []
    for i in range(len(all_fly_trials)):
        fly_trials_n = []
        trial_stdev = []
        for j in range(len(all_fly_trials[i])):
            vn = all_fly_trials[i][j][0:smallest_trial_fly]
            fly_trials_n.append(vn)

        all_fly_trials_n.append(np.nanmean((fly_trials_n), axis = 0))

        ###for ci
        mean_y = np.nanmean((all_fly_trials_n), axis = 0)
        wn = 0.05
        b, a = scipy.signal.butter(3, wn)
        mean_y = scipy.signal.filtfilt(b,a,mean_y)
        third = int(0.25*len(mean_y))
        baseline_shift = np.mean( mean_y[0:1800] )
        mean_y -= baseline_shift

        all_fly_trials_n_ci.append(mean_y)

    trial_median = np.nanmean(np.asarray(all_fly_trials_n), axis = 0)

    ci = get_95_confidence_intervals(all_fly_trials_n, iterations=100)

    return trial_median, ci

###############################################################################

def get_test_trial_n_median_ma081621(all_norm_test_trial_traces, all_test_trials, exp_start_frames, exp_end_frames):

    all_fly_trials = []
    for i in range(len(all_velos)):
        fly_trials = []
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial - 2]
            e = exp_end_frames[trial + 1]

            #get the velo
            v = all_velos[i][s:e]

            fly_trials.append(v)

        all_fly_trials.append(fly_trials)

    all_len_fly_trials = []
    for i in range(len(all_fly_trials)):
        for j in range(len(all_fly_trials[i])):
            lt = len(all_fly_trials[i][j])
            all_len_fly_trials.append(lt)

    smallest_trial_fly = min(all_len_fly_trials)
    #smallest_trial_fly = 4571

    all_fly_trials_n = []
    all_fly_trials_n_ci = []
    for i in range(len(all_fly_trials)):
        fly_trials_n = []
        trial_stdev = []
        for j in range(len(all_fly_trials[i])):
            vn = all_fly_trials[i][j][0:smallest_trial_fly]
            fly_trials_n.append(vn)

        all_fly_trials_n.append(np.nanmean((fly_trials_n), axis = 0))

        ###for ci
        mean_y = np.nanmean((all_fly_trials_n), axis = 0)
        wn = 0.05
        b, a = scipy.signal.butter(3, wn)
        mean_y = scipy.signal.filtfilt(b,a,mean_y)
        third = int(0.25*len(mean_y))
        baseline_shift = np.mean( mean_y[0:1800] )
        mean_y -= baseline_shift

        all_fly_trials_n_ci.append(mean_y)

    trial_median = np.nanmean(np.asarray(all_fly_trials_n), axis = 0)

    ci = get_95_confidence_intervals(all_fly_trials_n, iterations=100)

    return trial_median, ci

###############################################################################

def get_test_trial_median_ma081621_2(all_velos, all_test_trials, exp_start_frames, exp_end_frames):

    all_fly_trials = []
    for i in range(len(all_velos)):
        fly_trials = []
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial - 2]
            e = exp_end_frames[trial + 1]

            #get the velo
            v = all_velos[i][s:e]

            fly_trials.append(v)

        all_fly_trials.append(fly_trials)

    all_len_fly_trials = []
    for i in range(len(all_fly_trials)):
        for j in range(len(all_fly_trials[i])):
            lt = len(all_fly_trials[i][j])
            all_len_fly_trials.append(lt)

    smallest_trial_fly = min(all_len_fly_trials)
    #smallest_trial_fly = 4571

    all_fly_trials_n = []
    all_fly_trials_n_ci = []
    for i in range(len(all_fly_trials)):
        fly_trials_n = []
        trial_stdev = []
        for j in range(len(all_fly_trials[i])):
            vn = all_fly_trials[i][j][0:smallest_trial_fly]
            fly_trials_n.append(vn)

        all_fly_trials_n.append(np.nanmedian((fly_trials_n), axis = 0))

        ###for ci
        mean_y = np.nanmedian((all_fly_trials_n), axis = 0)
        wn = 0.05
        b, a = scipy.signal.butter(3, wn)
        mean_y = scipy.signal.filtfilt(b,a,mean_y)
        third = int(0.25*len(mean_y))
        baseline_shift = np.median( mean_y[0:1800] )
        mean_y -= baseline_shift

        all_fly_trials_n_ci.append(mean_y)

    trial_median = np.nanmedian(np.asarray(all_fly_trials_n), axis = 0)

    ci = get_95_confidence_intervals_median(all_fly_trials_n, iterations=100)

    return trial_median, ci
###############################################################################

def get_one_bootstrapped_mean(data):
    indices = np.random.randint(0, len(data), len(data))
    selection = [data[i] for i in indices]
    return np.mean(selection, axis=0)

def get_95_confidence_intervals(data, iterations=100):
    bootstrapped_data = np.mean(data, axis=0)
    for iteration in range(iterations):
        bootstrapped_data = np.vstack( (bootstrapped_data, get_one_bootstrapped_mean(data)) )
    bootstrapped_data.sort(axis=0)

    index_hi = int(iterations*0.975)
    index_lo = int(iterations*0.025)

    return bootstrapped_data[index_lo, :], bootstrapped_data[index_hi, :]

def get_one_bootstrapped_median(data):
    indices = np.random.randint(0, len(data), len(data))
    selection = [data[i] for i in indices]
    return np.median(selection, axis=0)

def get_95_confidence_intervals_median(data, iterations=100):
    bootstrapped_data = np.median(data, axis=0)
    for iteration in range(iterations):
        bootstrapped_data = np.vstack( (bootstrapped_data, get_one_bootstrapped_median(data)) )
    bootstrapped_data.sort(axis=0)

    index_hi = int(iterations*0.975)
    index_lo = int(iterations*0.025)

    return bootstrapped_data[index_lo, :], bootstrapped_data[index_hi, :]
###############################################################################
#get drop time
def get_all_decay_times(all_velos, test_trials, exp_start_frames, exp_end_frames, at_t, fps):

    all_drop_times = []
    for i in range(len(all_velos)):
        fly_drop_times = []
        for j in range(len(test_trials[i])):

            trial = test_trials[i][j]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial+1]
            e = exp_end_frames[trial+1]

            #get the velo
            v = all_velos[i][s:e]

            drop_v = np.nanmedian(v[at_t*fps-fps:(at_t+10)*fps+fps])

            fly_drop_times.append(drop_v)
        all_drop_times.append(np.nanmedian(fly_drop_times))
    return all_drop_times
###############################################################################

#normalization of test trial angular velocity
def normalize_trial_traces(ang_velo, pat_ang_velo, dark_median_ang_velo):
    n_ang_velo = (ang_velo - dark_median_ang_velo)/(pat_ang_velo - dark_median_ang_velo)

    return n_ang_velo

def get_all_norm_test_trial_traces(all_velos, all_test_trials, exp_start_frames, exp_end_frames, all_pattern_velos):

    all_norm_test_trial_traces = []
    for i in range(len(all_velos)):
        norm_test_trial_traces = []
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            pat_ang_velo = all_pattern_velos[i][trial]

            #get start/end frames of trials to plot
            s = exp_start_frames[trial]
            e = exp_end_frames[trial]
            #get start/end frames of trials to plot
            s_d = exp_start_frames[trial-1]
            e_d = exp_end_frames[trial-1]

            #get the velo
            v = all_velos[i][s:e]
            v_d = all_velos[i][s_d:e_d]

            median_v_d1, median_v_d2, median_v_d3  = get_3_med_ang_vels(v_d)

            n_v = normalize_trial_traces(v, pat_ang_velo, median_v_d3)

            norm_test_trial_traces.append(n_v)
        all_norm_test_trial_traces.append(norm_test_trial_traces)

    return all_norm_test_trial_traces

###############################################################################


###############################################################################
###############################################################################
###############################################################################

def get_test_trial_2medians_ma022622(all_velos, all_test_trials, all_start_frames, all_end_frames):

    all_fly_trials = []
    for i in range(len(all_velos)):
        fly_trials = []
        for j in range(len(all_test_trials[i])):

            trial = all_test_trials[i][j]

            #get start/end frames of trials to plot
            s = all_start_frames[i][trial - 1]
            e = all_end_frames[i][trial + 2]

            #get the velo
            v = all_velos[i][s:e]

            fly_trials.append(v)

        all_fly_trials.append(fly_trials)

    all_len_fly_trials = []
    for i in range(len(all_fly_trials)):
        for j in range(len(all_fly_trials[i])):
            lt = len(all_fly_trials[i][j])
        all_len_fly_trials.append(lt)

    smallest_trial_fly = min(all_len_fly_trials)

    all_fly_trials_n = []
    for i in range(len(all_fly_trials)):
        fly_trials_n = []
        trial_stdev = []
        for j in range(len(all_fly_trials[i])):
            vn = all_fly_trials[i][j][0:smallest_trial_fly]
            fly_trials_n.append(vn)
        all_fly_trials_n.append(np.nanmean(np.asarray(fly_trials_n), axis = 0))

    trial_median = np.nanmean(np.asarray(all_fly_trials_n), axis = 0)
    ci = get_95_confidence_intervals(all_fly_trials_n, iterations=100)

    return trial_median, ci
