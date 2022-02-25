#import packages
import numpy as np
import warnings
from fnmatch import fnmatch
import h5py

##########################################################################################

#import functions
#used to import data to python notebooks

def get_topics_from_hdf5_to_dict(datapaths):
#gets data from hdf5 file to arrays and generates a list per topic
#each list len is the number of files in dataset
#it then saves all the topics in a dictionary

    # all_params_ts = []
    # all_data_params = []

    all_ros_ts = []
    all_ts_bag = []
    all_elapsed_time = []
    all_trial_index = []
    all_trial_elapsed_time = []
    all_panels_running = []
    all_panels_started = []
    all_panels_stopped = []

    all_magnotether_angle = []
    all_magnotether_ros_tstamps = []
    all_magnotether_bag_tstamps = []

    all_ledpanels_ros_tstamps = []
    all_ledpanels_command = []
    all_ledpanels_1 = []
    all_ledpanels_2 = []
    all_ledpanels_3 = []
    all_ledpanels_4 = []
    all_ledpanels_5 = []
    all_ledpanels_6 = []

    for i in range(len(datapaths)):
        f = h5py.File(datapaths[i], "r")

        #parameters topic
    #     params_ts = np.asarray(f['data_params_ros_tstamps'])
    #     data_params = np.asarray(f['data_params'])

        #virtual_desert topic
        ros_ts = np.asarray(f['ros_tstamps'])
        ts_bag = np.asarray(f['tstamps'])
        elapsed_time = np.asarray(f['elapsed_time'])
        trial_index = np.asarray(f['current_trial_index'])
        trial_elapsed_time = np.asarray(f['trial_e_time'])

        panels_running = np.asarray(f['panels_action_running'])
        panels_started = np.asarray(f['panels_action_started'])
        panels_stopped = np.asarray(f['panels_action_stopped'])

        #magnotether_angle topic
        magnotether_angle = np.asarray(f['magnotether_angle'])
        magnotether_ros_tstamps = np.asarray(f['magnotether_ros_tstamps'])
        magnotether_bag_tstamps = np.asarray(f['magnotether_tstamps'])

        #ledpanels topic
        ledpanels_ros_tstamps = np.asarray(f['ledpanels_ros_tstamps'])
        ledpanels_command = np.asarray(f['ledpanels_panels_command'])
        ledpanels_1 = np.asarray(f['ledpanels_panels_arg1'])
        ledpanels_2 = np.asarray(f['ledpanels_panels_arg2'])
        ledpanels_3 = np.asarray(f['ledpanels_panels_arg3'])
        ledpanels_4 = np.asarray(f['ledpanels_panels_arg4'])
        ledpanels_5 = np.asarray(f['ledpanels_panels_arg5'])
        ledpanels_6 = np.asarray(f['ledpanels_panels_arg6'])

    #     all_params_ts.append(params_ts)
    #     all_data_params.append(data_params)

        all_ros_ts.append(ros_ts)
        all_ts_bag.append(ts_bag)
        all_elapsed_time.append(elapsed_time)
        all_trial_index.append(trial_index)
        all_trial_elapsed_time.append(trial_elapsed_time)
        all_panels_running.append(panels_running)
        all_panels_started.append(panels_started)
        all_panels_stopped.append(panels_stopped)
        all_magnotether_angle.append(magnotether_angle)
        all_magnotether_ros_tstamps.append(magnotether_ros_tstamps)
        all_magnotether_bag_tstamps.append(magnotether_bag_tstamps)
        all_ledpanels_1.append(ledpanels_1)
        all_ledpanels_2.append(ledpanels_2)
        all_ledpanels_3.append(ledpanels_3)
        all_ledpanels_4.append(ledpanels_4)
        all_ledpanels_5.append(ledpanels_5)
        all_ledpanels_6.append(ledpanels_6)
        all_ledpanels_command.append(ledpanels_command)
        all_ledpanels_ros_tstamps.append(ledpanels_ros_tstamps)

    all_data = {}
    all_data["all_bag_ts"] = all_ts_bag
    all_data["all_ros_ts"] = all_ros_ts
    all_data["all_elapsed_time"] = all_elapsed_time
    all_data["all_trial_index"] = all_trial_index
    all_data["all_trial_elapsed_time"] = all_trial_elapsed_time
    all_data["all_panels_running"] = all_panels_running
    all_data["all_panels_started"] = all_panels_started
    all_data["all_panels_stopped"] = all_panels_stopped
    all_data["all_magnotether_angle"] = all_magnotether_angle
    all_data["all_magnotether_ros_tstamps"] = all_magnotether_ros_tstamps
    all_data["all_magnotether_bag_tstamps"] = all_magnotether_bag_tstamps
    all_data["all_ledpanels_1"] = all_ledpanels_1
    all_data["all_ledpanels_2"] = all_ledpanels_2
    all_data["all_ledpanels_3"] = all_ledpanels_3
    all_data["all_ledpanels_4"] = all_ledpanels_4
    all_data["all_ledpanels_5"] = all_ledpanels_5
    all_data["all_ledpanels_6"] = all_ledpanels_6
    all_data["all_ledpanels_command"] = all_ledpanels_command
    all_data["all_ledpanels_ros_tstamps"] = all_ledpanels_ros_tstamps

    return  all_data
