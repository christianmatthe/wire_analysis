########## imports and setup
import numpy as np
import matplotlib.pyplot as plt
import os
import time
#from datetime import datetime, tzinfo
import datetime as dt
from scipy.interpolate import interp1d
import dill
#from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from scipy.special import binom
from numpy.linalg import inv
import json

# import functions from other file
# import Voltage_base_analysis as vba
# import calib_analysis as ca
# import pressure_analysis as pa

from .Voltage_base_analysis import (load_data, save_data, select_date_indices,
                                    prep_data_slowdash, )

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)

# DEPRECATED
# plot_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep 
#             + "output/flow_on_off/")
# data_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep + "data/"
# os.makedirs(plot_dir, exist_ok=True)
#######################

def date_picker(date_list, start_date, end_date):
    return np.array([int(i) for i,date in enumerate(date_list)
              if  start_date < date < end_date])

def make_index_chunks(data_dict, start_date, time_interval, n_cycles):
    chunk_dict={}
    for i in range(n_cycles):
        chunk_dict[f"{i}"] = {}
        for j, state in enumerate(["off", "on"]):
            date_indices = date_picker(data_dict["dates"]
                            , start_date + (2 * i + j) * time_interval
                            , start_date + (2* i + j + 1) * time_interval)
            # chunk_dict[f"{i}{state}"] = date_indices
            chunk_dict[f"{i}"][state] = date_indices
    return chunk_dict

# def plot_colored(data_dict, chunk_dict, plotname):

#     fig = plt.figure(0, figsize=(8,6.5))
#     ax1=plt.gca()
#     for i,key in enumerate(chunk_dict.keys()):
#         # labels = ["off", "on"]
#         # if i < 1:
#         #     label = labels[i%2]
#         # else:
#         #     label = "_nolegend_"
#         for j, state in enumerate(["off", "on"]):
#             if state == "off":
#                 color = "C0"
#             if state == "on":
#                 color = "C1"
#             if i < 1:
#                 label = state
#             else:
#                 label = "_nolegend_"
#             ax1.plot(data_dict["dates"][chunk_dict[key][state]],
#                     data_dict["voltage"][chunk_dict[key][state]],
#                     ".",# markersize=1,
#                     color = color,
#                     label = label
#                     )
#     # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
#     #         label = f"moving average {mavg_len}s")

#     plt.xticks(rotation = 45)

#     ax1.set_xlabel(r"Time")
#     ax1.set_ylabel(r"Voltage [mV]")

#     plt.grid(True)
#     plt.legend(shadow=True)
#     plt.tight_layout()
#     format_im = 'png' #'pdf' or png
#     dpi = 300
#     plt.savefig(plot_dir + plotname + "_colored"
#                 + '.{}'.format(format_im),
#                 format=format_im, dpi=dpi)
#     ax1.cla()

# def plot_colored_R_hack(data_dict, chunk_dict, plotname):

#     fig = plt.figure(0, figsize=(8,6.5))
#     ax1=plt.gca()
#     for i,key in enumerate(chunk_dict.keys()):
#         # labels = ["off", "on"]
#         # if i < 1:
#         #     label = labels[i%2]
#         # else:
#         #     label = "_nolegend_"
#         for j, state in enumerate(["off", "on"]):
#             if state == "off":
#                 color = "C0"
#             if state == "on":
#                 color = "C1"
#             if i < 1:
#                 label = state
#             else:
#                 label = "_nolegend_"
#             ax1.plot(data_dict["dates"][chunk_dict[key][state]],
#                     data_dict["voltage"][chunk_dict[key][state]],
#                     ".",# markersize=1,
#                     color = color,
#                     label = label
#                     )
#     # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
#     #         label = f"moving average {mavg_len}s")

#     plt.xticks(rotation = 45)

#     ax1.set_xlabel(r"Time")
#     ax1.set_ylabel(r"Resistance [$\Omega$]")

#     plt.grid(True)
#     plt.legend(shadow=True)
#     plt.tight_layout()
#     format_im = 'png' #'pdf' or png
#     dpi = 300
#     plt.savefig(plot_dir + plotname + "_colored_R_hack"
#                 + '.{}'.format(format_im),
#                 format=format_im, dpi=dpi)
#     ax1.cla()

# def calc_avg(data_dict, chunk_dict):
#     avg_dict = {"v" : [], "v_err" : []}
#     states = ["off", "on"]
#     # ncut = 60
#     # new assymetric cuut to prioritize mostly equilibrated  section
#     # 3 min fronnt cut "optimized" for 4 min measurement runs
#     n_cut_front = 60 * 3
#     n_cut_back = 30
#     keys = chunk_dict.keys()
#     for i,key in enumerate(keys):
#         for j, state in enumerate(states):
#             vs = data_dict["voltage"][chunk_dict[key][state]][
#                             n_cut_front:-n_cut_back]
#             v_avg = np.average(vs)
#             v_err = np.std(vs)

#             avg_dict["v"].append(v_avg)
#             avg_dict["v_err"].append(v_err)
#     # convert to nummpy  arrays
#     for key in avg_dict:
#         avg_dict[key] = np.array(avg_dict[key])
#     return avg_dict

# def calc_avg_for_drift_10(data_dict, chunk_dict):
#     """
#     The idea is to cut 4 x 2 min averages from every 10 min measurement.
#     ncut has to be such, that section can be cut with even spacing: 
#     n_cut = len_section/n_sections
#     """
#     avg_dict = {"v" : [], "v_err" : []}
#     states = ["off", "on"]

#     n_sec = 4
#     n_cut = 30
#     keys = chunk_dict.keys()
#     for i,key in enumerate(keys):
#         for j, state in enumerate(states):
#             vs = data_dict["voltage"][chunk_dict[key][state]][n_cut:-n_cut]
#             #v_avg_1 = np.average(vs[len(vs)//2])
#             v_avg_list = [np.average(
#                         vs[i*len(vs)//n_sec:(i+1)*len(vs)//n_sec])
#                         for i in range(n_sec)
#                         ]
#             v_err = np.std(vs)

#             for v_avg in v_avg_list:
#                 avg_dict["v"].append(v_avg)
#                 avg_dict["v_err"].append(v_err)
#     # convert to nummpy  arrays
#     for key in avg_dict:
#         avg_dict[key] = np.array(avg_dict[key])
#     return avg_dict

# def calc_diffs(avg_dict):
#     diff_dict = {"dv" : [], "dv_err" : []}
#     for j in range(1, len(avg_dict["v"])-1):
#         # subtract average of adjacent
#         dv = avg_dict["v"][j]-(avg_dict["v"][j-1] +  avg_dict["v"][j-2]) / 2
#         dv_err = np.sqrt(avg_dict["v_err"][j] **2 
#                         +(avg_dict["v_err"][j-1] / 2) **2
#                         +(avg_dict["v_err"][j+1] / 2) **2
#                         )
#         diff_dict["dv"].append(dv)
#         diff_dict["dv_err"].append(dv_err)
#             # convert to nummpy  arrays
#     for key in diff_dict:
#         diff_dict[key] = np.array(diff_dict[key])
#     return diff_dict

# def plot_diffs(diff_dict, plotname, exclude_list = [],
#                 labels = ["on","off"], color_increment = 0
#                 ):

#     fig = plt.figure(0, figsize=(8,6.5))
#     ax1=plt.gca()
#     # introduce plotcounter to stop legend when appropriate
#     plot_counter = 0
#     for i in range(len(diff_dict["dv"])):
#         if i in exclude_list:
#             pass
#         else:
#             if plot_counter < 2:
#                 label = labels[plot_counter]
#             else:
#                 label = "_nolegend_"
#             ax1.errorbar(i,
#                     diff_dict["dv"][i],
#                     yerr = diff_dict["dv_err"][i],
#                     fmt = ".",# markersize=1,
#                     color = f"C{(i + color_increment)%2}" ,
#                     label = label
#                     )
#             plot_counter += 1
        
#     # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
#     #         label = f"moving average {mavg_len}s")

#     #plt.xticks(rotation = 45)

#     ax1.set_xlabel(r"Cycle")
#     ax1.set_ylabel(r"Voltage difference [mV]")

#     plt.grid(True)
#     plt.legend(shadow=True)
#     plt.tight_layout()
#     format_im = 'png' #'pdf' or png
#     dpi = 300
#     plt.savefig(plot_dir + plotname + "_diffs"
#                 + '.{}'.format(format_im),
#                 format=format_im, dpi=dpi)
#     ax1.cla()

# def plot_diffs_R_hack(diff_dict, plotname, exclude_list = [],
#                         labels = ["on","off"], color_increment = 0
#                         ):

#     fig = plt.figure(0, figsize=(8,6.5))
#     ax1=plt.gca()
#     plot_counter = 0
#     for i in range(len(diff_dict["dv"])):
#         if i in exclude_list:
#             pass
#         else:
#             if plot_counter < 2:
#                 label = labels[plot_counter]
#             else:
#                 label = "_nolegend_"
#             ax1.errorbar(i,
#                     diff_dict["dv"][i],
#                     yerr = diff_dict["dv_err"][i],
#                     fmt = ".",# markersize=1,
#                     color = f"C{(i + color_increment)%2}" ,
#                     label = label
#                     )
#             plot_counter += 1
        
#     # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
#     #         label = f"moving average {mavg_len}s")

#     #plt.xticks(rotation = 45)

#     ax1.set_xlabel(r"Cycle")
#     ax1.set_ylabel(r"Resistance difference [$\Omega$]")

#     plt.grid(True)
#     plt.legend(shadow=True)
#     plt.tight_layout()
#     format_im = 'png' #'pdf' or png
#     dpi = 300
#     plt.savefig(plot_dir + plotname + "_diffs_R_hack"
#                 + '.{}'.format(format_im),
#                 format=format_im, dpi=dpi)
#     ax1.cla()
            
# def avg_diff(diff_dict, exclude_list = []):
#     index_list_on = [i for i in range(0,len(diff_dict["dv"]),2)
#                   if not (i in exclude_list)]
#     index_list_off = [i for i in range(1,len(diff_dict["dv"]),2)
#                 if not (i in exclude_list)]
#     print(index_list_on, index_list_off)
#     dv_avg_on = np.average(diff_dict["dv"][index_list_on])
#     dv_avg_off = np.average(diff_dict["dv"][index_list_off])
#     dv_err_avg_on = (np.sqrt(np.sum(diff_dict["dv_err"][index_list_on] **2)) / 
#                      np.sqrt(len(index_list_on)))
#     dv_err_avg_off =(np.sqrt(np.sum(diff_dict["dv_err"][index_list_off] **2))/ 
#                      np.sqrt(len(index_list_off)))
#     return dv_avg_on, dv_avg_off, dv_err_avg_on, dv_err_avg_off

# # implement polynomial drift filter as in https://arxiv.org/pdf/1009.1894.pdf
# # Note there will be someissues with too few data points
# def A_drift_filter(N,p):
#     A = [  
#         [
#         (1/(2*p))*binom(p,j-i)*(-1)**j if (0 <=(j-i) <= p) else 0
#         for j in range(0,N)
#         ]
#         for i in range(0,N-p)
    
#         ]
#     return np.asmatrix(A)

# def A_drift_filter_10_min_sectioned(N,p):
#     """
#     (j//4)%2 flips the sign every 4 iterations since every  measurement is cut
#     into 4 to determine drift.

#     Idea seems broken
#     """
#     A = [  
#         [
#         (1/(2*p))*binom(p,j-i)*(-1)**(j//4)%2 if (0 <=(j-i) <= p) else 0
#         for j in range(0,N)
#         ]
#         for i in range(0,N-p)
    
#         ]
#     return np.asmatrix(A)

# def U_drift_filter(avg_dict):
#     U = np.asmatrix(avg_dict["v"]).T
#     return U

# def drift_filter(avg_dict,p=3):
#     U = U_drift_filter(avg_dict)
#     N = len(U)
#     A = A_drift_filter(N,p)
#     #print("U: ", U)
#     #print("A: ", A)
#     Y = np.asmatrix(A)*np.asmatrix(U)
#     X = np.asmatrix([1 for i in range(N-p)]).T
#     mu = X.T* inv(A*A.T) * Y /(X.T * inv(A*A.T) * X)
#     # Calculate Errors
#     # W  =X.T* inv(A*A.T) /(X.T * inv(A*A.T) * X)
#     # print(np.asmatrix(np.cov(Y,Y)), type(np.asmatrix(np.cov(Y,Y))))
#     # print(W, type(W))
#     print(mu[0,0])
#     print(Y.shape, X.shape,inv(A*A.T).shape, ((Y - mu[0,0]*X).T).shape )
#     s2 = ((Y - mu[0,0]*X).T * inv(A*A.T) * (Y - mu[0,0]*X))/(N-p-1)
#     # # Doesn't work as I expected
#     # var = W * np.asmatrix(np.cov(Y)) * W.T
#     # sig = np.sqrt(np.abs(var))
#     var = s2 / (X.T * inv(A*A.T) * X)
#     sig = np.sqrt(var)



#     return [mu, Y, sig]

# def drift_filter_10_min_sectioned(avg_dict,p=3):
#     U = U_drift_filter(avg_dict)
#     N = len(U)
#     A = A_drift_filter_10_min_sectioned(N,p)
#     # print("U: ", U)
#     # print("A: ", A)
#     Y = np.asmatrix(A)*np.asmatrix(U)
#     X = np.asmatrix([1 for i in range(N-p)]).T
#     mu = X.T*(A*A.T)**-1 * Y /(X.T * (A*A.T)**-1 * X)
#     return mu


# def make_HABS_power(file_path):
#     # Function to load a file conntaining HBS current andVoltage outputs and
#     # writing power to a dict
#     with open(file_path) as f:
#         raw_dict  = json.load(f)

#     dates_c = np.array([pa.str_to_dt(ls[0])
#                 for ls in raw_dict["habs_current_output"]])
#     print("dates_c", dates_c)
#     dates_v = np.array([pa.str_to_dt(ls[0])
#                 for ls in raw_dict["habs_voltage_output"]])
#     # move to German timezone
#     utc_offset = 2
#     dates_c = np.array([(
#             date.astimezone(dt.timezone(
#                 dt.timedelta(hours=utc_offset))) 
#                 + dt.timedelta(hours=utc_offset))
#             for date in dates_c])
#     dates_v = np.array([(
#             date.astimezone(dt.timezone(
#                 dt.timedelta(hours=utc_offset))) 
#                 + dt.timedelta(hours=utc_offset))
#             for date in dates_v])
#     current = np.array([float(ls[1]) 
#                     for ls in raw_dict["habs_current_output"]])
#     voltage = np.array([float(ls[1]) 
#                     for ls in raw_dict["habs_voltage_output"]])

#     # only pick values where both Current and voltage are recorded
#     P_dict = {"dates": [], "current": [], "voltage": [], "power":[]}
#     for i_c,date in enumerate(dates_c):
#         print(i_c)
#         if date in dates_v:
#             i_v = np.where(dates_v == date)[0]
#             P_dict["dates"].append(date)
#             P_dict["current"].append(current[i_c])
#             P_dict["voltage"].append(current[i_v])
#             power = current[i_c] * current[i_v]
#             P_dict["power"].append(power)
    
#     for key in P_dict.keys():
#         P_dict[key] = np.array(P_dict[key])
    
#     return P_dict

# def make_HABS_power_1_per_min(file_path):
#     # Function to load a file conntaining HBS current andVoltage outputs and
#     # writing power to a dict
#     with open(file_path) as f:
#         raw_dict  = json.load(f)
#         data = raw_dict

#     #  Old date import
#     # # only take every 60th element of raw dict
#     # dates_c = np.array([pa.str_to_dt(ls[0])
#     #             for ls in raw_dict["habs_current_output"]])
#     # dates_v = np.array([pa.str_to_dt(ls[0])
#     #             for ls in raw_dict["habs_voltage_output"]])
#     # # move to German timezone
#     # utc_offset = 2
#     # dates_c = np.array([(
#     #         date.astimezone(dt.timezone(
#     #             dt.timedelta(hours=utc_offset))) 
#     #             + dt.timedelta(hours=utc_offset))
#     #         for date in dates_c])
#     # dates_v = np.array([(
#     #         date.astimezone(dt.timezone(
#     #             dt.timedelta(hours=utc_offset))) 
#     #             + dt.timedelta(hours=utc_offset))
#     #         for date in dates_v])

#     # new date import from slowdash
#     t_series_c = (np.array(data["habs_current_output"]["t"]) 
#                + float(data["habs_current_output"]["start"]))
#     dates_c = mpl.dates.num2date(t_series_c/86400)
#     t_series_v = (np.array(data["habs_voltage_output"]["t"]) 
#                + float(data["habs_voltage_output"]["start"]))
#     dates_v = mpl.dates.num2date(t_series_v/86400)

#     # OLD
#     # current = np.array([float(ls[1]) 
#     #                 for ls in raw_dict["habs_current_output"]
#     #                 ])
#     # voltage = np.array([float(ls[1]) 
#     #                 for ls in raw_dict["habs_voltage_output"]
#     #                 ])

#     current = np.array(data["habs_current_output"]["x"])
#     # voltage = np.array(data["habs_voltage_output"]["x"])
#     # print(type(current[1]))
#     # print("voltage:", voltage)

#     def v_interpolate(data):
#         t_arr = (np.array(data["habs_voltage_output"]["t"]) 
#                + float(data["habs_voltage_output"]["start"]))
#         f_int = interp1d(t_arr,
#                         data["habs_voltage_output"]["x"],
#                         kind = "linear",fill_value=(np.NaN,np.NaN),
#                         bounds_error= False
#                 )
#         return f_int
#     v_int = v_interpolate(data)
#     voltage = v_int(t_series_c)

#     # only pick values where both Current and voltage are recorded
#     P_dict = {"dates": [], "current": [], "voltage": [], "power":[]}
#     # # take only every 60th to speed up
#     # for i_c,date in enumerate(dates_c[::60]):
#     #     #print(i_c)
#     #     if date in dates_v:
#     #         i_v = np.where(dates_v == date)[0]
#     #         P_dict["dates"].append(date)
#     #         P_dict["current"].append(current[i_c*60])
#     #         P_dict["voltage"].append(voltage[i_v])
#     #         power = current[i_c*60] * voltage[i_v]
#     #         P_dict["power"].append(power)

#     P_dict["dates"] = dates_c
#     P_dict["current"] = current
#     P_dict["voltage"] = voltage
#     power = np.array(current) * np.array(voltage)
#     P_dict["power"] = power
    
#     for key in P_dict.keys():
#         P_dict[key] = np.array(P_dict[key])
    
#     return P_dict

def save_dict(dict, filename):
        with open(filename + ".pkl", "wb") as f:
            dill.dump(dict, f)

def load_dict(filename):
        with open(filename + ".pkl", "rb") as f:
            dict = dill.load(f)
        return dict

def power_to_temperature(P):
    # From old calibration
    #  https://discourse.project8.org/t/mainz-atomic-cracker-habs/290
    T = (270.9304338715183 * P**(1/2.627253832618772) 
         + 311.9698446637183 + (-0.31725280518758997*P) )
    return T


# TODO GOAL Implement piecewise background fitting of equilibrated sections in
# "Off" state
# Then Fit offset to "ON" state 

# Simplest approach ABA Fit 0th order to A's then calc average offset to B
# Steps:
# 1. Select Data
    # a) Section by time interval
    # b) Chop off edges and equilibrization time

class Extractor():
    """
    "The Extractor class contains functions and parameters necessary to
    extract flow ON-OFF Resistance/Voltage offset.
    "

Parameters
    ----------
    data_dict :  `dict`
        dictionary containing the data
    
    """

    def __init__(self,
                 data_dict,
                 start_date, 
                 end_date,
                 run_name,
                 plot_dir,
                 time_interval = dt.timedelta(minutes = 4),
                 utc_offset = 1,
                 initial_state = 0,
                 front_crop = dt.timedelta(minutes = 2),
                 rear_crop = dt.timedelta(minutes  =0, seconds  =30),
                 denoise = False
                 ):

        #no defaults
        self.data_dict = data_dict
        self.start_date = start_date
        self.end_date = end_date
        self.run_name = run_name

        # with  defaults
        self.time_interval = time_interval
        self.utc_offset = utc_offset
        self.initial_state = initial_state
        self.front_crop = front_crop
        self.rear_crop = rear_crop

        # Make output dirs
        self.plot_dir = plot_dir
        # Deliberately commented out to breakk bad functionality
        # self.data_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep 
        #             + "data/")
        os.makedirs(self.plot_dir, exist_ok=True)



        # Constants 
        self.state_names = ["off", "on"]

        # (pseudo)slots
        self.basic_ABA_fit_dict = {}
        self.quad_ABA_fit_dict = {}
        self.lin_exp_ABA_fit_dict = {}
        self.exp_back_ABA_fit_dict = {}
        self.fit_results = {}

        #denoise?
        if denoise == True:
            self.denoise()

    def denoise(self):
        print("checking if denoised file exists")
        # Check if data was already denoised
        dataname = self.run_name + "_data_denoised"
        if os.path.isfile(data_dir + dataname + ".pkl"):
            ##### SIGALRM only works in Unix
            # import signal
            # TIMEOUT = 30 # number of seconds your want for timeout

            # def interrupted(signum, frame):
            #     # "called when read times out"
            #     print('interrupted!')
            # signal.signal(signal.SIGALRM, interrupted)

            # def timeout_input(prompt):
            #     try:
            #             print('You have 30 seconds to respond')
            #             var = input(prompt)
            #             return var
            #     except:
            #             # timeout
            #             return "no"

            # # set alarm
            # signal.alarm(TIMEOUT)
            # s = timeout_input("A denoised dataset already exists for this run_name."
            #     + "Do you want to load it? (type 'y')")
            # # disable the alarm after success
            # signal.alarm(0)

            s = input("A denoised dataset already exists for this run_name."
                + "Do you want to load it? (type 'y')")
            if s == "y":
                print(f"loading {dataname}")
                self.data_dict = load_data(dataname)
                return
        
        print(f"denoising into file {dataname}")
        # Problem: This seems to be massively slowy: 5-10mins for 1day 
        # data
        # Note this  denoises the ENTIRE data set not the slice we are 
        # actually
        # looking at.
        # It goes a lot faster on the 2nd sectionn, so it likely only 
        # does it 
        # once and "data_dict" is permanently altered
        v_series = self.data_dict["voltage"]
        # Remove downward dip systematic noise pattern I call "fuzz"
        # If point is lower than both its neighbors by an average
        # of more than 2µV: (arbitrary cut)
        # Cut point before and 2 after 
        # Will result in ~50% loss of data for expected noise
        dip_mins = [i for i in range(1,len(v_series)-1)
                    if v_series[i] < ((v_series[i-1] + v_series[i+1])/2 
                                    - 2e-6)
                        ]
        #print(dip_mins)
        fuzz_mask = [[i-1, i, i +1, i+2] 
                    if i+2 < len(v_series) 
                    else [i-1, i, i +1]
                    for i in dip_mins
                    ]
        # flatten fuzz mask and remove duplicates
        fuzz_mask = list(set(np.concatenate(fuzz_mask)))
        # invert fuzz mask, now only keepps points that we want to keep
        fuzz_mask = [i for i in range(1, len(v_series)-1) 
                    if i not in fuzz_mask]
        #print(fuzz_mask)

        # adjust data_dict
        for key in self.data_dict.keys():
            self.data_dict[key] = self.data_dict[key][fuzz_mask]

        # bodge: save data dict
        # denoising should probably be a separate step outside th extractor?
        save_data(dataname, self.data_dict)



    def slice_dict(self):
        """"
        Slice data by time intervals and sort into on and off states 

        NOTE: This seems slower than  entirely necessary. Date selection?
        """
        initial_state = self.initial_state
        #sliced_dict = {name:[] for name in self.state_names}
        sliced_dict = {}
        # Make time interval edges
        n_slices = int(np.rint((self.end_date - self.start_date)
                            / self.time_interval))
        for n in range(n_slices):
            interval_start = self.start_date + n * self.time_interval
            if (interval_start + self.time_interval) < self.end_date:
                interval_end = interval_start + self.time_interval
            else:
                interval_end = self.end_date
            start,end = select_date_indices(self.data_dict,
                                                interval_start,
                                                interval_end)
            v_series = self.data_dict["voltage"][start:end] 
            dates = self.data_dict["dates"][start:end]
            # move to German timezone
            dates = np.array([
                    date.astimezone(dt.timezone(
                        dt.timedelta(hours=self.utc_offset)))
                    for date in dates])
            # build sclided dict in 2 separate pieces, one "OFF" , one "ON"
            # Is this even a good idea?
            # sliced_dict[self.state_names[(n + initial_state)%2]][n//2] = {
            #     "voltage":v_series, "dates":dates

            # }
            #just put slices in a list?
            # sliced_dict.append({"voltage":v_series, "dates":dates,
            #  "state": self.state_names[(n + initial_state)%2] })

            sliced_dict[n] = {"voltage":v_series, "dates":dates,
             "state": self.state_names[(n + initial_state)%2] }
        
        #pass to object
        self.sliced_dict = sliced_dict
        return sliced_dict

    def crop_slice(self,
                   slice, 
                   front_crop = dt.timedelta(minutes = 2),
                   rear_crop = dt.timedelta(minutes  =0, seconds  =30),
                    ):
        """"
        crop equilibrization times off front and rear of data sets 
        """
        interval_start = slice["dates"][0] + front_crop
        interval_end = slice["dates"][-1] - rear_crop
        start,end = select_date_indices(slice,
                                                interval_start,
                                                interval_end)
        v_series = slice["voltage"][start:end] 
        dates = slice["dates"][start:end]
        # make copy of slice dictionary
        crop = slice.copy()
        # replace changed entries
        crop["voltage"] = v_series
        crop["dates"] = dates
        return crop

    def crop_dict(self,
                ):

        cropped_dict = {}
        for n in self.sliced_dict.keys():
            slice = self.sliced_dict[n]
            crop = self.crop_slice(slice, front_crop=self.front_crop,
                            rear_crop = self.rear_crop,
                            )
            cropped_dict[n] = crop

        #pass to object
        self.cropped_dict = cropped_dict
        return cropped_dict
    
    def lin_exp_crop_dict(self,
                front_crop_A = dt.timedelta(minutes = 2),
                front_crop_B = dt.timedelta(minutes  =0, seconds  =30),
                rear_crop = dt.timedelta(minutes  =0, seconds  =30),
                ):

        cropped_dict = {}
        for n in self.sliced_dict.keys():
            if n %2 == 0: # Hack assumes As are even slices
                slice = self.sliced_dict[n]
                crop = self.crop_slice(slice, front_crop=front_crop_A,
                                rear_crop = rear_crop,
                                )
                cropped_dict[n] = crop
            if n %2 == 1:
                slice = self.sliced_dict[n]
                crop = self.crop_slice(slice, front_crop=front_crop_B,
                                rear_crop = rear_crop,
                                )
                cropped_dict[n] = crop

        #pass to object
        self.cropped_dict = cropped_dict
        return cropped_dict

    def basic_ABA_fit(self,
        slice_keys = [0,1,2],
        ):
        """"
        fit 1st order baseline through "A"'s and fit offset to "B"  
        """
        ti_sec = self.time_interval.total_seconds()
        def fit_func(t, c0, c1,
                 B_offset,   
                ):
            # separate function by on or off state
            # identify  B state
            f = np.piecewise(t,
                            [np.logical_and(ti_sec < t, t < 2 * ti_sec),
                            np.logical_not(np.logical_and(
                                ti_sec < t, t < 2 * ti_sec))],
                            [lambda t: (c1*t + c0) + B_offset,
                            lambda t: (c1*t + c0)]
                            )
            # if ti_sec < t < 2 * ti_sec:
            #     f = (c1*t + c0) + B_offset
            # else:
            #     f = (c1*t + c0)
            return f

        # convert dates to seconds since start of procedure
        start_date = self.sliced_dict[slice_keys[0]]["dates"][0]
        pre_con = (
            [np.array(self.cropped_dict[key]["dates"]) - start_date
            for key in slice_keys]
            )
        con = np.concatenate(pre_con)
        time_deltas = np.array([delta.total_seconds() for delta in con])
        # make continuous v_series
        v_series = np.array(np.concatenate(
            [self.cropped_dict[key]["voltage"] for key in slice_keys]))


        popt, pcov = curve_fit(fit_func, time_deltas, v_series,
                                        p0 = [80
                                                ,0
                                                ,0]
                                        #,sigma = P_err
                                        #, absolute_sigma = True
                                        )
        perr = np.sqrt(np.diag(pcov))
        # print("popt: ", popt)
        # print("perr: ", perr)

        #define lin space for the whole ABA section
        t_space = np.linspace(0* ti_sec,
                              3* ti_sec,
                              num=int(3* ti_sec))
        fit_series = np.array(
                    [fit_func(t, popt[0],popt[1],popt[2]) for t in t_space])
        # Save results in dictionary. First slice key becomes the key
        self.basic_ABA_fit_dict[slice_keys[0]] = {"popt": popt, "perr": perr,
                                   "start_date": start_date,
                                    "t_space": t_space,
                                    "fit_series": fit_series,
                                    "method": "basic"
        }
        return self.basic_ABA_fit_dict

    def quad_ABA_fit(self,
        slice_keys = [0,1,2],
        ):
        """"
        fit 1st order baseline through "A"'s and fit offset to "B"  
        """
        ti_sec = self.time_interval.total_seconds()
        def fit_func(t, c0, c1, c2,
                 B_offset,   
                ):
            # separate function by on or off state
            # identify  B state
            f = np.piecewise(t,
                            [np.logical_and(ti_sec < t, t < 2 * ti_sec),
                            np.logical_not(np.logical_and(
                                ti_sec < t, t < 2 * ti_sec))],
                            [lambda t: (c2*t**2 + c1*t + c0) + B_offset,
                            lambda t: (c2*t**2 +c1*t + c0)]
                            )
            return f

        # convert dates to seconds since start of procedure
        start_date = self.sliced_dict[slice_keys[0]]["dates"][0]
        pre_con = (
            [np.array(self.cropped_dict[key]["dates"]) - start_date
            for key in slice_keys]
            )
        con = np.concatenate(pre_con)
        time_deltas = np.array([delta.total_seconds() for delta in con])
        # make continuous v_series
        v_series = np.array(np.concatenate(
            [self.cropped_dict[key]["voltage"] for key in slice_keys]))


        popt, pcov = curve_fit(fit_func, time_deltas, v_series,
                                        p0 = [80
                                                ,0
                                                ,0
                                                ,0
                                                ]
                                        #,sigma = P_err
                                        #, absolute_sigma = True
                                        )
        perr = np.sqrt(np.diag(pcov))
        # print("popt: ", popt)
        # print("perr: ", perr)

        #define lin space for the whole ABA section
        t_space = np.linspace(0* ti_sec,
                              3* ti_sec,
                              num=int(3* ti_sec))
        fit_series = np.array(
                    [fit_func(t, popt[0],popt[1],popt[2],popt[3]) 
                     for t in t_space])
        # Save results in dictionary. First slice key becomes the key 
        self.quad_ABA_fit_dict[slice_keys[0]] = {"popt": popt, "perr": perr,
                                   "start_date": start_date,
                                    "t_space": t_space,
                                    "fit_series": fit_series,
                                    "method": "quad",
        }
        return self.quad_ABA_fit_dict

    def exp_back_ABA_fit(self,
        slice_keys = [0,1,2],
        front_crop_A = dt.timedelta(minutes = 2),
        front_crop_B = dt.timedelta(minutes  =0, seconds  =30),
        rear_crop = dt.timedelta(minutes  =0, seconds  =30),
        ):
        """"
        fit 1st order baseline + exponential through long term thermalization  
        """
        # first written 13 Jul 2022
        ti_sec = self.time_interval.total_seconds()
        t_front_crop_A = front_crop_A.total_seconds()
        t_front_crop_B = front_crop_B.total_seconds()
        def fit_func(t, c0, c1, A0, tau_A,
                 B_offset,   
                ):
            # fit exponential to background  temperature trend

            f = np.piecewise(t,
                [np.logical_and(ti_sec < t, t < 2 * ti_sec),
                t <= ti_sec , 
                t >= 2 * ti_sec],
                [lambda t: ((c1*t + c0) 
                            + A0 * (-np.exp(
                                -(t)/tau_A)))
                            + B_offset ,

                lambda t: ((c1*t + c0) 
                            + A0 * (-np.exp(
                                -(t )/tau_A))),
                lambda t: ((c1*t + c0)
                            + A0 * (-np.exp(
                                -(t )/tau_A)))
                ]
                )
            return f
        
        # define function for  just the  virtual background
        def v_eq_func(t, c0, c1):
            return (c1*t + c0)

        # convert dates to seconds since start of procedure
        start_date = self.sliced_dict[slice_keys[0]]["dates"][0]
        pre_con = (
            [np.array(self.cropped_dict[key]["dates"]) - start_date
            for key in slice_keys]
            )
        con = np.concatenate(pre_con)
        time_deltas = np.array([delta.total_seconds() for delta in con])
        # make continuous v_series
        v_series = np.array(np.concatenate(
            [self.cropped_dict[key]["voltage"] for key in slice_keys]))

        # t, c0, c1, A0,  tau_A,  B_offset,   
        popt, pcov = curve_fit(fit_func, time_deltas, v_series,
                                        p0 = [80, 0
                                                , 1
                                                , 800
                                                , 1e-4]
                                        #,sigma = P_err
                                        #, absolute_sigma = True
                                        )
        perr = np.sqrt(np.diag(pcov))
        # print("popt: ", popt)
        # print("perr: ", perr)

        #define lin space for the whole ABA section
        t_space = np.linspace(0* ti_sec,
                              3* ti_sec,
                              num=int(3* ti_sec))
        fit_series = np.array(
                    [fit_func(t, popt[0],popt[1],popt[2],popt[3],popt[4]
                    ) for t in t_space])
        v_eq_series = np.array(
                    [v_eq_func(t, popt[0],popt[1]
                    ) for t in t_space])
        # Save results in dictionary. First slice key becomes the key
        self.exp_back_ABA_fit_dict[slice_keys[0]] = {"popt": popt, "perr": perr,
                                   "start_date": start_date,
                                    "t_space": t_space,
                                    "fit_series": fit_series,
                                    "virtual_eq_series": v_eq_series,
                                    "method": "exp_back",
        }
        return self.exp_back_ABA_fit_dict

    def lin_exp_ABA_fit(self,
        slice_keys = [0,1,2],
        front_crop_A = dt.timedelta(minutes = 2),
        front_crop_B = dt.timedelta(minutes  =0, seconds  =30),
        rear_crop = dt.timedelta(minutes  =0, seconds  =30),
        ):
        """"
        fit 1st order baseline through "A"'s and fit offset to "B"  
        """
        # first written 28 Jun 2022
        ti_sec = self.time_interval.total_seconds()
        t_front_crop_A = front_crop_A.total_seconds()
        t_front_crop_B = front_crop_B.total_seconds()
        def fit_func(t, c0, c1, A0, A1, tau_A, tau_B,
                 B_offset,   
                ):
            # separate function by on or off state
            # Need essentially 2 separate fits:
            # 1. Fit exponentials to the A state to determine a 
            #    "virtual baseline" y0
            # 2. Fit exponential to B state choosing t0 such that f(t=t0) = y0
            # In combination this zeros the measurement on the right time and
            # base voltage, eliminating these ambiguous fit  parameters in the
            # B fit

            # Note: This procedure seems overcomplicated and is likely not
            #       robust

            # A_0 and A_1 are the amplitude of the exponential for the
            # A fits

            # Below is the hideously complicated version of getting B_offset as
            # a direct fit parameter
            
            # t0 (t_B here for symmetry) requires this bloated definiton to 
            # eliminate redundant parameters and impose the desired 
            # relationship
            dt0_A = t_front_crop_A
            A_diff = -A0 * (1 - np.exp(-(ti_sec - dt0_A)/tau_A))
            dt0_B = - tau_B * np.log(1-(A_diff/B_offset))

            f = np.piecewise(t,
                [np.logical_and(ti_sec < t, t < 2 * ti_sec),
                t <= ti_sec , 
                t >= 2 * ti_sec],
                [lambda t: ((c1*t + c0) 
                            + B_offset * (1 - np.exp(
                                -(t - (ti_sec - dt0_B))/tau_B))),
                                # sign before dt0_B has to be negative since
                                # it offsets in the opposite direction as dt0_A

                lambda t: ((c1*t + c0) 
                            + A0 * (-np.exp(
                                -(t - (0 * ti_sec + dt0_A))/tau_A))),
                lambda t: ((c1*t + c0)
                            + A1 * (-np.exp(
                                -(t - (2 * ti_sec + dt0_A))/tau_A)))
                ]
                )
            return f

        # define function for  just the  virtual background
        def v_eq_func(t, c0, c1):
            return (c1*t + c0)

        # convert dates to seconds since start of procedure
        start_date = self.sliced_dict[slice_keys[0]]["dates"][0]
        pre_con = (
            [np.array(self.cropped_dict[key]["dates"]) - start_date
            for key in slice_keys]
            )
        con = np.concatenate(pre_con)
        time_deltas = np.array([delta.total_seconds() for delta in con])
        # make continuous v_series
        v_series = np.array(np.concatenate(
            [self.cropped_dict[key]["voltage"] for key in slice_keys]))

        # t, c0, c1, A0, A1, tau_A, tau_B, B_offset,   
        popt, pcov = curve_fit(fit_func, time_deltas, v_series,
                                        p0 = [80, 0
                                                , 1e-5, 1e-5
                                                , 80,200
                                                , 1e-4]
                                        , bounds = (
                                                [-np.inf, -2* 1e-8
                                                , 1e-8, 1e-8
                                                , 1, 1
                                                , -np.inf],
                                                [ np.inf, +2* 1e-8
                                                , 1e-3, 1e-3
                                                , np.inf, np.inf
                                                , np.inf]
                                        )
                                        #,sigma = P_err
                                        #, absolute_sigma = True
                                        )
        perr = np.sqrt(np.diag(pcov))
        # print("popt: ", popt)
        # print("perr: ", perr)

        #define lin space for the whole ABA section
        t_space = np.linspace(0* ti_sec,
                              3* ti_sec,
                              num=int(3* ti_sec))
        fit_series = np.array(
                    [fit_func(t, popt[0],popt[1],popt[2],popt[3],popt[4]
                                ,popt[5]
                                , popt[6]
                    ) for t in t_space])
        v_eq_series = np.array(
                    [v_eq_func(t, popt[0],popt[1]
                    ) for t in t_space])
        # Save results in dictionary. First slice key becomes the key
        self.lin_exp_ABA_fit_dict[slice_keys[0]] = {"popt": popt, "perr": perr,
                                   "start_date": start_date,
                                    "t_space": t_space,
                                    "fit_series": fit_series,
                                    "virtual_eq_series": v_eq_series,
                                    "method": "lin_exp",
        }
        return selquad_ABA_fit_dictf.lin_exp_ABA_fit_dict


    def basic_ABA_fit_all_B(self,
        ):
        # perform fit for all available ABA connections with B centered
        for i in range(1,len(self.sliced_dict.keys()),2):   
            self.basic_ABA_fit(slice_keys = [j for j in range(i-1,i+2)])
        self.write_means_to_fit_results(self.basic_ABA_fit_dict)
        return
    
    def quad_ABA_fit_all_B(self,
        ):
        # perform fit for all available ABA connections with B centered
        for i in range(1,len(self.sliced_dict.keys()),2):   
            self.quad_ABA_fit(slice_keys = [j for j in range(i-1,i+2)])
        self.write_means_to_fit_results(self.quad_ABA_fit_dict)
        return

    def exp_back_ABA_fit_all_B(self,
            front_crop_A = dt.timedelta(minutes = 2),
            front_crop_B = dt.timedelta(minutes  =0, seconds  =30),
            rear_crop = dt.timedelta(minutes  =0, seconds  =30),
        ):
        # perform fit for all available ABA connections with B centered
        for i in range(1,len(self.sliced_dict.keys()),2):   
            self.exp_back_ABA_fit(slice_keys = [j for j in range(i-1,i+2)],
            front_crop_A = front_crop_A,
            front_crop_B = front_crop_B,
            rear_crop = rear_crop,
            )
        self.write_means_to_fit_results(self.exp_back_ABA_fit_dict)
        return

    def lin_exp_ABA_fit_all_B(self,
            front_crop_A = dt.timedelta(minutes = 2),
            front_crop_B = dt.timedelta(minutes  =0, seconds  =30),
            rear_crop = dt.timedelta(minutes  =0, seconds  =30),
        ):
        # perform fit for all available ABA connections with B centered
        for i in range(1,len(self.sliced_dict.keys()),2):   
            self.lin_exp_ABA_fit(slice_keys = [j for j in range(i-1,i+2)],
            front_crop_A = front_crop_A,
            front_crop_B = front_crop_B,
            rear_crop = rear_crop,
            )
        self.write_means_to_fit_results(self.lin_exp_ABA_fit_dict)
        return


    def plot_all_basic_ABA_fit(self,
                        plot_path,
                        method = "basic",
                        plot_v_eq = False,
                        **kwargs
                            ):

        fig = plt.figure(0, figsize=(8,6.5))
        ax1=plt.gca()
        # plot data
        for key,slice in self.sliced_dict.items():
            if key == 0:
                label = "data"
            else:
                label = "_nolegend_"
            ax1.plot(slice["dates"],slice["voltage"]*1000,".",
                markersize=8,
                color = "C0",
                label = label)
        
        data_ylims = ax1.get_ylim()
            
        # plot cropped data
        for key,slice in self.cropped_dict.items():
            if key == 0:
                w_mean = self.fit_results[method]["w_mean"]
                w_std = self.fit_results[method]["w_std"]
                label = ("fit_data" + f", offset w_mean: {w_mean:.3e}" 
                        + r"$\pm$" + f"{w_std:.1e}")
            else:
                label = "_nolegend_"
            ax1.plot(slice["dates"],slice["voltage"]*1000,".",
                markersize=8,
                color = "C1",
                label = label)

        # plot fits
        #color iterator
        i_c = 1
        if method == "exp_back":
            # TODO fix custom cropped data highlighting
            for key,fit in self.exp_back_ABA_fit_dict.items():
                i_c +=1
                dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                v_series = fit["fit_series"]
                ax1.plot(dates,v_series*1000,
                        "-",
                        #markersize=4,
                        linewidth=2,
                        alpha=1,
                        color = f"C{i_c}",
                        label = (f"c0: {fit['popt'][0]:.1e},"
                            + f"c1: {fit['popt'][1]:.1e},"
                            + f"A0: {fit['popt'][2]:.1e},"
                            + f"tau_A: {fit['popt'][3]:.1e},"
                            + f"B_offset: {fit['popt'][-1]:.3e}")
                    )
                # Plot virtual background
                if plot_v_eq == True:
                    dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                    v_series = fit["virtual_eq_series"]
                    ax1.plot(dates,v_series*1000,
                            "-",
                            #markersize=4,
                            linewidth=2,
                            alpha=1,
                            color = f"C{i_c}",
                            label = "_nolegend_"
                        )
                    ax1.plot(dates,(v_series + fit['popt'][-1])*1000 ,
                            "-",
                            #markersize=4,
                            linewidth=2,
                            alpha=1,
                            color = f"C{i_c}",
                            label = "_nolegend_"
                        )
            if plot_v_eq == True:
                data_ylims = ax1.get_ylim()
                    
            plt.legend(shadow=True,loc='lower left', bbox_to_anchor=(0, 1),
                  fontsize=8)
            ax1.set_ylim(data_ylims)


        if method == "lin_exp":
            # TODO fix custom cropped data highlighting
            for key,fit in self.lin_exp_ABA_fit_dict.items():
                i_c +=1
                dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                v_series = fit["fit_series"]
                ax1.plot(dates,v_series*1000,
                        "-",
                        #markersize=4,
                        linewidth=2,
                        alpha=1,
                        color = f"C{i_c}",
                        label = (f"c0: {fit['popt'][0]:.1e},"
                            + f"c1: {fit['popt'][1]:.1e},"
                            + f"A0: {fit['popt'][2]:.1e},"
                            + f"A1: {fit['popt'][3]:.1e},"
                            + f"tau_A: {fit['popt'][4]:.1e},"
                            + f"tau_B: {fit['popt'][5]:.1e},"
                            + f"B_offset: {fit['popt'][-1]:.3e}")
                    )
                # Plot virtual background
                if plot_v_eq == True:
                    dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                    v_series = fit["virtual_eq_series"]
                    ax1.plot(dates,v_series*1000,
                            "-",
                            #markersize=4,
                            linewidth=2,
                            alpha=1,
                            color = f"C{i_c}",
                            label = "_nolegend_"
                        )
                    ax1.plot(dates,(v_series + fit['popt'][-1])*1000 ,
                            "-",
                            #markersize=4,
                            linewidth=2,
                            alpha=1,
                            color = f"C{i_c}",
                            label = "_nolegend_"
                        )
            if plot_v_eq == True:
                data_ylims = ax1.get_ylim()
                    
            plt.legend(shadow=True,loc='lower left', bbox_to_anchor=(0, 1),
                  fontsize=8)
            ax1.set_ylim(data_ylims)
        if method == "quad":
            for key,fit in self.quad_ABA_fit_dict.items():
                i_c +=1
                dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                v_series = fit["fit_series"]
                ax1.plot(dates,v_series*1000,
                        "-",
                        #markersize=4,
                        linewidth=2,
                        alpha=1,
                        color = f"C{i_c}",
                        label = (f"c0: {fit['popt'][0]:.1e},"
                            + f"c1: {fit['popt'][1]:.1e},"
                            + f"c2: {fit['popt'][2]:.1e},"
                            + f"B_offset: {fit['popt'][3]:.3e}")
                    )
            plt.legend(shadow=True,loc='lower left', bbox_to_anchor=(0, 1),
                  fontsize=14)
        if method == "basic":
            for key,fit in self.basic_ABA_fit_dict.items():
                i_c +=1
                dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                v_series = fit["fit_series"]
                ax1.plot(dates,v_series*1000,
                        "-",
                        #markersize=4,
                        linewidth=2,
                        alpha=1,
                        color = f"C{i_c}",
                        label = (f"c0: {fit['popt'][0]:.1e},"
                            + f"c1: {fit['popt'][1]:.1e},"
                            + f"B_offset: {fit['popt'][2]:.3e}")
                    )
            plt.legend(shadow=True,loc='lower left', bbox_to_anchor=(0, 1))


        plt.xticks(rotation = 45)

        ax1.set_xlabel(r"Time")
        ax1.set_ylabel(r"Voltage [mV]")

        plt.grid(True)
        plt.tight_layout()
        format_im = 'png' #'pdf' or png
        dpi = 300
        if plot_v_eq == True:
            plt.savefig(plot_path + "_veq" + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        else:
            plt.savefig(plot_path + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        ax1.cla()
        return

    def plot_all_ABA_fit_paper(self,
                        plot_path,
                        method = "quad",
                        plot_v_eq = False,
                        **kwargs
                            ):

        fig = plt.figure(0, figsize=(8,6.5))
        ax1=plt.gca()
        # plot data
        for key,slice in self.sliced_dict.items():
            if key == 0:
                label = "data"
                date0 = slice["dates"][0]
            else:
                label = "_nolegend_"
            time_diff_min = [(slice["dates"][i] - date0).total_seconds()/60 
                    for i in range(len(slice["dates"]))]
            ax1.plot(time_diff_min,slice["voltage"]*1000,".",
                markersize=8,
                color = "C0",
                label = label)
        
        data_ylims = ax1.get_ylim()
            
        # plot cropped data
        for key,slice in self.cropped_dict.items():
            if key == 0:
                w_mean = self.fit_results[method]["w_mean"]
                w_std = self.fit_results[method]["w_std"]
                label = ("fit data, Average:" + "\n"
                         + r"$\Delta R = $ "
                         +  f"{w_mean* 1e6:.1f}" 
                         + r"$\pm$" + f"{w_std * 1e6:.1f}" 
                         + r"m$\Omega$"
                        )
            else:
                label = "_nolegend_"
            # ax1.plot(slice["dates"],slice["voltage"]*1000,".",
            #     markersize=8,
            #     color = "C1",
            #     label = label)
            time_diff_min = [(slice["dates"][i] - date0).total_seconds()/60 
                             for i in range(len(slice["dates"]))]
            ax1.plot(time_diff_min,slice["voltage"]*1000,".",
                markersize=8,
                color = "C1",
                label = label)
            

        # plot fits
        #color iterator
        i_c = 1

        if method == "quad":
            for key,fit in self.quad_ABA_fit_dict.items():
                i_c +=1
                # dates = (fit["start_date"]
                #         + np.array(
                #         [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                #         )
                dates = (fit["start_date"]
                        + np.array(
                        [dt.timedelta(seconds = t) for t in fit["t_space"]]) 
                        )
                time_diff_min = [(dates[i] - date0).total_seconds()/60 
                    for i in range(len(dates))]
                v_series = fit["fit_series"]
                ax1.plot(time_diff_min,v_series*1000,
                        "-",
                        #markersize=4,
                        linewidth=2,
                        alpha=1,
                        color = f"C{i_c}",
                        label = (f"{i_c - 1}: "
                                r"$\Delta R$= "
                                 + f"{fit['popt'][3] * 1e6:.1f}"
                                 + r"m$\Omega$")
                    )
            plt.legend(shadow=True,loc='lower left', bbox_to_anchor=(0, 1),
                  fontsize=14,
                  ncol = 2
                  )
        
        plt.xticks(rotation = 45)

        #set x tick spacing to 5 mins
        tick_spacing = 5
        import matplotlib.ticker as ticker
        ax1.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

        ax1.set_xlabel(r"Time [min]")
        ax1.set_ylabel(r"Resistance [$\Omega$]")

        plt.grid(True)
        plt.tight_layout()
        format_im = 'png' #'pdf' or png
        dpi = 300
        if plot_v_eq == True:
            plt.savefig(plot_path + "_veq" + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        else:
            plt.savefig(plot_path + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        ax1.cla()
        return

    def plot_all_quad_ABA_fit(self,
                        plot_path,
                        **kwargs
                            ):

        self.plot_all_basic_ABA_fit(
                        plot_path,
                        method = "quad",
                        **kwargs
                            )
        return

    def weighted_mean_calculation(self,values, errors):
        """
        Values and errors msut be numpy arrays of equivalent shape
        """
        n = len(values)
        weights = 1/ (errors**2)
        # Unweighted
        mean = np.sum(values) / n
        std = np.sqrt(np.sum((values - mean)**2)/(n-1))
        #Weighted
        w_mean = np.sum(weights * values) / np.sum(weights)
        # weighted standat deviation
        w_std = np.sqrt(np.sum(weights * (values - w_mean)**2)
                        / (((n-1)/n) * np.sum(weights)  )
        )
        # the std of mean will liekly be a gross underestimate because it will
        # assume fit errors to be  genuine and independent
        std_of_mean = np.sqrt(1/(np.sum(weights)))
        
        #Sigma Clipping
        for i in range(3): 
            #Weighted
            w_mean_clipped = np.sum(weights * values) / np.sum(weights)
            # weighted standat deviation
            w_std_clipped = np.sqrt(np.sum(weights * (values - w_mean)**2)
                            / (((n-1)/n) * np.sum(weights)  )
            )
            
            residuals = np.abs(values-w_mean_clipped)
            kappa = 2
            mask = residuals < kappa * w_std_clipped
            #clip values and weights
            values = values[mask]
            weights = weights[mask]
            n = len(values)

        #Weighted
        w_mean_clipped = np.sum(weights * values) / np.sum(weights)
        # weighted standat deviation
        w_std_clipped = np.sqrt(np.sum(weights * (values - w_mean)**2)
                        / (((n-1)/n) * np.sum(weights)  )
        )
        

        return (mean, std, w_mean, w_std, std_of_mean,
                w_mean_clipped,w_std_clipped)

    def write_means_to_fit_results(self,
                          dic,
                          ):
        # assumes Measured quantity to be final parameter
        values = np.array([dic[i]["popt"][-1] for i  in dic.keys()])
        errors = np.array([dic[i]["perr"][-1] for i  in dic.keys()])
        res = self.weighted_mean_calculation(values, errors)
        # HACK add c0 parameter (resistance) (for calibration base line)
        c0_values = np.array([dic[i]["popt"][0] for i  in dic.keys()])
        c0_errors = np.array([dic[i]["perr"][0] for i  in dic.keys()])
        res_c0 = self.weighted_mean_calculation(c0_values, c0_errors)
        keys = ["values", "errors", "mean", "std",
                 "w_mean", "w_std", "std_of_mean",
                 "w_mean_clipped","w_std_clipped",
                 "c0_w_mean", "c0_w_std"]
        # TODO Q: shoudl the error not be the std_of_mean
        vals = [values, errors, res[0], res[1], res[2], res[3], res[4],
                res[5], res[6],
                res_c0[2], res_c0[3]]
        print([(keys[i] + ": ", f"{vals[i]:.3e}") for i in range(2,len(keys))])
        self.fit_results[dic[0]["method"]] = {}
        for i,key in enumerate(keys):
            self.fit_results[dic[0]["method"]][key] = vals[i]
        return
    
    def basic_extraction(self):
        self.slice_dict()
        # crop slices
        self.crop_dict()
        self.basic_ABA_fit_all_B()
        # make run_dir
        run_dir = os.path.dirname(self.plot_dir + f"{self.run_name}" 
                                 + os.sep)
        os.makedirs(run_dir, exist_ok=True)
        self.plot_all_basic_ABA_fit(plot_path= run_dir + os.sep +
                                        "basic_ABA_fit_all_B")
        return
    
    def quad_extraction(self):
        self.slice_dict()
        # crop slices
        self.crop_dict()
        self.quad_ABA_fit_all_B()
        # make run_dir
        run_dir = os.path.dirname(self.plot_dir + f"{self.run_name}" 
                                 + os.sep)
        print("cwd: ", os.getcwd())
        print("run_dir: ", run_dir)
        os.makedirs(run_dir, exist_ok=True)
        self.plot_all_quad_ABA_fit(plot_path= run_dir + os.sep +
                                        "quad_ABA_fit_all_B")
        return

    def method_extraction(self, method = "basic", **kwargs):
        self.slice_dict()
        # crop slices
        self.crop_dict()
        # make run_dir
        run_dir = os.path.dirname(self.plot_dir + f"{self.run_name}" 
                                    + os.sep)
        self.run_dir = run_dir
        os.makedirs(run_dir, exist_ok=True)
        if method == "exp_back":
            self.lin_exp_crop_dict(kwargs["front_crop_A"],
                                   kwargs["front_crop_B"],
                                   kwargs["rear_crop"] )
            self.exp_back_ABA_fit_all_B(kwargs["front_crop_A"],
                                   kwargs["front_crop_B"],
                                   kwargs["rear_crop"])
            self.plot_all_basic_ABA_fit(
                        plot_path= run_dir + os.sep +
                                        f"{method}_ABA_fit_all_B",
                        method = method,
                        **kwargs
                        )
        if method == "lin_exp":
            self.lin_exp_crop_dict(kwargs["front_crop_A"],
                                   kwargs["front_crop_B"],
                                   kwargs["rear_crop"] )
            self.lin_exp_ABA_fit_all_B(kwargs["front_crop_A"],
                                   kwargs["front_crop_B"],
                                   kwargs["rear_crop"])
            self.plot_all_basic_ABA_fit(
                        plot_path= run_dir + os.sep +
                                        f"{method}_ABA_fit_all_B",
                        method = method,
                        **kwargs
                        )
        if method == "basic":
            self.basic_ABA_fit_all_B()
            self.plot_all_basic_ABA_fit(plot_path= run_dir + os.sep +
                                        "basic_ABA_fit_all_B")
        if method == "quad":
            self.quad_ABA_fit_all_B()
            self.plot_all_quad_ABA_fit(plot_path= run_dir + os.sep +
                                            "quad_ABA_fit_all_B")
        return
#############
######## END Extractor

def make_tex_table(
                       col_headers,
                       row_headers,
                       data_array,
                      ):
        """
        converts data array into tex table

        col_headers :  (n,) list of strings
            column headers
        row_headers :  (m,) list of strings
            row headers
        data_array : (m,n) list of strings
        """
        n = len(col_headers)
        m = len(row_headers)
        begin_tab = (r"\begin{tabular}" + r"{|c|" + n*r"c|" + "}" + r"\hline"
                    + "\n")
        end_tab = r"\hline" + "\n" + r"\end{tabular}"

        # make colum header line by concatenating headers with spacers
        col_head_tab = " & "
        for i,header in enumerate(col_headers):
            col_head_tab += header
            if i < n -1:
                col_head_tab += r" & "
            else:
                # add line end
                col_head_tab += r" \\ \hline " + "\n"

        # make table body by inserting row header and then concatenating
        # data array entries
        body_tab = ""
        for i, header  in enumerate(row_headers):
            body_tab += header
            body_tab += r" & "
            for j in range(n):
                body_tab += data_array[i][j]
                if j < m -1:
                    body_tab += r" & "
                else:
                    # add line end
                    body_tab += r" \\ " + "\n"

        tab_string = begin_tab + col_head_tab+ body_tab + end_tab
        return tab_string

def fudge_power_plot(data_sets, plot_path, T_lst= [295,1310,2350],
                     flow_list = [0.2, 1.0, 20.0], i_split = 17):

    # Run once clipped once unclipped see 
    for clip in [False, True]:
         #setup plot arrays
        v_mean_arr =np.transpose(
                np.array([[ext.fit_results["quad"]["w_mean"] 
                # if isinstance(ext, Extractor)
                # else "N.A."
                for ext in row
                ] 
                for row in data_sets
                ]))
        v_err_arr =np.transpose(
                np.array([[ext.fit_results["quad"]["w_std"] 
                # if isinstance(ext, Extractor)
                # else "N.A."
                for ext in row
                ] 
                for row in data_sets
                ]))
        if clip:
            #setup plot arrays
            v_mean_arr =np.transpose(
                    np.array([[ext.fit_results["quad"]["w_mean_clipped"] 
                    # if isinstance(ext, Extractor)
                    # else "N.A."
                    for ext in row
                    ] 
                    for row in data_sets
                    ]))
            v_err_arr =np.transpose(
                    np.array([[ext.fit_results["quad"]["w_std_clipped"] 
                    # if isinstance(ext, Extractor)
                    # else "N.A."
                    for ext in row
                    ] 
                    for row in data_sets
                    ]))

        print("v_mean_arr:", v_mean_arr)
        #fudge the conversion (-1 for additional self heating ddue to meas current)
        µW_per_ohm = 7.4 * 1000 - 1 * 1000 
        p_arr= µW_per_ohm * v_mean_arr
        p_err_arr = µW_per_ohm * v_err_arr
        T_lst_plot = T_lst[0:len(p_arr[0])]
        print("p_arr",p_arr)
        print("p_arr_err",p_err_arr)
        print("len(T_lst_plot)",len(T_lst_plot))
        print("len p_arr",len(p_arr[0]))
        #plot
        fig = plt.figure(0, figsize=(8,6.5))
        ax1=plt.gca()
        for i,flow_const_list in enumerate(flow_list):
            #bodge different color on the way down
            if True:
                i_split = i_split
                ax1.errorbar(T_lst_plot[:i_split], p_arr[i][:i_split],
                    yerr = p_err_arr[i][:i_split],
                    label = f"{flow_list[i]} sccm, up",
                    color = "C0",
                    linewidth = 0.5,
                    elinewidth= 4  
                    )
                ax1.errorbar(T_lst_plot[i_split:], p_arr[i][i_split:], 
                    yerr = p_err_arr[i][i_split:],
                    label = f"{flow_list[i]} sccm, down",
                    color = "C1",
                    linewidth = 0.5,
                    elinewidth= 1.5  
                    )
                continue
            else:
                # Original version without splitting
                ax1.errorbar(T_lst_plot, p_arr[i], yerr = p_err_arr[i],
                            label = f"{flow_list[i]} sccm",
                            linewidth = 0.5,
                            elinewidth= 3  
                                )


        ax1.set_xlabel(r"Source Temperature [K]")
        ax1.set_ylabel(r"Detected Power [µW]")

        plt.grid(True)
        plt.legend(shadow=True)
        plt.tight_layout()
        format_im = 'png' #'pdf' or png
        dpi = 600
        if clip:
            plt.savefig(plot_path
                        + '_clipped.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        else:
            plt.savefig(plot_path
                        + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        ax1.cla()
    return



def fudge_power_plot_fit(data_sets, plot_path, T_lst= [295,1310,2350],
                     flow_list = [0.2, 1.0, 20.0], i_split = 17):

    # Run once clipped once unclipped see 
    #for clip in [False, True]:
    for clip in [False]:
         #setup plot arrays
        v_mean_arr =np.transpose(
                np.array([[ext.fit_results["quad"]["w_mean"] 
                # if isinstance(ext, Extractor)
                # else "N.A."
                for ext in row
                ] 
                for row in data_sets
                ]))
        v_err_arr =np.transpose(
                np.array([[ext.fit_results["quad"]["w_std"] 
                # if isinstance(ext, Extractor)
                # else "N.A."
                for ext in row
                ] 
                for row in data_sets
                ]))
        if clip:
            #setup plot arrays
            v_mean_arr =np.transpose(
                    np.array([[ext.fit_results["quad"]["w_mean_clipped"] 
                    # if isinstance(ext, Extractor)
                    # else "N.A."
                    for ext in row
                    ] 
                    for row in data_sets
                    ]))
            v_err_arr =np.transpose(
                    np.array([[ext.fit_results["quad"]["w_std_clipped"] 
                    # if isinstance(ext, Extractor)
                    # else "N.A."
                    for ext in row
                    ] 
                    for row in data_sets
                    ]))

        #fudge the conversion (-1 for additional self heating ddue to meas current)
        µW_per_ohm = 7.4 * 1000 - 1 * 1000 
        p_arr, p_err_arr = µW_per_ohm * v_mean_arr, µW_per_ohm * v_err_arr
        T_lst_plot = T_lst[0:len(p_arr[0])]
        #print("p_arr",p_arr)
        #print("len(T_lst_plot)",len(T_lst_plot))
       #print("len p_arr",len(p_arr[0]))

        # # fit all powers below 1500K
        # T_mask = np.array(T_lst_plot)<1500
        # fit all powers below 1300K
        T_mask = np.array(T_lst_plot)<1300
        #print("T_mask", T_mask)
        T_masked = np.array(T_lst_plot)[T_mask]
        p_masked,p_err_masked = p_arr[0][T_mask], p_err_arr[0][T_mask]
        def fit_func(x,m,b):
            return m*x + b
        popt, pcov = curve_fit(fit_func, T_masked, p_masked
                                        # p0 = [80, 0
                                        #         , 1e-5, 1e-5
                                        #         , 80,200
                                        #         , 1e-4]
                                        # , bounds = (
                                        #         [-np.inf, -2* 1e-8
                                        #         , 1e-8, 1e-8
                                        #         , 1, 1
                                        #         , -np.inf],
                                        #         [ np.inf, +2* 1e-8
                                        #         , 1e-3, 1e-3
                                        #         , np.inf, np.inf
                                        #         , np.inf]
                                        # )
                                        ,sigma = p_err_masked
                                        , absolute_sigma = True
                                        )
        T_linspace = np.linspace(np.min(T_lst_plot),np.max(T_lst_plot))
        p_fit = [fit_func(T,popt[0],popt[1]) for T in T_linspace]
        #print("p_fit",p_fit)

        p_gas = np.array([fit_func(T,popt[0],popt[1]) for T in T_lst_plot])
        p_residuals = p_arr - p_gas
        T_room = 296 # should really be T_chamber
        p_thermal_beam_gas = p_gas - fit_func(T_room ,popt[0],popt[1])
        E_rec_per_E_thermal_factor = (4.4634/2 # Ev per recomb atom
                                    /(1/(1.16045e4) # EV per kelvin
                                    * (np.array(T_lst_plot) - T_room))
                                        )
        #print("E_rec_per_E_thermal_factor", E_rec_per_E_thermal_factor)

        # #cracking_eff_bodge THIS ONE IS WRONG:
        # p_theory = 6.8 #µW prediction
        # cracking_eff_bodge = p_residuals / p_theory
        # cracking_eff_bodge_err =  p_err_arr / p_theory

        
        #cracking_eff_bodge comparison to Thermal gas method
        # assumes eta_rec and eta_thermal to be the same
        # E_rec_per_E_thermal_factor = 14.1 # BODGE Bodge 
        # (only for deltaT = 1950)
        cracking_eff_bodge = ((p_residuals / p_thermal_beam_gas)/
                                E_rec_per_E_thermal_factor)
        cracking_eff_bodge_err =  ((p_err_arr / p_thermal_beam_gas)/
                                E_rec_per_E_thermal_factor)

        #plot
        fig = plt.figure(0, figsize=(8,6.5))
        ax1=plt.gca()
        for i,flow_const_list in enumerate(flow_list):
            ax1.plot(T_linspace, p_fit, color = "C4", ls = "-"
                    ,label = "P_gas_fit" 
                    )
            #bodge different color on the way down
            if True:
                ax1.errorbar(T_lst_plot[:i_split], p_arr[i][:i_split],
                    yerr = p_err_arr[i][:i_split],
                    label = f"{flow_list[i]} sccm, up",
                    color = "C0",
                    linewidth = 0.5,
                    elinewidth= 4  
                    )
                ax1.errorbar(T_lst_plot[i_split:], p_arr[i][i_split:], 
                    yerr = p_err_arr[i][i_split:],
                    label = f"{flow_list[i]} sccm, down",
                    color = "C1",
                    linewidth = 0.5,
                    elinewidth= 1.5  
                    )
                continue
            else:
                # Original version without splitting
                ax1.errorbar(T_lst_plot, p_arr[i], yerr = p_err_arr[i],
                            label = f"{flow_list[i]} sccm",
                            linewidth = 0.5,
                            elinewidth= 3  
                                )


        ax1.set_xlabel(r"Source Temperature [K]")
        ax1.set_ylabel(r"Detected Power [µW]")

        plt.grid(True)
        plt.legend(shadow=True)
        plt.tight_layout()
        format_im = 'png' #'pdf' or png
        dpi = 600
        if clip:
            plt.savefig(plot_path
                        + '_clipped.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        else:
            plt.savefig(plot_path
                        + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        ax1.cla()


        #plot cracking eff_bodge
        fig = plt.figure(0, figsize=(8,6.5))
        ax1=plt.gca()
        for i,flow_const_list in enumerate(flow_list):
            #bodge different color on the way down
            # Don't  plot 0 temperature point in cracking eff, becasue it is
            # not well constrained
            cut = 0
            if True:
                ax1.errorbar(T_lst_plot[cut:i_split], 
                    cracking_eff_bodge[i][cut:i_split],
                    yerr = cracking_eff_bodge_err[i][cut:i_split],
                    label = f"{flow_list[i]} sccm, up",
                    color = "C0",
                    linewidth = 0.5,
                    elinewidth= 4  
                    )
                ax1.errorbar(T_lst_plot[i_split:], 
                    cracking_eff_bodge[i][i_split:], 
                    yerr = cracking_eff_bodge_err[i][i_split:],
                    label = f"{flow_list[i]} sccm, down",
                    color = "C1",
                    linewidth = 0.5,
                    elinewidth= 1.5  
                    )
                continue
            else:
                # Original version without splitting
                ax1.errorbar(T_lst_plot, p_arr[i], yerr = p_err_arr[i],
                            label = f"{flow_list[i]} sccm",
                            linewidth = 0.5,
                            elinewidth= 3  
                                )


        ax1.set_xlabel(r"Source Temperature [K]")
        ax1.set_ylabel(r"Cracking efficiency bodge")

        plt.grid(True)
        plt.legend(shadow=True)
        plt.tight_layout()
        format_im = 'png' #'pdf' or png
        dpi = 600
        if clip:
            plt.savefig(plot_path + "_CEB"
                        + '_clipped.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        else:
            plt.savefig(plot_path + "_CEB"
                        + '.{}'.format(format_im),
                        format=format_im, dpi=dpi)
        ax1.cla()
        

        #### HUGE HACK ###
        (directory , filename) = os.path.split(plot_path) 
        print("directory", directory)
        # Save raw power
        p_arr, p_err_arr
        power_dict = {"T": T_lst_plot,
                "power": p_arr[0],
                "power_err": p_err_arr[0]
                }
        #print("power_dict", power_dict)
        save_dict(power_dict, directory + os.sep +  "power_dict")
        # Save Subtracted power
        p_arr, p_err_arr
        power_recomb_dict = {"T": T_lst_plot,
                "power_recomb": p_residuals[0],
                "power_err": p_err_arr[0]
                }
        #print("power_recomb_dict", power_recomb_dict)
        save_dict(power_recomb_dict, directory + os.sep +"power_recomb_dict")
        # save cracking efficiency list:
        CEB_dict = {"T": T_lst_plot,
                "cracking_eff_bodge": cracking_eff_bodge[0],
                "cracking_eff_bodge_err": cracking_eff_bodge_err[0]
                }
        #print("CEB_dict", CEB_dict)
        save_dict(CEB_dict, directory + os.sep +"CEB_dict")



    return (T_lst_plot, cracking_eff_bodge,
             cracking_eff_bodge_err)

# Out of date see fudge_power_plot for working version
def fudge_power_per_sccm_plot(data_sets, plot_path, T_lst= [295,1310,2350],
                     flow_list = [0.2, 1.0, 20.0]):

    #setup plot arrays
    v_mean_arr =np.transpose(np.array([[ext.fit_results["quad"]["w_mean"] 
             if isinstance(ext, Extractor)
             else "N.A."
             for ext in row
             ] 
             for row in data_sets
            ]))
    v_err_arr =np.transpose(np.array([[ext.fit_results["quad"]["w_std"] 
             if isinstance(ext, Extractor)
             else "N.A."
             for ext in row
             ] 
             for row in data_sets
            ]))

    #fudge the conversion # even more fudged
    µW_per_ohm = 7 * 1000 *1000
    p_arr, p_err_arr = µW_per_ohm * v_mean_arr, µW_per_ohm * v_err_arr
    T_lst = T_lst[0:len(p_arr)]
    #plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    for i,flow_const_list in enumerate(flow_list):
        ax1.errorbar(T_lst, p_arr[i]/flow_list[i], yerr = p_err_arr[i],
                     label = f"{flow_list[i]} sccm",
                     linewidth = 0.5,
                     elinewidth= 3  
                        )


    ax1.set_xlabel(r"Source Temperature [K]")
    ax1.set_ylabel(r"Detected Power [nW/sccm]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 600
    plt.savefig(plot_path
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    return


###########
############defudging
def calib_factor(r_base):
    # HACK HACK HACK
    #popt from "calib_analysis" based on 2023-01-09_long_calib
    popt = [-4.02271476e-10,  8.90056085e-08, -1.09041548e-06,  1.37890429e-01,
    6.70824014e+01]

    r_from_p = lambda p: np.poly1d([*popt])(p) # p un µW
    k_from_deriv = lambda p: 1/(np.poly1d([*popt]).deriv()(p))
    # Goal is to make an inverse function that yields p_from_r
    from scipy.interpolate import interp1d
    p_space = np.linspace(0,300,301)
    r_space = r_from_p(p_space)

    p_from_r = lambda r: interp1d(r_space,p_space)(r)
    k = k_from_deriv(p_from_r(r_base))
    # reduce k by 1 µW/Ohm to eliminate self heating factor
    calib = k - 1
    return calib

def make_result_dict(ext_dict):
        v_mean_arr = np.array(
                    [val["extractor"].fit_results["quad"]["w_mean"] 
                    for val in ext_dict.values()]
                    )
        v_err_arr = np.array(
                    [val["extractor"].fit_results["quad"]["w_std"] 
                    for val in ext_dict.values()]
                    )
        # #print("v_mean_arr", v_mean_arr)
        # µW_per_ohm = 7.4 * 1000 - 1 * 1000
        # HACK 2024-04-16 Introduce recalibration based on base temperature of
        # wire
        r0_arr = np.array(
                    [val["extractor"].fit_results["quad"]["c0_w_mean"] 
                    for val in ext_dict.values()]
                    ) * 1000
        µW_per_ohm = calib_factor(r0_arr)
        p_arr = µW_per_ohm * 1000 * v_mean_arr
        p_err_arr = µW_per_ohm * 1000* v_err_arr
        result_dict = {}
        # result_dict["µW_per_ohm"] = np.array([µW_per_ohm 
        #                                  for i in range(len(v_mean_arr))])
        result_dict["µW_per_ohm"] = µW_per_ohm
        result_dict["v_mean_arr"] = v_mean_arr
        result_dict["v_err_arr"] = v_err_arr
        result_dict["p_arr"] = p_arr
        result_dict["p_err_arr"] = p_err_arr
        return result_dict

def power_plot(result_dict, plot_path,
                x_lst= [295,1310,2350],
                x_label = r"Source Temperature [K]",
                y_label = r"Detected Power [µW]",
                flow_list = [1.0],
                i_split = 17,):

    p_arr = result_dict["p_arr"]
    p_err_arr = result_dict["p_err_arr"]

    sort_args = np.argsort(x_lst)
    x_lst = x_lst
    x_lst_plot = x_lst[0:len(p_arr)]
    #plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    for i,flow_const_list in enumerate(flow_list):
        #bodge different color on the way down
        if True:
            i_split = i_split
            ax1.errorbar(x_lst_plot[:i_split], p_arr[:i_split],
                yerr = p_err_arr[:i_split],
                label = f"{flow_list[i]} sccm, up",
                color = f"C{2*(i//2)}",
                linewidth = 0.5,
                elinewidth= 4  
                )
            if len(x_lst_plot) > i_split:
                ax1.errorbar(x_lst_plot[i_split:], p_arr[i_split:], 
                    yerr = p_err_arr[i_split:],
                    label = f"{flow_list[i]} sccm, down",
                    color = f"C{2*(i//2) + 1}",
                    linewidth = 0.5,
                    elinewidth= 1.5  
                    )
            continue
        else:
            # Original version without splitting
            ax1.errorbar(T_lst_plot, p_arr[i], yerr = p_err_arr[i],
                        label = f"{flow_list[i]} sccm",
                        linewidth = 0.5,
                        elinewidth= 3  
                            )


    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 600
    plt.savefig(plot_path
                    + '.{}'.format(format_im),
                    format=format_im, dpi=dpi)
    ax1.cla()
    return

def power_plot_mult(results_array, plot_path,
                x_arr= [[295,1310,2350]],
                x_label = r"Source Temperature [K]",
                y_label = r"Detected Power [µW]",
                label_list = ["1.0 sccm"],
                i_split = 17,):
    #Results must be presented in format:
    # [
    #     [y_results_1 ,y_err_1], [y_results_2 ,y_err_2]
    # ]
    
    #plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    for i,label in enumerate(label_list):
        p_arr = results_array[i][0]
        p_err_arr = results_array[i][1]
        x_lst_plot = x_arr[i][0:len(p_arr)]
        #bodge different color on the way down
        ewidth = len(label_list) + 2 - i 
        ax1.errorbar(x_lst_plot, p_arr,
            yerr = p_err_arr,
            label = label,
            color = f"C{i}",
            linewidth = 0.5,
            elinewidth= ewidth
        )

    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label)

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 600
    plt.savefig(plot_path
                    + '.{}'.format(format_im),
                    format=format_im, dpi=dpi)
    ax1.cla()
    return

def sort_by_z_list(z_list, result_dict):
    arg_order = np.argsort(z_list)
    z_list_sorted = np.array(z_list)[arg_order]
    result_dict_sorted = {key:value[arg_order] 
                    for key,value in result_dict.items()}
    return z_list_sorted, result_dict_sorted

def T_from_TC(TC_val):
    # Bodge of  Temperature to TC PID value alignnment
    # https://discourse.project8.org/t/
    # mainz-habs-power-supply-tdk-lambda/291
    # Digitized 2 values along the fit from the plot
    TC_lst = [52.395, 710.479]
    T_Lst = [310.810, 2185.810]
    f = lambda x, a,b: a*x + b
    popt, pcov = curve_fit(
        f,xdata = TC_lst,ydata = [310.810, 2185.810])
    return f(TC_val, *popt)

def extract(sc_dir, HABS_dir, run_name, data_name,
                pid_list, label_list,
                global_start_date,
                n_mins, flow_resynch, flow_interval, 
                HABS_resynch, HABS_interval, HABS_front_cut
                ):
        #### HACK START Comment oot to laod instead of calculating
        #load wire data
        prep_data_slowdash("../SC_downloads/Wire/"+ data_name + ".json"
                        , run_name
                        )
        data_dict = load_data(data_name)

        # extraction dict:
        ext_dict = {}
        # fudge by cutting malfunctioning last current measurement
        for i,label in enumerate(label_list):
            start_time = time.time()
            if i == 0:
                ext = Extractor(data_dict, 
                    start_date = global_start_date + i  * HABS_interval
                                + flow_interval * 4, 
                                # cut  first 20 minutes that were not logged 
                    end_date = global_start_date + (i+1) * HABS_interval
                            - flow_interval ,

                    time_interval = flow_interval,
                    run_name = run_name + os.sep 
                            + "{}_{}".format(i,label),
                    utc_offset = 2,
                    initial_state = 0,
                    front_crop = dt.timedelta(minutes = 2,seconds  =30),
                    rear_crop = dt.timedelta(minutes  =0, seconds  =45)
                    )
            else:
                ext = Extractor(data_dict, 
                    start_date = global_start_date + i * HABS_interval
                                + HABS_front_cut, 
                    end_date = global_start_date + (i+1) * HABS_interval
                            - flow_interval ,
                            # cut last extraneous "B" cycle since we only care 
                            # fully completed ABA cycles
                    time_interval = flow_interval,
                    run_name = run_name + os.sep 
                            + "{}_{}".format(i,label),
                    utc_offset = 2,
                    initial_state = 0,
                    front_crop = dt.timedelta(minutes = 2,seconds  =30),
                    rear_crop = dt.timedelta(minutes  =0, seconds  =45)
                    )
            ext.basic_extraction()
            ext.quad_extraction()
            ext_dict[i] = {"label":label,
                        "extractor":ext
                            }
            run_time = (time.time()- start_time)/60
            print("time taken for step:", run_time, "min" )
            start_time = time.time()

        save_dict(ext_dict, self.plot_dir + run_name + os.sep 
                                    + "ext_dict")
        data_sets = [[ext_dict[i]["extractor"]] 
                    for i,z in enumerate(label_list)
                    ]
        save_dict(data_sets, self.plot_dir + run_name + os.sep 
                                    + "data_sets")
        return ext
        # ##### HACK END  load  instead        

####################################################################
if __name__ =="__main__": 
####################################################################
# Defudgin testbed: TODO Defudge
# Left as example
####################### 2023-06-30
# 2023-06-28_005sccm_TC_temp_cycle_jf+hg_wire
    #####  Power load tests
    sc_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep + "../" 
         + "SC_downloads/")
    HABS_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep + "../" 
         + "SC_downloads/HABS_power/")
    os.makedirs(HABS_dir, exist_ok=True)
    run_name = "2023-06-28_005sccm_TC_temp_cycle_jf+hg_wire"
    data_name = run_name    

    # # flow on-off cycles
    # # TODO reformat this into a function with save and load
    # # list of current settings
    pid_list = ([ (715/10) * ((i))  for i in range(0,11,1)]
            + [(715/10) * (10 - ((i)))  for i in range(0,11,1)]
            )
    label_list = pid_list
    print("label_list", label_list)

    global_start_date = dt.datetime(2023, 6, 28, 23, 37, 20,
                                tzinfo=dt.timezone(dt.timedelta(hours=1)))
    # Second annd milisecond fudges here are to realign unaccounted slow
    # controls delay
    #     #minutes per day shift
    n_mins = 5
    flow_resynch = 780 + ((n_mins/(24*60))*(5*60)*1000) 
    flow_interval = dt.timedelta(minutes = 5, milliseconds=flow_resynch)
    HABS_resynch =  flow_resynch * 1e-3 *24
    HABS_interval = dt.timedelta(minutes = 120, seconds=HABS_resynch )
    # cut first 10 minutes both to eliminate synchronization drift and the
    # fastest part of the tempperature change
    HABS_front_cut = dt.timedelta(minutes = 40)

    # ########
    # #######Comment oout if extraction has already run
    # extract(sc_dir, HABS_dir, run_name, data_name,
    #             pid_list, label_list,
    #             global_start_date,
    #             n_mins, flow_resynch, flow_interval, 
    #             HABS_resynch, HABS_interval, HABS_front_cut
    #             )

    # ######## Load data back             
    data_sets = load_dict(self.plot_dir + run_name + os.sep 
                                + "data_sets")
    # print("data_sets:",data_sets)

    ext_dict = load_dict(self.plot_dir + run_name + os.sep 
                                + "ext_dict")
    result_dict = make_result_dict(ext_dict)


    power_plot(result_dict,
                plot_path = self.plot_dir + run_name + os.sep 
                           + "Detected_P_over_pid",
                x_lst= label_list,
                x_label = r"PID Setpoint",
                y_label = r"Detected Power [µW]",
                flow_list = [0.05],
                i_split = 11,)
    # fudge_power_plot(data_sets, 
    #             plot_path= self.plot_dir + run_name + os.sep 
    #                        + "Detected_P_over_T{}".format(
    #                         denoise_string),
    #             flow_list = [10],
    #             T_lst = T_list
    #             )

    T_list = T_from_TC(np.array(pid_list))
    
    fudge_power_plot_fit(data_sets, 
                plot_path= self.plot_dir + run_name + os.sep 
                           + "Detected_P_over_T_fit",
                flow_list = [0.05],
                T_lst = T_list,
                i_split = 11,
                )

