# based on "Voltage_stability.py'

########## imports and setup
from calendar import day_abbr
from matplotlib.cbook import flatten
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from datetime import datetime, tzinfo
import datetime as dt
from scipy.interpolate import interp1d
import dill
import json
from decimal import Decimal
def remove_exponent(d):
    return d.quantize(Decimal(1)) if d == d.to_integral() else d.normalize()

#import pressure_analysis as pa

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)

#######################

# Plot with movinng average
def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w

def plot_with_avg(filename, plotname, plot_dir):
    data = np.loadtxt(filename, delimiter= ",")
    t_series, v_series = data.T
    mavg_len = 30 # In seconds
    w = mavg_len * 1 # ~1 measurements per second
    t_avg_series = moving_average(t_series, w)
    v_avg_series = moving_average(v_series, w)

    try:
            sigma= np.std(v_avg_series - v_series[w//2:-w//2])
    except:
            sigma= np.std(v_avg_series[:-1] - v_series[w//2:-w//2]) 

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.plot(t_series-t_series[0],v_series*1000,".",markersize=1,
            label = f"data, std vs mavg = {sigma * 1000:.1e}")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    ax1.set_xlabel(r"Time [s]")
    ax1.set_ylabel(r"Voltage [mV]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def plot_with_dates(filename,data_dir, plotname,plot_dir, start = 0, end = -1):
    if type(filename) == type("string"):
        file, ext = os.path.splitext(filename)
        if ext == ".txt":
            data = np.loadtxt(filename, delimiter= ",")
            t_series, v_series = data[start:end].T
            dates = mpl.dates.num2date(t_series/86400)
        if ext == ".pkl":
            data_dict = load_data(filename, data_dir)
            v_series = data_dict["voltage"]
            dates = data_dict["dates"]
        else:
            Exception("Invalid file extension. .txt or .pkl are valid")
    if type(filename) == type({"dict":[]}):
        data_dict = filename
        v_series = data_dict["voltage"][start:end]
        dates = data_dict["dates"][start:end]

    mavg_len = 30 # In seconds
    w = mavg_len * 1 # ~1 measurements per second
    #t_avg_series = moving_average(t_series, w)
    v_avg_series = moving_average(v_series, w)

    try:
            sigma= np.std(v_avg_series - v_series[w//2:-w//2])
    except:
            sigma= np.std(v_avg_series[:-1] - v_series[w//2:-w//2]) 

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.plot(dates,v_series*1000,".",markersize=1,
            label = f"data, std vs mavg = {sigma * 1000:.1e}")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"Voltage [mV]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_dates"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def select_date_indices(data_dict,start_date , end_date ):
    # # original: slow?
    # start  = np.argwhere(start_date < data_dict["dates"])[0][0]
    # end  = np.argwhere(data_dict["dates"] < end_date)[-1][0]

    # Faster? ~ factor 4 in exammple tested. Assumes ordered date list
    # look for first instance when date after start date
    for i,date in enumerate(data_dict["dates"]):
        if start_date > date:
            start = i
            continue
        else:
            start = i
            break
    # to look for end start from "start"
    # look for first instance when date after end date
    for i,date in enumerate(data_dict["dates"][start:]):
        if end_date > date:
            end = start + i
            continue
        else:
            end = start + i
            break

    return start, end

def plot_with_dates_select(data_dict, plotname, plot_dir,
                           start_date, end_date,
                           utc_offset = 1,
                           large_points = False,
                           defuzz = False,
                           ):
    if type(data_dict) == type({"dict":[]}):
        start,end = select_date_indices(data_dict,start_date, end_date)
        v_series = data_dict["voltage"][start:end]
        dates = data_dict["dates"][start:end]
        # move to German timezone
        dates = np.array([
                date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                for date in dates])
            
    if defuzz == True:
        # Remove downward dip systematic noise pattern I call "fuzz"
        # If point is lower than both its neighbors by an average
        # of more than 2µV: (arbitrary cut)
        # Cut point before and 2 after 
        # Will result in ~50% loss of data
        dip_mins = [i for i in range(1,len(v_series)-1)
                     if v_series[i] < (v_series[i-1] + v_series[i+1])/2 - 2e-6
                        ]
        print(dip_mins)
        fuzz_mask = [[i-1, i, i +1, i+2] 
                     if i+2 < len(v_series) 
                     else [i-1, i, i +1]
                     for i in dip_mins
                    ]
        # flatten fuzz mask and remove duplicates
        fuzz_mask = list(set(np.concatenate(fuzz_mask)))
        # invert fuzz mask
        fuzz_mask = [i for i in range(1, len(v_series)-1) 
                     if i not in fuzz_mask]
        print(fuzz_mask)



    else:
        Exception("Invalid file provided")

    mavg_len = 30 # In seconds
    w = mavg_len * 1 # ~1 measurements per second
    #t_avg_series = moving_average(t_series, w)
    v_avg_series = moving_average(v_series, w)

    try:
            sigma= np.std(v_avg_series - v_series[w//2:-w//2])
    except:
            sigma= np.std(v_avg_series[:-1] - v_series[w//2:-w//2]) 
        

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    if large_points == False:
        ax1.plot(dates,v_series*1000,".",markersize=1,
                label = f"data, std vs mavg = {sigma * 1000:.1e}")
    elif large_points == True:
        ax1.plot(dates,v_series*1000,"x",markersize=3,
                 ls = "-", linewidth = 0.5,
                label = f"data, std vs mavg = {sigma * 1000:.1e}")
    if defuzz == True:
        ax1.plot(dates[fuzz_mask],v_series[fuzz_mask]*1000,"x",markersize=3,
                        ls = "-", linewidth = 0.5,
                label = f"defuzzed data")
        
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"Voltage [mV]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_dates"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

#from https://www.kongsberg.com/globalassets/maritime/km-products/
# product-documents/tsiec751_ce.pdf
def T_to_R(T):
    """
    T in °C
    """
    T_list = np.linspace(-200,850, num = 1051)
    A = 3.9083e-3
    B = -5.775e-7
    C = -4.183e-12
    R_0 = 1000.0
    def fun(T):
        if T >= 0:
            R_T = R_0 * (1 + A*T + B*T**2)
        if T < 0:
            R_T = R_0 * (1 + A*T + B*T**2 + C*(T - 100)* T**3) 
        return R_T
    T_list = np.linspace(-200,850, num = 1051)
    R_list = np.vectorize(fun)(T_list)
    f_int = interp1d(T_list,
                    R_list
                    ,
                    kind = "cubic"
                    ,fill_value="extrapolate"
                    )
    return f_int(T)

def R_to_T(R):
    T_list = np.linspace(-200,850, num = 1051)
    R_list = np.vectorize(T_to_R)(T_list)
    f_int = interp1d(R_list,
                    T_list
                    ,
                    kind = "cubic"
                    ,fill_value="extrapolate"
                    )
    return f_int(R)

def plot_Pt_1000_select(data_dict, plotname,plot_dir, start_date, end_date,
                           utc_offset = 1
                           ):
    if type(data_dict) == type({"dict":[]}):
        start,end = select_date_indices(data_dict,start_date, end_date)
        dates = data_dict["dates"][start:end]
        # move to German timezone
        dates = np.array([
                date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                for date in dates])
            
    else:
        Exception("Invalid file provided")

    fig = plt.figure(0, figsize=(8,6.5))
    fig.set_size_inches(10.5, 6)
    ax1=plt.gca()
    ax1.plot(dates,data_dict["R_Pt_1000"][start:end],".",markersize=1,
            label = f"R_Pt_1000")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    secax_1 = ax1.secondary_yaxis(-0.12, functions=(R_to_T
                                , T_to_R))
    secax_1.set_ylabel('T [°C]')

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_Pt_1000"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def prep_data(filename, dataname, data_dir,start = 0, end = -1):
    data = np.loadtxt(filename, delimiter= ",")
    t_series, v_series = data[start:end].T[0:2] # select only entry 0 and 1 
    # for every datapoint
    if len(data[0]) == 3:
            t_series, v_series, r_series = data[start:end].T[0:3]
    dates = mpl.dates.num2date(t_series/86400)
    data_dict = {"dates":np.array(dates),
                "voltage":np.array(v_series)}
    #extract PT_1000 resistance if  it was recorded
    if len(data[0]) == 3:
        _, _, r_series = data[start:end].T[0:3]
        data_dict["R_Pt_1000"] = np.array(r_series)

    save_data(dataname, data_dir, data_dict)
    return data_dict

def prep_data_with_i(filename, dataname, data_dir,start = 0, end = -1):
    data = np.loadtxt(filename, delimiter= ",")
    t_series, v_series, i_series, i_set = data[start:end].T
    dates = mpl.dates.num2date(t_series/86400)
    data_dict = {"dates":dates, "voltage":v_series,
                 "i_series":i_series, "i_set":i_set}
    save_data(dataname + "_with_i", data_dict)
    return data_dict

def prep_data_slowdash(filename, dataname, data_dir,start = 0, end = -1, wire_str = "Wire1"):
    # This function loads data files in the slowdash json format, and converts
    # them to the same pickled data format as previous data
    with open(filename,"r") as f:
        data = json.load(f)
    t_series = (np.array(data[f"read_V_{wire_str}_WireDet_Source_MATS"]["t"]) 
               + float(data[f"read_V_{wire_str}_WireDet_Source_MATS"]["start"]))
    v_series = np.array(data[f"read_V_{wire_str}_WireDet_Source_MATS"]["x"]) 
    # sample board resistance on the same timeline
    def R_interpolate(data):
        t_arr = (np.array(data["read_ohm_Pt1000_WireDet_Source_MATS"]["t"]) 
               + float(data["read_ohm_Pt1000_WireDet_Source_MATS"]["start"]))
        f_int = interp1d(t_arr,
                        data["read_ohm_Pt1000_WireDet_Source_MATS"]["x"],
                        kind = "linear",fill_value=(np.NaN,np.NaN),
                        bounds_error= False
                )
        return f_int
    r_int = R_interpolate(data)
    r_series = r_int(t_series)
    dates = mpl.dates.num2date(t_series/86400)
    data_dict = {"dates":np.array(dates),
                "voltage":np.array(v_series)
                }
    data_dict["R_Pt_1000"] = np.array(r_series)

    save_data(dataname, data_dir, data_dict)
    return data_dict

def prep_data_calib(filename, dataname, data_dir, start = 0, end = -1):
    # This function loads data files in the slowdash json format, and converts
    # them to the same pickled data format as previous data
    with open(data_dir + filename,"r") as f:
        data = json.load(f)
    try:
        sep_str = "__"
        t_series = (np.array(data[f"read_V_Wire1_WireDet{sep_str}Source_MATS"]["t"]) 
                + float(data[f"read_V_Wire1_WireDet{sep_str}Source_MATS"]["start"]))
    except:
        sep_str = "_"
        t_series = (np.array(data[f"read_V_Wire1_WireDet{sep_str}Source_MATS"]["t"]) 
               + float(data[f"read_V_Wire1_WireDet{sep_str}Source_MATS"]["start"]))
    v_series = np.array(data[f"read_V_Wire1_WireDet{sep_str}Source_MATS"]["x"]) 
    # sample board resistance on the same timeline
    def R_interpolate(data):
        t_arr = (np.array(data[f"read_ohm_Pt1000_WireDet{sep_str}Source_MATS"]["t"]) 
               + float(data[f"read_ohm_Pt1000_WireDet{sep_str}Source_MATS"]["start"]))
        f_int = interp1d(t_arr,
                        data[f"read_ohm_Pt1000_WireDet{sep_str}Source_MATS"]["x"],
                        kind = "linear",fill_value=(np.NaN,np.NaN),
                        bounds_error= False
                )
        return f_int
    r_int = R_interpolate(data)
    r_series = r_int(t_series)
   
    def i_interpolate(data):
        t_arr = (np.array(data["goal_A_Wire1_WireDet_Source_MATS"]["t"]) 
               + float(data["goal_A_Wire1_WireDet_Source_MATS"]["start"]))
        f_int = interp1d(t_arr,
                        data["goal_A_Wire1_WireDet_Source_MATS"]["x"],
                        kind = "linear",fill_value=(np.NaN,np.NaN),
                        bounds_error= False
                )
        return f_int
    i_int = i_interpolate(data)
    i_series = i_int(t_series)
    # Make i_set after the fact, by rounding to nearest in list of setpoints
    i_set_arr = np.array([round(item,5)
                          for item in i_series])
    i_set_str = np.array(["{:.5f}".format(item) for item in i_set_arr])


    dates = mpl.dates.num2date(t_series/86400)
    data_dict = {"dates":np.array(dates),
                "voltage":np.array(v_series)
                }
    data_dict["R_Pt_1000"] = np.array(r_series)
    data_dict["i_series"] = np.array(i_series)
    data_dict["i_set"] = np.array(i_set_arr)
    data_dict["i_set_str"] = np.array(i_set_str)

    save_data(dataname, data_dir, data_dict)
    return data_dict


def v_range_adjust(data_dict, factor):
    data_dict["voltage"] = factor * data_dict["voltage"]

def save_data(dataname, data_dir, data_dict):
    with open(data_dir + dataname + ".pkl", "wb") as f:
            dill.dump(data_dict, f)
    return data_dict

def load_data(dataname, data_dir):
        with open(data_dir + dataname + ".pkl", "rb") as f:
            data_dict = dill.load(f)
        return data_dict

def join_data(dict_1, dict_2):
    dict_joined = {"dates": np.concatenate((dict_1["dates"], dict_2["dates"]))
                 , "voltage":np.concatenate((dict_1["voltage"], 
                   dict_2["voltage"]))}
    return dict_joined

# define temperature conversions:
###########
# Quick and dirty reference caibration from UHV comissioning (OUTDATED)
v_list =  np.array(
        [6.6861, 13.44219, 27.35371,42.18528,58.5537,77.40546,100.0043,
        128.5273,166.4199,2.18E+02,
        281.0631])
T_fudge = 1.25
T_list = np.array(
        [24.55613258,25.59051025,29.3261707,35.43186468,44.68369048,
        58.13071938,77.06963344,104.1148357,143.0549115,197.3745623,
        259.4811837])*T_fudge

P_list =np.array( 
        [0.668202148,2.687827725,10.93932306,25.30978432,46.8288427,
        77.37949821,119.9560579,179.8447807,266.136707,392.0810162,
        561.9376067])


#Fudge factor to covert calibration to base resistances of new wire
# No good data available:
R_fudge = 70
R_list = (70/66.9018101)*np.array(
[66.9018101,
67.22621035,
68.39778356,
70.31264376,
73.21418993,
77.4314305,
83.37102931,
91.85291219,
104.0652506,
121.1009123,
140.5787142,
]
)

# BODGE: Use R as v standin since current is constant at 1mA, 
# otherwise V gets quadratic scaling
v_list = R_list

def v_to_T(x):
        f_int = interp1d(v_list,
                T_list
                ,
                kind = "cubic"
                ,fill_value="extrapolate"
                )
        return f_int(x)

def T_to_v(x):
        f_int = interp1d(T_list,
                        v_list,
                        kind = "cubic"
                        ,fill_value="extrapolate"
                )
        return f_int(x)

# Fudge facotors to zero the axes to show deltas
P_fudge = 77.37949821 - 1.8 - 5.5
v_fudge = 77.40546
def v_to_P(x):
        f_int = interp1d(v_list,
                P_list - P_fudge
                ,
                kind = "cubic",fill_value="extrapolate"
                )
        return f_int(x)

def P_to_v(x):
        f_int = interp1d(P_list - P_fudge,
                        v_list,
                        kind = "cubic",fill_value="extrapolate"
                )
        return f_int(x)

#data_dict to T_diff
def data_to_T_dict(data_dict):
        #T_wire = np.array([v_to_T(v) for v in data_dict["voltage"]])
        #T_Pt_1000 = np.array([R_to_T(R) for R in data_dict["R_Pt_1000"]])
        T_wire = np.array(v_to_T(data_dict["voltage"]*1000))
        T_Pt_1000 = np.array(R_to_T(data_dict["R_Pt_1000"]))

        T_diff = T_wire - T_Pt_1000
        T_dict = {"dates":data_dict["dates"].copy(),
                     "T_wire" : T_wire,
                     "T_Pt_1000" : T_Pt_1000,
                     "T_diff" : T_diff
                        }
        return T_dict

def plot_T_diff(T_dict, plotname, plot_dir,
                start_date, end_date,
                           utc_offset = 1,
                           large_points = False,
                           ):

    start,end = select_date_indices(T_dict,start_date, end_date)
    dates = T_dict["dates"][start:end]
    T_diff = T_dict["T_diff"][start:end]
        

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    if large_points == False:
        ax1.plot(dates,T_diff,".",markersize=1,
                label = f"T_diff")
    elif large_points == True:
        ax1.plot(dates,T_diff,"x",markersize=3,
                 ls = "-", linewidth = 0.5,
                label = f"T_diff")
        
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"$\Delta$ Temperature [K]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_T_diff"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()


def plot_with_axes(data_dict, plotname, plot_dir,start_date, end_date,
                           utc_offset = 1):

    ###########
    # Quick and dirty reference caibration from UHV comissioning (OUTDATED)
    v_list =  np.array(
            [6.6861, 13.44219, 27.35371,42.18528,58.5537,77.40546,100.0043,
            128.5273,166.4199,2.18E+02,
            281.0631])
    T_list = np.array(
            [24.55613258,25.59051025,29.3261707,35.43186468,44.68369048,
            58.13071938,77.06963344,104.1148357,143.0549115,197.3745623,
            259.4811837])

    P_list =np.array( 
            [0.668202148,2.687827725,10.93932306,25.30978432,46.8288427,
            77.37949821,119.9560579,179.8447807,266.136707,392.0810162,
            561.9376067])


    #Fudge factor to covert calibration to base resistances of new wire
    # No good data available:
    R_fudge = 70
    R_list = (70/66.9018101)*np.array(
    [66.9018101,
    67.22621035,
    68.39778356,
    70.31264376,
    73.21418993,
    77.4314305,
    83.37102931,
    91.85291219,
    104.0652506,
    121.1009123,
    140.5787142,
    ]
    )

    # BODGE: Use R as v standin since current is constant at 1mA, 
    # otherwise V gets quadratic scaling
    v_list = R_list

    def v_to_T(x):
            f_int = interp1d(v_list,
                    T_list
                    ,
                    kind = "cubic"
                    ,fill_value="extrapolate"
                    )
            return f_int(x)

    def T_to_v(x):
            f_int = interp1d(T_list,
                            v_list,
                            kind = "cubic"
                            ,fill_value="extrapolate"
                    )
            return f_int(x)

    # Fudge facotors to zero the axes to show deltas
    P_fudge = 77.37949821 - 1.8 - 5.5
    v_fudge = 77.40546
    def v_to_P(x):
            f_int = interp1d(v_list,
                    P_list - P_fudge
                    ,
                    kind = "cubic",fill_value="extrapolate"
                    )
            return f_int(x)

    def P_to_v(x):
            f_int = interp1d(P_list - P_fudge,
                            v_list,
                            kind = "cubic",fill_value="extrapolate"
                    )
            return f_int(x)

    d_wire = 5e-6
    l_illuminated = 1.6e-2
    A_illuminated = d_wire * l_illuminated
    joules_per_electronvolt = 1.60218e-19
    energy_per_atom = (4.75/2) * joules_per_electronvolt

    def flux_to_v(flux):
            joules_per_electronvolt = 1.60218e-19
            energy_per_atom = (4.75/2) * joules_per_electronvolt
            return P_to_v(flux * A_illuminated * energy_per_atom * 1e6)

    def v_to_flux(v):
            joules_per_electronvolt = 1.60218e-19
            energy_per_atom = (4.75/2) * joules_per_electronvolt
            return v_to_P(v) /( A_illuminated * energy_per_atom * 1e6)

    ###########
    # select data (like in select_dates)
    if type(data_dict) == type({"dict":[]}):
        start,end = select_date_indices(data_dict,start_date, end_date)
        v_series = data_dict["voltage"][start:end]
        dates = data_dict["dates"][start:end]
        # move to German timezone
        dates = np.array([
                date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                for date in dates])
            
    else:
        Exception("Invalid file provided")

    mavg_len = 30 # In seconds
    w = mavg_len * 1 # ~1 measurements per second
    #t_avg_series = moving_average(t_series, w)
    v_avg_series = moving_average(v_series, w)

    try:
            sigma= np.std(v_avg_series - v_series[w//2:-w//2])
    except:
            sigma= np.std(v_avg_series[:-1] - v_series[w//2:-w//2]) 
    #############
    # make figure
    fig = plt.figure(0, figsize=(8,6.5))
    fig.set_size_inches(15, 6)
    ax1=plt.gca()
    ax1.set_aspect

    mavg_len = 30 # In seconds
    w = mavg_len * 1 # ~1 measurements per second
    #t_avg_series = moving_average(t_series, w)
    v_avg_series = moving_average(v_series, w)

    try:
            sigma= np.std(v_avg_series - v_series[w//2:-w//2])
    except:
            sigma= np.std(v_avg_series[:-1] - v_series[w//2:-w//2]) 

    ax1.plot(dates,v_series*1000,".",markersize=1,
            label = f"data, std vs mavg = {sigma * 1000:.1e}")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"Voltage [mV]")

    secax_1 = ax1.secondary_yaxis(-0.12, functions=(v_to_T, T_to_v))
    secax_1.set_ylabel('T [°C]')

    secax_2 = ax1.secondary_yaxis(-0.24, functions=(v_to_P, P_to_v))
    secax_2.set_ylabel(r'$\Delta$ Power [µW]')

    secax_3 = ax1.secondary_yaxis(-0.35, functions=(v_to_flux, flux_to_v))
    secax_3.set_ylabel(r'$\Delta$ hits [Atoms per s]')


    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_axes"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

#interpolate pressure, does not work. too naive
def p_interpolate(p_dict):
        t_arr = np.array([date.timestamp() for date in p_dict["dates"]])
        f_int = interp1d(t_arr,
                        p_dict["pressure"],
                        kind = "linear",fill_value="extrapolate"
                )
        return f_int

def make_corr_dict(diff_dict, p_dict):
        # sample pressure at voltage/temperature times
        p_int = p_interpolate(p_dict)
        t_arr = np.array([date.timestamp() for date in diff_dict["dates"]])
        p_arr = p_int(t_arr)
        corr_arr = (diff_dict["T_diff"] - (diff_dict["T_diff"][0] +0.09))/p_arr
        return {"dates" : diff_dict["dates"].copy(),
                "corr" : corr_arr}

def plot_T_diff_p_corr(corr_dict, plotname, plot_dir,
                       start_date, end_date,
                           utc_offset = 1,
                           large_points = False,
                           ):

    start,end = select_date_indices(corr_dict,start_date, end_date)
    dates = corr_dict["dates"][start:end]
    corr = corr_dict["corr"][start:end]
        

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    if large_points == False:
        ax1.plot(dates,corr,".",markersize=1,
                #label = f"T_diff"
                )
    elif large_points == True:
        ax1.plot(dates,corr,"x",markersize=3,
                 ls = "-", linewidth = 0.5,
                #label = f"T_diff"
                )
        
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_ylim(-1e4,5e3)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"T/p [K/mbar]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_T_diff_by_p"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

# plot wire Voltage and Pg60 pressur together
def plot_V_and_p(data_dict,p_dict, plotname, plot_dir,
                 start_date, end_date,
                           utc_offset = 1,
                           large_points = False,
                           p_dict_2 = None,
                           ):
                           
    start,end = select_date_indices(data_dict,start_date, end_date)
    v_series = data_dict["voltage"][start:end]
    dates = data_dict["dates"][start:end]
    # move to German timezone
    dates = np.array([
                date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                for date in dates])

    # make p_series
    start,end = select_date_indices(p_dict,start_date, end_date)
    p_series = p_dict["pressure"][start:end]
    dates_p = p_dict["dates"][start:end]
    # move to German timezone
    dates_p = np.array([
                date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                for date in dates_p])

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    if large_points == False:
        ax1.plot(dates,v_series*1000,".",markersize=1,
                 label = "Wire"
                #label = f"data, std vs mavg = {sigma * 1000:.1e}"
                )
    elif large_points == True:
        ax1.plot(dates,v_series*1000,"x",markersize=3,
                 ls = "-", linewidth = 0.5,
                #label = f"data, std vs mavg = {sigma * 1000:.1e}"
                )


    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
    ax1.set_ylabel(r"Voltage [mV]",color="C0")
    #plt.legend(shadow=True,loc = "upper left")

    ax1.grid(True)
    # twin object for two different y-axis on the sample plot
    ax2=ax1.twinx()
    # make a plot with different y-axis using second axis object
    ax2.plot(dates_p, p_series,color="C1",
             label = "PG60")
    ax2.set_ylabel("pressure [mbar]",color="C1")

    # make p_series_2
    if not (p_dict_2 == None):
        p_dict = p_dict_2
        start,end = select_date_indices(p_dict,start_date, end_date)
        p_series = p_dict["pressure"][start:end]
        dates_p = p_dict["dates"][start:end]
        # move to German timezone
        dates_p = np.array([
                        date.astimezone(dt.timezone(dt.timedelta(hours=utc_offset)))
                        for date in dates_p])
        
        ax2.plot(dates_p, p_series,color="C2",
             label = "PG9")
        #ax2.set_ylabel("pressure [mbar]",color="C1")

#     plt.xticks(rotation = 45)

#     ax1.set_xlabel(r"Time" + f" {dates[0].tzname()}" )
#     ax1.set_ylabel(r"Voltage [mV]",color="C0")


#     h, l = ax1.get_legend_handles_labels()
#     h2, l2 = ax2.get_legend_handles_labels()
#     #print("labels:", l)
#     l[1] = r"5$\cdot$error clipped points"
#     select = [0,1,3,4]
#     ax2.legend([h[i] for i in select], [l[i] for i in select],
#                #shadow = True,
#                framealpha = 0.5,
#                loc = "upper left",
#                fontsize = "small",
               
#                ncol = 2
#                )


    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_V_and_p"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    ax2.cla()
    plt.cla()
    plt.clf()


if __name__ =="__main__":
    pass