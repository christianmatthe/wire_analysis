########## imports and setup
#from this import d
#from tracemalloc import start
#from audioop import avg
import re
import numpy as np
import matplotlib.pyplot as plt
import os
import time
from datetime import datetime
from scipy.interpolate import interp1d
import dill
#from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from scipy.optimize import curve_fit
from decimal import Decimal

# import functions from other file
from . import Voltage_base_analysis as vba

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)

plot_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep + "output/"
data_dir = os.path.dirname(os.path.abspath(__file__)) + os.sep + "data/"
#######################

def pt1000_T_to_R(T):
    #in Celcius
    #where

    # R = resistance of sensor in ohms
    # R0 = resistance at zero degrees Celsius in ohms (100 for Pt100, and 1000 for Pt1000)
    # T = temperature in degrees Celsius
    A = 3.9083*10**-3
    B = -5.775*10**-7
    C = -4.183*10**-12
    # C = -4.183*10**-12 for T < 0 degrees Celsius
    # C = 0 for T >= 0 degrees Celsius
    R0 = 1000

    cond = (T < 0)
    R = np.piecewise(T, 
        [cond, ~cond],
        [
        R0 * (1 + A*T[cond] + B*T[cond]**2 -100*C*T[cond]**3 + C*T[cond]**4),
        R0 * (1 + A*T[~cond]  + B*T[~cond] **2) 
        ]
            )
    # if T < 0:
    #     R = R0 * (1 + A*T + B*T^2 -100*C*T^3 + C*T^4)
    # else:
    #     R = R0 * (1 + A*T + B*T^2)
    return R

def pt1000_interpolate():
    Tlst = np.linspace(-50,800,num=851)
    Rlst = pt1000_T_to_R(Tlst)
    #print("R = " , Rlst)
    return (lambda T: interp1d(Tlst,Rlst)(T), lambda R: interp1d(Rlst,Tlst)(R))


def make_index_dict(data_dict):
    """
    create a dictionary that specifies the index range for which a certain
    current was applied
    """
    data = data_dict
    if data_dict["i_set_str"] is not None:
        # for new "SlowDash" data
        unique_i_arr = np.unique(data["i_set_str"])
        #print("unique_i_arr: ", unique_i_arr)

        index_dict={}
        for current in unique_i_arr:
            index_dict[f"{current}"] = []
            # find all datapoints with a certain current
            index_section = np.where(data["i_set_str"] == current)[0]
            if len(index_section) < 8:
                del index_dict[f"{current}"]
                continue
            # cut index section at non consectutive indexes
            # create list of indexes where previous index is not consecutive
            non_consec = [i for i in range(1,len(index_section))
            if (index_section[i]
                != index_section[i-1] + 1)
            ]
            # prepend 0th instance
            non_consec = [0] + non_consec

            # print(non_consec)
            index_dict[f"{current}"] = [
                index_section[non_consec[i]:non_consec[i+1]]
                if i != len(non_consec) - 1
                else index_section[non_consec[i]:]
                for i in range(len(non_consec))
                ]

    else:
        unique_i_arr = np.unique(data["i_set"])
        #print("unique_i_arr: ", unique_i_arr)

        index_dict={}
        for current in unique_i_arr:
            index_dict[f"{current}"] = []
            # find all datapoints with a certain current
            index_section = np.where(data["i_set"] == current)[0]
            if len(index_section) < 5:
                del index_dict[f"{current}"]
                continue
            # cut index section at non consectutive indexes
            # create list of indexes where previous index is not consecutive
            non_consec = [i for i in range(1,len(index_section))
            if (index_section[i]
                != index_section[i-1] + 1)
            ]
            # prepend 0th instance
            non_consec = [0] + non_consec

            # print(non_consec)
            index_dict[f"{current}"] = [
                index_section[non_consec[i]:non_consec[i+1]]
                if i != len(non_consec) - 1
                else index_section[non_consec[i]:]
                for i in range(len(non_consec))
                ]
    return index_dict 

def make_index_dict_sd(data_dict):
    """
    create a dictionary that specifies the index range for which a certain
    current was applied
    """
    data = data_dict
    unique_i_arr = np.unique(data["i_set_str"])
    #print("unique_i_arr: ", unique_i_arr)

    index_dict={}
    for current in unique_i_arr:
        index_dict[f"{current}"] = []
        # find all datapoints with a certain current
        index_section = np.where(data["i_set_str"] == current)[0]
        if len(index_section) < 5:
            del index_dict[f"{current}"]
            continue
        # cut index section at non consectutive indexes
        # create list of indexes where previous index is not consecutive
        non_consec = [i for i in range(1,len(index_section))
        if (index_section[i]
             != index_section[i-1] + 1)
        ]
        # prepend 0th instance
        non_consec = [0] + non_consec

        # print(non_consec)
        index_dict[f"{current}"] = [
            index_section[non_consec[i]:non_consec[i+1]]
            if i != len(non_consec) - 1
            else index_section[non_consec[i]:]
            for i in range(len(non_consec))
            ]
    return index_dict 

def plot_calib(data_dict, plotname,index_arr = None):
    if index_arr is None:
        v_series = data_dict["voltage"]
        dates = data_dict["dates"]
    else:
        v_series = np.array(data_dict["voltage"])[index_arr]
        dates = np.array(data_dict["dates"])[index_arr]
    
    sigma = np.std(v_series[4:])

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.plot(dates,v_series*1000,".",# markersize=1,
            label = f"data, std = {sigma * 1000:.1e}")
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

def plot_calib_sectioned(data_dict, plotname, clip_leading = 4):
    index_dict = make_index_dict(data_dict)

    # Plot each section individually
    for key in index_dict.keys():
        for i,section in enumerate(index_dict[key]):
            v_series = np.array(data_dict["voltage"])[section][clip_leading:]
            dates = np.array(data_dict["dates"])[section][clip_leading:]

            sigma = np.std(v_series)

            fig = plt.figure(0, figsize=(8,6.5))
            ax1=plt.gca()
            ax1.plot(dates,v_series*1000,".",# markersize=1,
                    label = f"data, std = {sigma * 1000:.1e}")
            # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
            #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"Voltage [mV]")

    plt.grid(True)
    # plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_sec"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def plot_calib_section(data_dict, plotname, key_list, clip_leading = 4):
    index_dict = make_index_dict(data_dict)
    # index_dict = {key:index_dict[key] for key in key_list}
    # Plot each section individually
    for key in index_dict.keys():
        for i,section in enumerate(index_dict[key]):
            v_series = np.array(data_dict["voltage"])[section][clip_leading:]
            dates = np.array(data_dict["dates"])[section][clip_leading:]

            sigma = np.std(v_series)

            fig = plt.figure(0, figsize=(8,6.5))
            ax1=plt.gca()
            ax1.plot(dates,v_series*1000,".",# markersize=1,
                    label = f"{key}_{i}, std = {sigma * 1000:.1e}")
            # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
            #         label = f"moving average {mavg_len}s")

    plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Time")
    ax1.set_ylabel(r"Voltage [mV]")

    plt.grid(True)
    #plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "keys" + f"{key_list}"+"_sec"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def plot_calib_kappa(data_dict, plotname,index_arr = None):
    if index_arr is None:
        v_series = data_dict["voltage"]
        dates = data_dict["dates"]
    else:
        v_series = np.array(data_dict["voltage"])[index_arr]
        dates = np.array(data_dict["dates"])[index_arr]
    
    sigma = np.std(v_series[4:])

    v_series_masked = sigma_clip(v_series,masked = True)
    mask = np.logical_not(v_series_masked.mask)

    sigma = np.std(v_series[mask])

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.plot(dates[mask],v_series[mask]*1000,".",# markersize=1,
            label = f"data, std = {sigma * 1000:.1e}")
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


def R_err(i,i_err,v,v_err): 
    res = np.sqrt(
            ((1/i) * v_err )**2
            + ((v/(i**2)) * i_err)**2
                  )
    return res

def P_err(i,i_err,v,v_err): 
    res = np.sqrt(
            (i * v_err )**2
            + (v * i_err)**2
                  )
    return res

def raw_dict_to_avg(data_dict,  key_list = "all",
                    clip_leading = 4, clip_trailing =2, index_dict = "default"):
    if index_dict == "default":
        index_dict = make_index_dict(data_dict)
    if key_list == "all":
        key_list = index_dict.keys()

    #initialize dict 
    avg_dict = {
        "v" : [],
        "i" : [],
        "i_set" : [],
        "v_err" : [],
        "i_err" : [],
        "v_point_err" : [],
        "i_point_err" : [],
        "R_Pt_1000" : [],
        "R_Pt_1000_err" : [],
        "R_Pt_1000_point_err" : [],

    }
    # Add average date?

    # Split series and average each
    for key in key_list:
        for i,section in enumerate(index_dict[key]):
            v_series = np.array(data_dict["voltage"])[section][clip_leading:
                                                        -clip_trailing]
            i_series = np.array(data_dict["i_series"])[section][clip_leading:
                                                        -clip_trailing]
            i_set = np.array(data_dict["i_set"])[section][clip_leading:
                                                        -clip_trailing]
            dates = np.array(data_dict["dates"])[section][clip_leading:
                                                        -clip_trailing]
            R_Pt = np.array(data_dict["R_Pt_1000"])[section][clip_leading:
                                                        -clip_trailing]
            

            # convert unit: we want mA and mV
            v_series = v_series * 1000
            i_series = i_series * 1000

            n_meas = len(v_series)
            v_mean = np.average(v_series)
            i_mean = np.average(i_series)
            R_Pt_mean = np.average(R_Pt)
            i_set_mean = np.average(i_set)

            v_sigma = np.std(v_series)
            i_sigma = np.std(i_series)
            R_Pt_sigma = np.std(R_Pt)
            v_mean_err = v_sigma / np.sqrt(n_meas)
            i_mean_err = i_sigma / np.sqrt(n_meas)
            R_Pt_mean_err = R_Pt_sigma / np.sqrt(n_meas)

            avg_dict["v"].append(v_mean)
            avg_dict["i"].append(i_mean)
            avg_dict["i_set"].append(i_set_mean)
            avg_dict["v_point_err"].append(v_sigma) # stat error of 1 meas.
            avg_dict["i_point_err"].append(i_sigma)
            avg_dict["v_err"].append(v_mean_err) # error of the mean
            avg_dict["i_err"].append(i_mean_err)
            avg_dict["R_Pt_1000"].append(R_Pt_mean)
            avg_dict["R_Pt_1000_point_err"].append(R_Pt_sigma)
            avg_dict["R_Pt_1000_err"].append(R_Pt_mean_err)


    # convert to nummpy  arrays
    for key in avg_dict:
        avg_dict[key] = np.array(avg_dict[key])
    
    #make calculated quantities
    avg_dict["R"] = avg_dict["v"] / avg_dict["i"]
    avg_dict["P"] = avg_dict["v"] * avg_dict["i"]
    avg_dict["R_err"] = R_err(avg_dict["i"], avg_dict["i_err"],
                              avg_dict["v"], avg_dict["v_err"]
                                )
    avg_dict["P_err"] = P_err(avg_dict["i"], avg_dict["i_err"],
                              avg_dict["v"], avg_dict["v_err"]
                                )

    return avg_dict



def plot_R_vs_P(avg_dict, plotname, plot_dir = plot_dir):
    ############## Resistance vs power plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.errorbar(avg_dict["P"],
             avg_dict["R"],
             yerr = avg_dict["R_err"],
             xerr = avg_dict["P_err"],
             fmt = ".",
             # markersize=1,
             label = f"data")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

def plot_R_vs_P(avg_dict, plotname, plot_dir = plot_dir, inset = False):
    ############## Resistance vs power plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    ax1.errorbar(avg_dict["P"],
             avg_dict["R"],
             yerr = avg_dict["R_err"],
             xerr = avg_dict["P_err"],
             fmt = ".",
             #markersize=4,
             label = f"data")
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")
        #Add secondary axis with current
    def i_to_P(x):
            f_int = interp1d(avg_dict["i"][~np.isnan(avg_dict["i"])],
                    avg_dict["P"][~np.isnan(avg_dict["i"])]
                    #  ,
                    #  kind = "cubic"
                    ,fill_value="extrapolate"
                    )
            return f_int(x)

    def P_to_i(x):
            f_int = interp1d(avg_dict["P"][~np.isnan(avg_dict["i"])],
                            avg_dict["i"][~np.isnan(avg_dict["i"])]
                            #,kind = "cubic"
                            ,fill_value="extrapolate"
                            )
            out  = f_int(x.astype(float))
            return out
    secax = ax1.secondary_xaxis("top", functions=(P_to_i, i_to_P))
    secax.set_xlabel('Current [mA]')
    # secax.set_xticks([i/10 for i in range(0,30,5)])
    i_lst = avg_dict["i"][~np.isnan(avg_dict["i"])]
    xticks = [0] + list(np.round(np.linspace(i_lst[0]+i_lst[-1]*0.3 ,i_lst[-1]
                                             , 8)
                        , decimals=1))
    secax.set_xticks(xticks
                     )

    if inset == True:
        # dx = max(avg_dict["P"]) - min(avg_dict["P"]) 
        # x1 = min(avg_dict["P"]) - dx * 0.001
        # x2 = min(avg_dict["P"]) + dx * 0.01
        # dy = max(avg_dict["R"]) - min(avg_dict["R"]) 
        # y1 = min(avg_dict["R"]) - dy * 0.002
        # y2 = min(avg_dict["R"]) + dy * 0.08

        #BAsed on ordered list bottom 10
        dx = avg_dict["P"][9] - avg_dict["P"][0] 
        x1 = avg_dict["P"][0] - dx * 0.05
        x2 = avg_dict["P"][9] + dx * 0.05
        dy = avg_dict["R"][9] - avg_dict["R"][0] 
        y1 = avg_dict["R"][0] - dy * 0.05
        y2 = avg_dict["R"][9] + dy * 0.05

        axins = ax1.inset_axes(
            [0.62, 0.065, 0.36, 0.40],
            xlim=(x1, x2), ylim=(y1, y2), #xticklabels=[], yticklabels=[],
            )
        # axins = ax1.inset_axes(
        #     [1.06, 0.00, 1.07, 1.00],
        #     xlim=(x1, x2), ylim=(y1, y2), #xticklabels=[], yticklabels=[]
        #     )
        axins.errorbar(avg_dict["P"],
             avg_dict["R"],
             yerr = avg_dict["R_err"],
             xerr = avg_dict["P_err"],
             fmt = ".",
             # markersize=1,
             label = f"data")
        axins.grid(True)
        ax1.indicate_inset_zoom(axins, edgecolor="black",label = None)
        plotname = plotname + "_inset"
        secax_ins = axins.secondary_xaxis("top", functions=(P_to_i, i_to_P))

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi, bbox_inches='tight')
    ax1.cla()

def fit_base_R(avg_dict, plotname, plot_dir = plot_dir, inset = False):
    # possibly undercuts 2nd order correction effects (R is left uncorrected)
    # def R_func(i, m, R_0, i_offset):
    #     u = (R_0 * (i + i_offset))/(1 - m * (i + i_offset)**2)
    #     R = m* u*(i + i_offset) + R_0
    #     return R 
    def v_func(i, m, R_0, i_offset):
        v = (R_0 * (i + i_offset))/(1 - m * (i + i_offset)**2)
        return v

    # try a different idea
    # def R_func(i, m, R_0, i_offset):
    #     R = R_0 /(1 - m *(i+i_offset)**2)
    #     return R 

    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]
    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    #fit
    # #HACK:
    # def v_func(i, m, R_0, i_offset):
    #     i_offset = 0.000105
    #     v = (R_0 * (i + i_offset))/(1 - m * (i + i_offset)**2)
    #     return v
    # popt, pcov = curve_fit(v_func, i_arr, v_arr, p0 = [0.135,66.14],
    #                        sigma = v_err, absolute_sigma=True,
    #                        #bounds = ((0.1,60,-0.00035),(10,75, -0.00025))
    #                        bounds = ((0.1,60),(10,75))
    # )
    # #End HACK

    # #HACK test V offset
    # def v_func(i, m, R_0, i_offset, v_offset):
    #     v = (R_0 * (i + i_offset))/(1 - m * (i + i_offset)**2) + v_offset
    #     return v
    # popt, pcov = curve_fit(v_func, i_arr, v_arr, p0 = [0.135,66.14,-0.000105,0],
    #                     sigma = v_err, absolute_sigma=True,
    #                     #bounds = ((0.1,60,-0.00035),(10,75, -0.00025))
    #                 bounds = ((0.1,60,-0.0003, -0.01),(10,75, 0.0003, +0.01))
    # )
    # v_off = popt[3]
    # v_arr_off = v_arr + v_off
    # #HACK
    # v_arr = v_arr_off

    popt, pcov = curve_fit(v_func, i_arr, v_arr, p0 = [0.135,66.14,-0.000105],
                           sigma = v_err, absolute_sigma=True,
                           #bounds = ((0.1,60,-0.00035),(10,75, -0.00025))
                           bounds = ((0.1,60,-0.0005),(10,75, 0.0005))
    )


    offset = popt[2]
    # #HACK
    # offset = -0.000105
    # popt[2] = offset
    ###
    offset_err = pcov[2][2]
    i_arr_off = i_arr + offset
    i_err_off = np.sqrt(i_err **2 + offset_err **2)

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # old no errorbars
    # ax1.plot(i_arr_off * v_arr ,
    #          v_arr / i_arr_off,
    #          ".",# markersize=1,
    #          label = f"data, offset_{1000*offset:3f} [µA]")
    # # with errorbars
    ax1.errorbar(i_arr_off * v_arr ,
             v_arr / i_arr_off,
             yerr = R_err(i_arr_off, i_err_off, v_arr, v_err),
             xerr = P_err(i_arr_off, i_err_off, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = "data, " + r"$I_{\rm off}$=" + f"{1000*offset:.4f} [µA]")

    # Plot Fit
    xdata=p_arr
    plt.plot(xdata, v_func(i_arr, *popt)/i_arr_off, 'r-',
         label=('fit: m=%5.3f [Ohm/µW],' % popt[0]
                +'\n' +  r"    $R_0=$" "%5.3f [Ohm]," % popt[1]
                + '\n' + r"    $I_{\rm off}$=" + f"{1000*popt[2]:5.4f} [µA]"
                ))
    # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    #         label = f"moving average {mavg_len}s")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    ######
    #Add secondary axis with current
    def i_to_P(x):
            f_int = interp1d(avg_dict["i"][~np.isnan(avg_dict["i"])],
                    avg_dict["P"][~np.isnan(avg_dict["i"])]
                    #  ,
                    #  kind = "cubic"
                    ,fill_value="extrapolate"
                    )
            return f_int(x)

    def P_to_i(x):
            f_int = interp1d(avg_dict["P"][~np.isnan(avg_dict["i"])],
                            avg_dict["i"][~np.isnan(avg_dict["i"])]
                            #,kind = "cubic"
                            ,fill_value="extrapolate"
                            )
            out  = f_int(x.astype(float))
            return out
    secax = ax1.secondary_xaxis("top", functions=(P_to_i, i_to_P))
    secax.set_xlabel('Current [mA]')
    # secax.set_xticks([i/10 for i in range(0,30,5)])
    i_lst = avg_dict["i"][~np.isnan(avg_dict["i"])]
    xticks = [0] + list(np.round(np.linspace(i_lst[0]+i_lst[-1]*0.3 ,i_lst[-1]
                                             , 8)
                        , decimals=1))
    secax.set_xticks(xticks
                     )
    #HACK to set xtiks on low current data
    if inset == False:
        secax.set_xticks([0.0,0.03,0.06,0.09]
                     )

    if inset == True:
        # dx = max(avg_dict["P"]) - min(avg_dict["P"]) 
        # x1 = min(avg_dict["P"]) - dx * 0.001
        # x2 = min(avg_dict["P"]) + dx * 0.01
        # dy = max(avg_dict["R"]) - min(avg_dict["R"]) 
        # y1 = min(avg_dict["R"]) - dy * 0.002
        # y2 = min(avg_dict["R"]) + dy * 0.08

        #BAsed on ordered list bottom 10
        xlst = i_arr_off * v_arr
        ylst = v_arr / i_arr_off
        dx = xlst[9] - xlst[0] 
        x1 = xlst[0] - dx * 0.05
        x2 = xlst[9] + dx * 0.05
        dy = ylst[9] - ylst[0] 
        y1 = ylst[0] - dy * 0.05
        y2 = ylst[9] + dy * 0.05

        axins = ax1.inset_axes(
            [0.62, 0.065, 0.36, 0.40],
            xlim=(x1, x2), ylim=(y1, y2), #xticklabels=[], yticklabels=[],
            )
        # axins = ax1.inset_axes(
        #     [1.06, 0.00, 1.07, 1.00],
        #     xlim=(x1, x2), ylim=(y1, y2), #xticklabels=[], yticklabels=[]
        #     )
        axins.errorbar(i_arr_off * v_arr ,
             v_arr / i_arr_off,
             yerr = R_err(i_arr_off, i_err_off, v_arr, v_err),
             xerr = P_err(i_arr_off, i_err_off, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data with offset")
        axins.plot(xdata, v_func(i_arr, *popt)/i_arr_off, 'r-',
        #  label=(('fit: m=%5.3f [Ohm/µW],\n R_0=%5.3f [Ohm],' 
        #           + '\n i_off=%5.6f [mA]')
        #         % tuple(popt))
                )
        axins.grid(True)
        ax1.indicate_inset_zoom(axins, edgecolor="black",label = None)
        plotname = plotname + "_inset"
        secax_ins = axins.secondary_xaxis("top", functions=(P_to_i, i_to_P))


    ######

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    # troubleshoot plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # ax1.plot(i_arr ,
    #          v_arr / i_arr,
    #          ".",# markersize=1,
    #          label = f"data, offset_{1000*offset:3f} [µA]")
    ax1.errorbar(i_arr ,
             v_arr / i_arr,
             yerr = avg_dict["R_err"][~np.isnan(avg_dict["i"])],
             xerr = avg_dict["i_err"][~np.isnan(avg_dict["i"])],
             fmt = ".",# markersize=1,
             #label = f"data, offset_{1000*offset:3f} [µA]"
             label = f"data with offset")

    # Plot Fit
    xdata=i_arr
    plt.plot(xdata, v_func(i_arr, *popt)/i_arr, 'r-',
         label=(('fit: m=%5.3f [Ohm/µW],\n R_0=%5.3f [Ohm],' 
                  + '\n i_off=%5.6f [mA]')
                % tuple(popt[0:3]))
                )
    ax1.set_xlabel(r"current [mA]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + "trouble"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    return popt, pcov

def basic_R_over_P_calib(avg_dict, plotname, plot_dir = plot_dir):
    # fit 0.9mA to 1.1mA range directly with a
    def fit_func(P, m, R_0):
        R = R_0 + m * P
        return R


    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]
    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)
    # #fit
    popt, pcov = curve_fit(fit_func, p_arr, r_arr, p0 = [0.135,66.7],
                           sigma = R_err(i_arr, i_err, v_arr, v_err),
                           absolute_sigma=True,
                           bounds = ((0.05,60),(1.0,75))
    )

    #########Plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    plt.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW],
         R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
             popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
                  ))
                )

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    #Add secondary axis with current
    def i_to_P(x):
            f_int = interp1d(avg_dict["i"][~np.isnan(avg_dict["i"])],
                    avg_dict["P"][~np.isnan(avg_dict["i"])]
                    #  ,
                    #  kind = "cubic"
                    ,fill_value="extrapolate"
                    )
            return f_int(x)

    def P_to_i(x):
            f_int = interp1d(avg_dict["P"][~np.isnan(avg_dict["i"])],
                            avg_dict["i"][~np.isnan(avg_dict["i"])]
                            #,kind = "cubic"
                            ,fill_value="extrapolate"
                            )
            out  = f_int(x.astype(float))
            return out
    secax = ax1.secondary_xaxis("top", functions=(P_to_i, i_to_P))
    secax.set_xlabel('Current [mA]')
    # secax.set_xticks([i/10 for i in range(0,30,5)])
    i_lst = avg_dict["i"][~np.isnan(avg_dict["i"])]
    xticks = [0] + list(np.round(np.linspace(i_lst[0]+i_lst[-1]*0.3 ,i_lst[-1]
                                             , 8)
                        , decimals=1))
    secax.set_xticks(xticks
                     )



    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #################### plot with residuals
    residuals = r_arr - fit_func(p_arr, *popt)
    chi2 = np.sum((residuals/r_err)**2)
    dof = len(residuals) - len(popt)
    chi2_red = chi2/dof

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    gs.update(#wspace=0.05
            hspace = 0.005
        )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ### ax1
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = r"data, $\chi^2_{red}$"+" = {:2.2f} ".format(chi2_red))

    # Plot Fit
    xdata=p_arr
    ax1.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [$\Omega$/µW],
         R$_0$={:5.3f} $\pm$ {:2.1e} [$\Omega$]'''.format(
             popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
                  ))
                )

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    secax = ax1.secondary_xaxis("top", functions=(P_to_i, i_to_P))
    secax.set_xlabel('Current [mA]')
    # secax.set_xticks([i/10 for i in range(0,30,5)])
    i_lst = avg_dict["i"][~np.isnan(avg_dict["i"])]
    xticks = [0] + list(np.round(np.linspace(i_lst[0]+i_lst[-1]*0.3 ,i_lst[-1]
                                             , 8)
                        , decimals=1))
    secax.set_xticks(xticks
                     )

    ax1.grid(True)
    ax1.legend(shadow=True, fontsize = 13)
    # ax1.tight_layout()

    #ax 2

    ax2.errorbar(p_arr ,
             residuals,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    ax2.plot(xdata, 0.0 * fit_func(p_arr, *popt), 'r-',
        #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW],
        #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
        #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
        #           ))
                )
    ax2.set_xlabel(r"Power [µW]")
    ax2.set_ylabel(r"Residuals [$\Omega$]")

    ax2.grid(True)

    #make custom pruning of uppper tick (do not plot ticks in upper 10%)
    #so that ax2 tick does nto interfere with  ax1 tick
    ax2.locator_params(axis="y", min_n_ticks = 3
                       )
    y_loc = ax2.yaxis.get_majorticklocs()
    #print("y_loc: ", y_loc)
    #print("y_loc[1:-2]: ", y_loc[1:-2])
    #print("ylim: ", ax2.get_ylim())
    y2_min, y2_max = ax2.get_ylim()
    y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
    #print("y_loc: ", y_loc)
    ax2.set_yticks(y_loc)

    fig.tight_layout()
    fig.subplots_adjust(left=0.2)

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_with_residuals"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    fig.clf()
    plt.close()

    ######## Make Pt_1000 Plot REALLY SHOULD BE ITS OWN FUNCTIONIN A CLASS
    R_Pt = avg_dict["R_Pt_1000"][~np.isnan(avg_dict["i"])]
    R_Pt_err = avg_dict["R_Pt_1000_err"][~np.isnan(avg_dict["i"])]

    r = np.corrcoef(residuals, R_Pt)
    print("Correlation matrix:", r)

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    print(R_Pt_err)
    print(R_Pt)
    ax1.errorbar(p_arr ,
             R_Pt,
             yerr = R_Pt_err,
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data" + " , correlation {:.3f}".format(r[0,1]))

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance Pt1000 [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname+ "Pt_1000"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    ##### Calculate Correlation: Should also be separate function

    return popt, pcov

def plot_poly_k(avg_dict, plotname):
    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]

    #HACK remove duplicates should they exist:
    i_arr = np.unique(i_arr)
    v_arr = np.unique(v_arr)
    v_err = np.unique(v_err)
    i_err = np.unique(i_err)


    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)
    #########
    # Fit original data with a 4th order polynomial, and then derive.
    # fit_func = lambda x,c0,c1,c2,c3,c4: (
    #     c4*x**4 + c3*x**3 + c2*x**2 + c1*x**1 + c0)

    fit_func = lambda x,c4,c3,c2,c1,c0: np.poly1d([c4, c3, c2, c1, c0])(x)
    popt, pcov = curve_fit(fit_func, p_arr, r_arr,  p0 = [0,0,0,0.135,66.7],
                           sigma = R_err(i_arr, i_err, v_arr, v_err),
                        #    absolute_sigma=True,
                        #    bounds = ((0.05,60,-1),(1.0,75,10))
    )
    #########Plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    plt.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label="fit"
                )

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #################### plot with residuals
    residuals = r_arr - fit_func(p_arr, *popt)
    chi2 = np.sum((residuals/r_err)**2)
    dof = len(residuals) - len(popt)
    chi2_red = chi2/dof

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    gs.update(#wspace=0.05
            hspace = 0.005
        )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ### ax1
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = r"data, $\chi^2_{red}$"+" = {:2.2f} ".format(chi2_red))

    # Plot Fit
    xdata=p_arr
    ax1.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label="fit"
                )

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    ax1.grid(True)
    ax1.legend(shadow=True, fontsize = 13)
    # ax1.tight_layout()

    #ax 2

    ax2.errorbar(p_arr ,
             residuals,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    ax2.plot(xdata, 0.0 * fit_func(p_arr, *popt), 'r-',
        #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW],
        #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
        #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
        #           ))
                )
    ax2.set_xlabel(r"Power [µW]")
    ax2.set_ylabel(r"Residuals [$\Omega$]")

    ax2.grid(True)

    #make custom pruning of uppper tick (do not plot ticks in upper 10%)
    #so that ax2 tick does nto interfere with  ax1 tick
    ax2.locator_params(axis="y", min_n_ticks = 3
                       )
    y_loc = ax2.yaxis.get_majorticklocs()
    #print("y_loc: ", y_loc)
    #print("y_loc[1:-2]: ", y_loc[1:-2])
    #print("ylim: ", ax2.get_ylim())
    y2_min, y2_max = ax2.get_ylim()
    y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
    #print("y_loc: ", y_loc)
    ax2.set_yticks(y_loc)

    fig.tight_layout()
    fig.subplots_adjust(left=0.2)

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_with_residuals"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    fig.clf()
    plt.close()

    ##########################
    def interval_m(i,j):
        m = (r_arr[j]-r_arr[i])/(p_arr[j]-p_arr[i])
        return m
    def interval_k(i,j):
        return 1/interval_m(i,j)

    def interval_k_err(i,j): 
        # TODO Check this
        delta_p = p_arr[j]-p_arr[i]
        delta_r = r_arr[j]-r_arr[i]
        s_p = np.sqrt(p_err[i]**2 + p_err[j]**2)
        s_r = np.sqrt(r_err[i]**2 + r_err[j]**2)
        k_err = np.sqrt(
            ((1 / delta_r) * s_p)**2
            + ((delta_p / delta_r ** 2) * s_r)**2
        )
        return k_err
    
    k_arr = np.array([interval_k(i,i+1) for i in range(len(p_arr)-1)])
    k_errs = np.array([interval_k_err(i,i+1) for i in range(len(p_arr)-1)])

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(r_arr[:-1] ,
             k_arr,
             yerr = k_errs,
             xerr = R_err(i_arr, i_err, v_arr, v_err)[:-1],
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Resistance [$\Omega$]")
    ax1.set_ylabel(r"k [µW/$\Omega$]")


    # Plot derivative based k
    xdata=r_arr
    plt.plot(xdata, 1/(np.poly1d([*popt]).deriv()(p_arr)), 'r-',
         label="k from deriv"
                )
    # Aditionally plot k derived from  fit to power over R rather than other


    
    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "deriv_k_over_R"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    return popt,pcov

def plot_poly_k_P_over_R(avg_dict, plotname):
    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]

    #HACK remove duplicates should they exist:
    i_arr = np.unique(i_arr)
    v_arr = np.unique(v_arr)
    v_err = np.unique(v_err)
    i_err = np.unique(i_err)


    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)
    #########
    # Fit original data with a 4th order polynomial, and then derive.
    # fit_func = lambda x,c0,c1,c2,c3,c4: (
    #     c4*x**4 + c3*x**3 + c2*x**2 + c1*x**1 + c0)

    fit_func = lambda x,c4,c3,c2,c1,c0: np.poly1d([c4, c3, c2, c1, c0])(x)
    popt, pcov = curve_fit(fit_func, r_arr, p_arr,  #p0 = [0,0,0,0.135,66.7],
                           sigma = P_err(i_arr, i_err, v_arr, v_err),
                        #    absolute_sigma=True,
                        #    bounds = ((0.05,60,-1),(1.0,75,10))
    )
    #########Plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(r_arr ,
             p_arr,
             yerr = P_err(i_arr, i_err, v_arr, v_err),
             xerr = R_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    plt.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label="fit"
                )

    #plt.xticks(rotation = 45)

    ax1.set_ylabel(r"Power [µW]")
    ax1.set_xlabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #################### plot with residuals
    residuals = p_arr - fit_func(r_arr, *popt)
    chi2 = np.sum((residuals/p_err)**2)
    dof = len(residuals) - len(popt)
    chi2_red = chi2/dof

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    gs.update(#wspace=0.05
            hspace = 0.005
        )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ### ax1
    ax1.errorbar(r_arr ,
             p_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = r"data, $\chi^2_{red}$"+" = {:2.2f} ".format(chi2_red))

    # Plot Fit
    xdata=r_arr
    ax1.plot(xdata, fit_func(r_arr, *popt), 'r-',
         label="fit"
                )

    #plt.xticks(rotation = 45)

    ax1.set_ylabel(r"Power [µW]")
    ax1.set_xlabel(r"Resistance [$\Omega$]")

    ax1.grid(True)
    ax1.legend(shadow=True, fontsize = 13)
    # ax1.tight_layout()

    #ax 2

    ax2.errorbar(r_arr ,
             residuals,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=r_arr
    ax2.plot(xdata, 0.0 * fit_func(r_arr, *popt), 'r-',
        #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW],
        #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
        #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
        #           ))
                )
    ax2.set_xlabel(r"Resistance [$\Omega$]")
    ax2.set_ylabel(r"Residuals [µW]")

    ax2.grid(True)

    #make custom pruning of uppper tick (do not plot ticks in upper 10%)
    #so that ax2 tick does nto interfere with  ax1 tick
    ax2.locator_params(axis="y", min_n_ticks = 3
                       )
    y_loc = ax2.yaxis.get_majorticklocs()
    #print("y_loc: ", y_loc)
    #print("y_loc[1:-2]: ", y_loc[1:-2])
    #print("ylim: ", ax2.get_ylim())
    y2_min, y2_max = ax2.get_ylim()
    y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
    #print("y_loc: ", y_loc)
    ax2.set_yticks(y_loc)

    fig.tight_layout()
    fig.subplots_adjust(left=0.2)

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_with_residuals"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    fig.clf()
    plt.close()

    ##########################
    def interval_m(i,j):
        m = (r_arr[j]-r_arr[i])/(p_arr[j]-p_arr[i])
        return m
    def interval_k(i,j):
        return 1/interval_m(i,j)

    def interval_k_err(i,j): 
        # TODO Check this
        delta_p = p_arr[j]-p_arr[i]
        delta_r = r_arr[j]-r_arr[i]
        s_p = np.sqrt(p_err[i]**2 + p_err[j]**2)
        s_r = np.sqrt(r_err[i]**2 + r_err[j]**2)
        k_err = np.sqrt(
            ((1 / delta_r) * s_p)**2
            + ((delta_p / delta_r ** 2) * s_r)**2
        )
        return k_err
    
    k_arr = np.array([interval_k(i,i+1) for i in range(len(p_arr)-1)])
    k_errs = np.array([interval_k_err(i,i+1) for i in range(len(p_arr)-1)])

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(r_arr[:-1] ,
             k_arr,
             yerr = k_errs,
             xerr = R_err(i_arr, i_err, v_arr, v_err)[:-1],
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Resistance [$\Omega$]")
    ax1.set_ylabel(r"k [µW/$\Omega$]")


    # Plot derivative based k
    xdata=r_arr
    plt.plot(xdata, (np.poly1d([*popt]).deriv()(xdata)), 'r-',
         label="k from deriv"
                )
    # Aditionally plot k derived from  fit to power over R rather than other


    
    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "deriv_k_over_R"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    # Aditionally plot k derived from  fit to power over R rather than other
    # way around


    return popt,pcov


def plot_interval_k(avg_dict, plotname):

    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]

    #HACK remove duplicates should they exist:
    i_arr = np.unique(i_arr)
    v_arr = np.unique(v_arr)
    v_err = np.unique(v_err)
    i_err = np.unique(i_err)


    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)

    def interval_m(i,j):
        m = (r_arr[j]-r_arr[i])/(p_arr[j]-p_arr[i])
        return m
    def interval_k(i,j):
        return 1/interval_m(i,j)

    def interval_k_err(i,j): 
        # TODO Check this
        delta_p = p_arr[j]-p_arr[i]
        delta_r = r_arr[j]-r_arr[i]
        s_p = np.sqrt(p_err[i]**2 + p_err[j]**2)
        s_r = np.sqrt(r_err[i]**2 + r_err[j]**2)
        k_err = np.sqrt(
            ((1 / delta_r) * s_p)**2
            + ((delta_p / delta_r ** 2) * s_r)**2
        )
        return k_err

    

    # #########
    # def fit_func(P, m, R_0):
    #     R = R_0 + m * P
    #     return R
    m_arr = np.array([interval_m(i,i+1) for i in range(len(p_arr)-1)])
    k_arr = np.array([interval_k(i,i+1) for i in range(len(p_arr)-1)])
    k_errs = np.array([interval_k_err(i,i+1) for i in range(len(p_arr)-1)])

    #########Plot k
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(p_arr[:-1] ,
             k_arr,
             yerr = k_errs,
             xerr = P_err(i_arr, i_err, v_arr, v_err)[:-1],
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"k [µW/$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #########Plot k over R
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(r_arr[:-1] ,
             k_arr,
             yerr = k_errs,
             xerr = R_err(i_arr, i_err, v_arr, v_err)[:-1],
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Resistance [$\Omega$]")
    ax1.set_ylabel(r"k [µW/$\Omega$]")
    # Plot mmovinng average:
    # def moving_average(a, n=7):
    #     ret = np.cumsum(a, dtype=float)
    #     ret[n:] = ret[n:] - ret[:-n]
    #     return ret[n - 1:] / n

    def weighted_moving_average(a,errs, n=4):
        # double sided n: the point itself and n forward and back
        # Weight with 1 / variance
        ret = [np.average(a[i-n : i+n+1], weights = 1/(errs[i-n : i+n+1])**2) 
               for i in range(n,len(a)-n)]
        return np.array(ret)
    #overwrite function
    moving_average = weighted_moving_average
    k_avg_arr = moving_average(k_arr, errs = k_errs)
    r_avg_arr = moving_average(r_arr[:-1] , errs= r_err[:-1])
    ax1.plot(r_avg_arr, k_avg_arr, 'r-',
         label="weighted moving average"
                )

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_over_R"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #########Plot k over R with B-spline
    # NOTE also consider "UnivariateSpline" for the same job.
    def B_spline_fit(x, y, degree = 3):
        k = degree
        from scipy.interpolate import make_lsq_spline, BSpline
        # 7 - 2 = 5 uniform knots inside data range with 2 endpoint knnots
        t  = np.linspace(x[0], x[-1], num=7, endpoint=True)
        #t = np.arange[-1, 0, 1]
        # concatenate with internal knots
        # add k boundary knots 
        # to get k + 1 boundary knots on both sides to make "regular"
        # following: 
        # https://scipy.github.io/devdocs/reference/generated/
        # scipy.interpolate.make_lsq_spline.html
        t = np.r_[(x[0],)*(k),
                t,
                (x[-1],)*(k)]
        spl = make_lsq_spline(x, y, t, k) 
        return spl

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(r_arr[:-1] ,
             k_arr,
             yerr = k_errs,
             xerr = R_err(i_arr, i_err, v_arr, v_err)[:-1],
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Resistance [$\Omega$]")
    ax1.set_ylabel(r"k [µW/$\Omega$]")

    moving_average = weighted_moving_average
    k_avg_arr = moving_average(k_arr, errs = k_errs)
    r_avg_arr = moving_average(r_arr[:-1] , errs= r_err[:-1])
    ax1.plot(r_avg_arr, k_avg_arr, 'r-',
         label="weighted moving average"
                )
    
    xs = np.linspace(r_arr[:-1][0], r_arr[:-1][-1], num=50, endpoint=True)
    # print("r_arr[:-1]",r_arr[:-1])

    #r_arr is not sorted (likely duee to noise) but it has to be 
    r_sorted = np.array(sorted(r_arr[:-1]))
  
    k_sorted = k_arr[[np.argwhere(r_arr == r)[0][0] for r in r_sorted]]
    #print("r_sorted", r_sorted)


    xs = np.linspace(r_sorted[0], r_sorted[-1], num=50, endpoint=True)
    spl = B_spline_fit(r_sorted, k_sorted)
    plt.plot(xs, spl(xs), 'g-', lw=3, label='LSQ B-spline')
    #TODO add weights to spline fit

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_over_R_B-spline"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    
    #########Plot m
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(p_arr[:-1] ,
             m_arr,
            #  yerr = R_err(i_arr, i_err, v_arr, v_err),
            #  xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"m [$\Omega$/µW]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_m"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    # return B_spline
    return spl

    

def P14_R_over_P_calib(avg_dict, plotname, plot_dir = plot_dir):
    # fit 0.9mA to 1.1mA range directly

    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]
    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)
    # #fit
    #   fit_func
    def fit_func(P, m, R_0, a):
        R = R_0 + m * P + a* P**(1/4)
        return R

    popt, pcov = curve_fit(fit_func, p_arr, r_arr, p0 = [0.135,66.7, 0.0],
                           sigma = R_err(i_arr, i_err, v_arr, v_err),
                           absolute_sigma=True,
                           bounds = ((0.05,60,-1),(1.0,75,10))
    )


    #########Plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    label=(("fit: m={:5.4f} $\pm$ {:2.1e} [$\Omega$/µW],".format(
             popt[0], np.sqrt(pcov[0][0]))
          +" \n R_0={:5.3f} $\pm$ {:2.1e} [$\Omega$/µW]".format(
              popt[1], np.sqrt(pcov[1][1])) 
          +" \n a={:5.3f} $\pm$ {:2.1e} [$\Omega$/(µW)^(1/4)]".format(
              popt[2], np.sqrt(pcov[2][2])) 
                  ))
    plt.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label=label
                )

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    #################### plot with residuals
    residuals = r_arr - fit_func(p_arr, *popt)
    chi2 = np.sum((residuals/r_err)**2)
    dof = len(residuals) - len(popt)
    chi2_red = chi2/dof

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
    gs.update(#wspace=0.05
            hspace = 0.005
        )

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ### ax1
    ax1.errorbar(p_arr ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = r"data, $\chi^2_{red}$"+" = {:2.2f} ".format(chi2_red))

    # Plot Fit
    xdata=p_arr
    label=(("fit: m={:5.4f} $\pm$ {:2.1e} [$\Omega$/µW],".format(
             popt[0], np.sqrt(pcov[0][0]))
          +" \n R_0={:5.3f} $\pm$ {:2.1e} [$\Omega$/µW]".format(
              popt[1], np.sqrt(pcov[1][1])) 
          +" \n a={:5.3f} $\pm$ {:2.1e} [$\Omega$/(µW)^(1/4)]".format(
              popt[2], np.sqrt(pcov[2][2])) 
                  ))
    ax1.plot(xdata, fit_func(p_arr, *popt), 'r-',
         label=label)

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    ax1.grid(True)
    ax1.legend(shadow=True, fontsize = 13)
    # ax1.tight_layout()

    #ax 2

    ax2.errorbar(p_arr ,
             residuals,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")

    # Plot Fit
    xdata=p_arr
    ax2.plot(xdata, 0.0 * fit_func(p_arr, *popt), 'r-',
        #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW],
        #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
        #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
        #           ))
                )
    ax2.set_xlabel(r"Power [µW]")
    ax2.set_ylabel(r"Residuals [$\Omega$]")

    ax2.grid(True)

    #make custom pruning of uppper tick (do not plot ticks in upper 10%)
    #so that ax2 tick does nto interfere with  ax1 tick
    ax2.locator_params(axis="y", min_n_ticks = 3
                       )
    y_loc = ax2.yaxis.get_majorticklocs()
    #print("y_loc: ", y_loc)
    #print("y_loc[1:-2]: ", y_loc[1:-2])
    #print("ylim: ", ax2.get_ylim())
    y2_min, y2_max = ax2.get_ylim()
    y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max - y2_min)*0.1]
    #print("y_loc: ", y_loc)
    ax2.set_yticks(y_loc)

    fig.tight_layout()
    fig.subplots_adjust(left=0.2)

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_with_residuals"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()
    fig.clf()
    plt.close()

    ######## Make Pt_1000 Plot REALLY SHOULD BE ITS OWN FUNCTIONIN A CLASS
    R_Pt = avg_dict["R_Pt_1000"][~np.isnan(avg_dict["i"])]
    R_Pt_err = avg_dict["R_Pt_1000_err"][~np.isnan(avg_dict["i"])]

    r = np.corrcoef(residuals, R_Pt)
    print("Correlation matrix:", r)

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    print(R_Pt_err)
    print(R_Pt)
    ax1.errorbar(p_arr ,
             R_Pt,
             yerr = R_Pt_err,
             xerr = P_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data" + " , correlation {:.3f}".format(r[0,1]))

    #plt.xticks(rotation = 45)

    ax1.set_xlabel(r"Power [µW]")
    ax1.set_ylabel(r"Resistance Pt1000 [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname+ "Pt_1000"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    ##### Calculate Correlation: Should also be separate function

    return popt, pcov

def plot_avg_over_index(avg_dict, plotname):

    i_arr = avg_dict["i"][~np.isnan(avg_dict["i"])]
    #print("np.isnan(avg_dict['i'])",np.isnan(avg_dict["i"]))
    v_arr = avg_dict["v"][~np.isnan(avg_dict["i"])]
    v_err = avg_dict["v_err"][~np.isnan(avg_dict["i"])]
    i_err = avg_dict["i_err"][~np.isnan(avg_dict["i"])]
    p_arr = i_arr * v_arr
    r_arr = v_arr / i_arr
    r_err = R_err(i_arr, i_err, v_arr, v_err)
    p_err = P_err(i_arr, i_err, v_arr, v_err)
    #fit
    #########Plot
    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar([i for i in range(len(r_arr))] ,
             r_arr,
             yerr = R_err(i_arr, i_err, v_arr, v_err),
             fmt = ".",
             # markersize=1,
             label = f"data")


    ax1.set_xlabel(r"Index")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    ######### Plot just 0.9 to 1mA corresponding
    r_arr_cut = r_arr[18:36]

    fig = plt.figure(0, figsize=(8,6.5))
    ax1=plt.gca()
    # # with errorbars
    ax1.errorbar([i for i in range(len(r_arr_cut))] ,
             r_arr_cut,
             yerr = R_err(i_arr, i_err, v_arr, v_err)[18:36],
             fmt = ".",
             # markersize=1,
             label = f"data")


    ax1.set_xlabel(r"Index")
    ax1.set_ylabel(r"Resistance [$\Omega$]")

    plt.grid(True)
    plt.legend(shadow=True)
    plt.tight_layout()
    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(plot_dir + plotname + "_cut"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
    ax1.cla()

    return

if __name__ =="__main__": 
################## 2022-12-08 New Automated Calib run 2022-12-09_calib_7
    run_name = "2023-01-09_calib_long"
    vba.prep_data_calib("../SC_downloads/Wire/2023-01-09_calib_long.json"
                     , run_name
                     )
    data_dict = vba.load_data(run_name)

    #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
    out_dir = os.sep + run_name + os.sep
    os.makedirs(plot_dir + out_dir, exist_ok=True)

    ##### See if old scripts will just work
    avg_dict = raw_dict_to_avg(data_dict)
    avg_dict_low = raw_dict_to_avg(data_dict,
            key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
                        #+  [str(0.1*i) for i in range(2,3)]
                        )

    plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
    plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

    fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
    #calib around 1mA
    avg_dict_09to11 = raw_dict_to_avg(data_dict,
                    key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
                                for i in range(1,10)] 
                             + ["0.00" + "{:0>3d}".format(100 + 1*i) 
                                for i in range(1,10)] 
                             )
    avg_dict_09to116 = raw_dict_to_avg(data_dict,
                    key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
                                for i in range(1,10)] 
                             + ["0.00" + "{:0>3d}".format(100 + 1*i) 
                                for i in range(1,10)] 
                             + ["0.00" + "{:0>3d}".format(100 + 1*i) 
                                for i in range(1,61)] 
                             )

    avg_dict_14to16 = raw_dict_to_avg(data_dict,
                    key_list =["0.00" + "{:0>3d}".format(100 + 1*i) 
                                for i in range(40,61)] 
                             )
    avg_dict_ref = raw_dict_to_avg(data_dict,
                    key_list = ["0.00" + "{:0>3d}".format(100) 
                                for i in range(1)] 
                             )

    # plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
    # basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
    # basic_R_over_P_calib(avg_dict_09to116, out_dir +"R_vs_P_fit_09to116")
    # basic_R_over_P_calib(avg_dict_14to16, out_dir +"R_vs_P_fit_14to16")
    # out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
    # os.makedirs(plot_dir + out_dir_2, exist_ok=True)
    # P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")
    # P14_R_over_P_calib(avg_dict_09to116, out_dir_2 +"R_vs_P_fit_09to116")
    # P14_R_over_P_calib(avg_dict_14to16, out_dir_2 +"R_vs_P_fit_14to16")

    # plot_calib_sectioned(data_dict, out_dir + "plot_calib")

    #####
    #Plot interval slopes
    #plot_interval_k(avg_dict_09to116, out_dir +"interval_k_09to116")
    plot_poly_k(avg_dict_09to116, out_dir +"poly_fit_09to116")
    plot_poly_k_P_over_R(avg_dict_09to116, out_dir +"PoR_poly_fit_09to116")

# ################## 2022-12-08 New Automated Calib run 2022-12-09_calib_7
#     run_name = "2022-12-09_calib_7"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-09_calib_7.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_09to116 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 10*i) 
#                                 for i in range(1,7)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     basic_R_over_P_calib(avg_dict_09to116, out_dir +"R_vs_P_fit_09to116")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")
#     P14_R_over_P_calib(avg_dict_09to116, out_dir_2 +"R_vs_P_fit_09to116")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")

# ################## 2022-12-08 New Automated Calib run 2022-12-08_calib_6
#     run_name = "2022-12-08_calib_6"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-08_calib_6.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")

#     ################## 2022-12-08 New Automated Calib run 2022-12-08_calib_5
#     run_name = "2022-12-08_calib_5"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-08_calib_5.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")

# ################## 2022-12-08 New Automated Calib run 2022-12-08_calib_4
#     run_name = "2022-12-08_calib_4"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-08_calib_4.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")


# ################## 2022-12-08 New Automated Calib run 3
#     run_name = "2022-12-08_calib_3"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-08_calib_3.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")


#     ################## 2022-12-05 New Automated Calib run 2
#     run_name = "2022-12-05_calib_2"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-05_calib_2.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     #print('data_dict["R_Pt_1000"][10:-10]',data_dict["R_Pt_1000"][10:-10])
#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#             key_list = ["0.00" + "{:0>3d}".format(1*i) for i in range(1,10)] 
#                         #+  [str(0.1*i) for i in range(2,3)]
#                         )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )
#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")


#     ############## 2022-12-05 New Automated Calib
#     run_name = "2022-12-05_calib_1"
#     vba.prep_data_calib("../SC_downloads/Wire/2022-12-05_calib_1.json"
#                      , run_name
#                      )
#     data_dict = vba.load_data(run_name)

#     out_dir = os.sep + run_name + os.sep
#     os.makedirs(plot_dir + out_dir, exist_ok=True)

#     ##### See if old scripts will just work
#     avg_dict = raw_dict_to_avg(data_dict)
#     avg_dict_low = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(1*i) 
#                                 for i in range(1,10)] 
#                              #+  [str(0.1*i) for i in range(2,3)]
#                              )
#     avg_dict_to_99 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(10*i) 
#                                 for i in range(1,10)]
#                              + ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 


#                              #+  [str(0.1*i) for i in range(2,3)]
#                              )

#     #calib around 1mA
#     avg_dict_09to11 = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(90 + 1*i) 
#                                 for i in range(1,10)] 
#                              + ["0.00" + "{:0>3d}".format(100 + 1*i) 
#                                 for i in range(1,10)] 
#                              )

#     plot_R_vs_P(avg_dict, plotname = out_dir + "R_vs_P")
#     plot_R_vs_P(avg_dict_low, plotname = out_dir +"R_vs_P_to_0.1")
#     plot_R_vs_P(avg_dict_09to11, plotname = out_dir +"R_vs_P_09to11")

#     fit_base_R(avg_dict_low, out_dir +"R_vs_P_offset_fit")
#     fit_base_R(avg_dict_to_99, out_dir +"R_vs_P_offset_fit_to_99")
#     fit_base_R(avg_dict_09to11, out_dir +"R_vs_P_offset_fit_09to11")
#     basic_R_over_P_calib(avg_dict_09to11, out_dir +"R_vs_P_fit_09to11")
#     out_dir_2 = out_dir + os.sep + "4_root_P" + os.sep
#     os.makedirs(plot_dir + out_dir_2, exist_ok=True)
#     P14_R_over_P_calib(avg_dict_09to11, out_dir_2 +"R_vs_P_fit_09to11")

#     plot_calib_sectioned(data_dict, out_dir + "plot_calib")

#     avg_dict_ref = raw_dict_to_avg(data_dict,
#                     key_list = ["0.00" + "{:0>3d}".format(100) 
#                                 for i in range(1)] 
#                              )

#     plot_avg_over_index(avg_dict_ref, out_dir +"R_over index_fit_i_ref")

    ###################
    # vba.prep_data_with_i("../logs/calib_0.01-2.5_high_range_2.txt"
    #                  , "calib_0.01-2.5_high_range_2"
    #                  )
    # data_dict = vba.load_data("calib_0.01-2.5_high_range_2" + "_with_i")
    # #print(data_dict["i_set"])
    # # index_dict = make_index_dict(data_dict)
    # # print(index_dict)

    # # plot_calib_sectioned(data_dict, plotname = f"calib/0.01-2.5"
    # #             )
    # # plot_calib_section(data_dict, plotname = f"calib/0.01-2.5", 
    # #                    key_list = ["1.0"]
    # #             )
    # # plot_calib_section(data_dict, plotname = f"calib/0.01-2.5", 
    # #                    key_list = ["1.1"]
    # #             )

    # # # Plot all individually
    # # index_dict = make_index_dict(data_dict)
    # # for key in index_dict.keys():
    # #     plot_calib_section(data_dict, plotname = f"calib/0.01-2.5", 
    # #                     key_list = [key]
    # #             )

    # # # Plot just low currents
    # # plot_calib_section(data_dict, plotname = f"calib/0.01-0.1", 
    # #                     key_list = [str(0.01*i) for i in range(1,10)]
    # #             )

    # avg_dict = raw_dict_to_avg(data_dict)
    # #print(avg_dict)
    # avg_dict_low = raw_dict_to_avg(data_dict,
    #                 key_list = [str(0.01*i) for i in range(1,10)] 
    #                          #+  [str(0.1*i) for i in range(2,3)]
    #                          )

    # plot_R_vs_P(avg_dict, plotname = "11_Mar_R_vs_P")
    # plot_R_vs_P(avg_dict_low, plotname = "11_Mar_R_vs_P_to_0.1")

    # fit_base_R(avg_dict_low, "R_vs_P_offset_fit")

    ############## Resistance vs power plot
    # offset = -0.01* (65.75 +( 0.46)- 64.25)/64.25
    # i_arr = avg_dict_low["i"] + offset

    # plotname = f"11_Mar_R_vs_P_to_0.1_offset_{offset:3f}"


    # fig = plt.figure(0, figsize=(8,6.5))
    # ax1=plt.gca()
    # ax1.plot(i_arr* 1000 * avg_dict_low["v"] ,
    #          1000 * avg_dict_low["v"]/(i_arr),
    #          ".",# markersize=1,
    #          label = f"data, offset_{1000*offset:3f} [µA]")
    # # ax1.plot(t_avg_series-t_series[0],v_avg_series*1000,
    # #         label = f"moving average {mavg_len}s")

    # #plt.xticks(rotation = 45)

    # ax1.set_xlabel(r"Power [µW]")
    # ax1.set_ylabel(r"Resistance [$\Omega$]")

    # plt.grid(True)
    # plt.legend(shadow=True)
    # plt.tight_layout()
    # format_im = 'png' #'pdf' or png
    # dpi = 300
    # plt.savefig(plot_dir + plotname
    #             + '.{}'.format(format_im),
    #             format=format_im, dpi=dpi)
    # ax1.cla()


    ############## Old style
    # for key in index_dict.keys():
    #     print("key: ", key)
    #     # vba.plot_with_dates(data_dict, plotname = f"calib/{key}mA", 
    #     #                 start = index_dict[key][0], end = index_dict[key][-1]
    #     #                 )

    #     plot_calib(data_dict, plotname = f"calib/{key}mA", 
    #                     index_arr= index_dict[key]
    #                     )

    # key = "1.0"
    # plot_calib_kappa(data_dict, plotname = f"calib/{key}mA_kappa", 
    #                     index_arr= index_dict[key]
    #                     )

