# Based on 2025-03-05_6T-point-method_combinations.py
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit
from itertools import chain, combinations
from tqdm import tqdm

import wire_analysis as wa
from wire_analysis.utils import (load_json_dict,save_json_dict,
                                  load_extractor_dict_json)
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient, TC_to_T_Hack, Cv)
from wire_analysis.flow_on_off_cycle_analysis import (
                                        sort_by_z_list)

#plot Options
# Adjusted for large size in multi figure array
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 20
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)
# plt.rcParams.update({
#     "text.usetex": True,
#     "font.family": "sans-serif",
#     "font.sans-serif": "Helvetica",
# })
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": "Computer Modern",
    "font.weight":"bold"
})

base_dir = (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\\"
            + r"scripts\\")
work_dir_list = ["2025-03-06_multi-point(interp)_method_10sccm",
                "2025-03-05_6-point_method_1sccm",
                "2025-05-10_0.2sccm_cut_data_dch",
                "2025-05-10_0.05sccm_Cv _cut_data_dch",
                "2025-07-11_1Temp_leffs",
                "2025-07-11_1Temp_leffs"
                
]
flow_ind_lst = ["10sccm",
                "1sccm",
                "02sccm",
                "005sccm",
                "001sccm",
                "0002sccm"
                ]

for i_w,wd in enumerate(work_dir_list):
    # Run through all the analysis dicts in question
    work_dir = base_dir + wd +"/run_dicts/"
    out_dir = "./output/"
    os.makedirs(out_dir, exist_ok=True)

    # function for extracting data lists
    def p_data_plot_dict(filename):
        beamfit = wa.Beamfit(
                    run_dict_path = work_dir + filename)
        rd = beamfit.run_dict
        ######################################
        # Load from newly reformatted result dict files
        extractor_dict_unsorted = load_extractor_dict_json(
                                rd["extractor_dict_path"])
        # ##############
        z_array_unsorted = np.array(rd["z_list_unsorted"])

        # Sort these:
        z_arr, extractor_dict = sort_by_z_list(z_array_unsorted,
                                            extractor_dict_unsorted)

        P_arr = extractor_dict["p_arr"]
        P_err_arr = extractor_dict["p_err_arr"]
        plot_dict = {}
        plot_dict["z_arr"] = z_arr
        plot_dict["p_arr"] = P_arr
        plot_dict["p_err_arr"] = P_err_arr
        # #HACK, artificially inflate errors to see what happens
        # plot_dict["p_err_arr"] = P_err_arr + 0.05
        return plot_dict

    # plot data into shared plot   
    #######
    # TC_lst = [200, 300, 390, 475]
    indicator_list = ["720TC", "715TC", "475TC", "390TC", "300TC",
                      "200TC",  "0A"]
    # for TC_lst 44 is a HACK for 0 A ("room temp" =  44 = 297.7K)
    TC_lst = [720,715, 475, 390, 300, 200, 44] 
    # # Rearrange to make color code uniform HACK
    # indicator_list = ["720TC", "715TC", "390TC",  "0A",
    #                   "475TC", "300TC", "200TC", ]
    # # for TC_lst 44 is a HACK for 0 A ("room temp" =  44 = 297.7K)
    # TC_lst = [720,715, 390, 44, 
    #            475, 300, 200] 
    # T_lst = [TC_to_T_Hack(TC) for TC in TC_lst]

    #find subselection of indicators that exissts for this flow:
    filename_list = []
    T_lst = []
    for j, indicator in enumerate(indicator_list):
        filename = flow_ind_lst[i_w] + "_" + indicator + ".json" 
        if os.path.isfile((work_dir + filename)):
            filename_list.append(filename)
            T_lst.append(TC_to_T_Hack(TC_lst[j]))
        else:
            print(filename, "does not exists in ", wd, work_dir)
            pass

    pd_dict = {}
    #HACK color list to make same color for same temperature
    if flow_ind_lst[i_w] == "1sccm": #Just becasue 1sccm has more datasets
        color_list = ["C0", "C3", "C1", "C4", "C5", "C2"]
    else:
        color_list = ["C0", "C1", "C2"]
    # Start plot:
    plotname = flow_ind_lst[i_w] + "_" + "multi_run_power"
    fig = plt.figure(0, figsize=(8,5.5), dpi =300)
    ax1=plt.gca()
    x_label = r"$z_{pos}$ [mm]"
    for i,filename in enumerate(filename_list):
        print("filename:", filename)
        pd = p_data_plot_dict(filename)
        ## HACK Data cuts for 005 measurement
        if flow_ind_lst[i_w] == "005sccm":
            cut_list = [0,5,10,15]
            for key in pd.keys():
                pd[key] = np.delete(pd[key],cut_list)

        pd_dict[indicator_list[i]] = pd
        pd_dict[indicator_list[i]]["index"] = i
        pd_dict[indicator_list[i]]["T"] = T_lst[i]
        ### ax1
        # ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"],
        #  fmt = ".",
        #             label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))
        ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"],
                      fmt = ".",
            label = (r"$P_{\rm{meas}}(\rm{T}$ $\approx$ "+f"{T_lst[i]:.0f}K)"), 
            markersize =10,
            color = color_list[i]
            )


    ax1.set_ylabel(r"Detected Power [µW]")
    ax1.set_xlabel(x_label)

    ax1.grid(True)
    ax1.legend(shadow=True, fontsize = 13)
    # ax1.tight_layout()


    fig.tight_layout()

    format_im = 'png' #'pdf' or png
    dpi = 300
    plt.savefig(out_dir + plotname
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi, bbox_inches="tight")
    # plt.show()
    ax1.cla()
    fig.clf()
    plt.close()