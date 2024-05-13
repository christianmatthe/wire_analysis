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


# Change working patch to directory layer above
# import sys
#sys.path.insert(0,(os.path.dirname(os.path.abspath(__file__)) + os.sep /..)
# cwd = os.getcwd()    
# print(cwd + os.sep + ".." + os.sep)         
# os.chdir(cwd + os.sep + ".." + os.sep)
cwd = os.getcwd()
print("cwd:", os.getcwd())

# # import functions from other file
# import sys
# sys.path.insert(0,cwd)
# import Voltage_base_analysis as vba
# import calib_analysis as ca
# import pressure_analysis as pa
#from flow_on_off_cycle_analysis_2 import *
#HACK horriblle frankenstein extraction. Use new extraction with old code mix
from wire_analysis.flow_on_off_cycle_analysis import *

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)

plot_dir = (cwd + os.sep 
            + "output/flow_on_off/")
# data_dir = cwd + os.sep + "data/"
os.makedirs(plot_dir, exist_ok=True)
#######################

if __name__ =="__main__": 
# ########################################
########## 2023-04-24  2023-04-21_1sccm_15A_TC_z-scan_jf_wire
    sc_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep + "../" 
         + "SC_downloads/")
    HABS_dir = (os.path.dirname(os.path.abspath(__file__)) + os.sep + "../" 
         + "SC_downloads/HABS_power/")
    os.makedirs(HABS_dir, exist_ok=True)
    run_name = "2023-04-21_1sccm_15A_TC_z-scan_jf_wire_recalib"
    data_name = "2023-04-21_1sccm_15A_TC_z-scan_jf_wire"


    # z_list = [-11,-10,-9,-8,-7,-6,-5,-4,-3.5,-3,
    #                -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2,
    #                2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 7, 
    #                 8,
    #                 9, 10, 11, 12, 14, 16, 18, 20
    #                ]
    # #arbitrarily chose 30mm as distanc eof cracker to wire plane
    # # Measured as 35.17mm in CAD drawing, not actual distance
    # dist = 35.17
    # alpha_list = [np.tan(z/dist) * 180/np.pi for  z in z_list]

    # global_start_date = dt.datetime(2023, 4, 21,21, 29, 56,
    #                             tzinfo=dt.timezone(dt.timedelta(hours=1)))
    # # Second annd milisecond fudges here are to realign unaccounted slow
    # # controls delay
    # # an additional 1/(24*60)*(5*60)*1000 = 208.3 milliseconds
    # # to account for 1 min per day residual drift
    # flow_resynch = 780 + (1/(24*60)*(5*60)*1000) 
    # flow_interval = dt.timedelta(minutes = 5, milliseconds=flow_resynch)
    # z_resynch =  flow_resynch * 1e-3 *24
    # z_interval = dt.timedelta(minutes = 120, seconds=z_resynch )
    # # cut first 10 minutes both to eliminate synchronization drift and the
    # # fastest part of the tempperature change
    # z_front_cut = dt.timedelta(minutes = 20)
    # # eliminate last 5 minutes where the z-translater tennds to desynch
    # z_rear_cut = dt.timedelta(minutes = 10)
    # denoise = False
    # denoise_stringt = ""
    # if denoise == True:
    #     denoise_string = "_denoised"
    # else:
    #     denoise_string = ""

    # z_list_unsorted = z_list

    # #  #### HACK START Comment oot to laod instead of calculating
    # # #load wire data
    # # vba.prep_data_slowdash("../SC_downloads/Wire/"+ data_name + ".json"
    # #                  , run_name
    # #                  )
    # # data_dict = vba.load_data(run_name)

    # # # extraction dict:
    # # ext_dict = {}
    # # # fudge by cutting malfunctioning last z measurement
    # # for i,z in enumerate(z_list):
    # #     start_time = time.time()
    # #     if i == 0:
    # #         ext = Extractor(data_dict, 
    # #             start_date = global_start_date + i  * z_interval
    # #                          + flow_interval * 4, 
    # #                          # cut  first 20 minutes that were not logged 
    # #             end_date = global_start_date + (i+1) * z_interval
    # #                        - flow_interval ,

    # #             time_interval = flow_interval,
    # #             run_name = run_name + os.sep 
    # #                        + "{}_{:.2f}mm{}".format(i,z,denoise_string),
    # #             plot_dir = plot_dir,
    # #             utc_offset = 2,
    # #             initial_state = 0,
    # #             front_crop = dt.timedelta(minutes = 2,seconds  =30),
    # #             rear_crop = dt.timedelta(minutes  =0, seconds  =45),
    # #             denoise = denoise
    # #             )
    # #     else:
    # #         ext = Extractor(data_dict, 
    # #             start_date = global_start_date + i * z_interval
    # #                         + z_front_cut, 
    # #             end_date = global_start_date + (i+1) * z_interval
    # #                        - flow_interval - z_rear_cut ,
    # #                        # cut last extraneous "B" cycle since we only care 
    # #                        # fully completed ABA cycles
    # #             time_interval = flow_interval,
    # #             run_name = run_name + os.sep 
    # #                        + "{}_{:.2f}mm{}".format(i,z,denoise_string),
    # #             plot_dir = plot_dir,
    # #             utc_offset = 2,
    # #             initial_state = 0,
    # #             front_crop = dt.timedelta(minutes = 2,seconds  =30),
    # #             rear_crop = dt.timedelta(minutes  =0, seconds  =45),
    # #             denoise = False
    # #             )
    # #     # ext.basic_extraction()
    # #     ext.quad_extraction()
    # #     ext_dict[i] = {"z":z,
    # #                    "extractor":ext
    # #                     }
    # #     run_time = (time.time()- start_time)/60
    # #     print("time taken for step:", run_time, "min" )
    # #     start_time = time.time()

    
    # # save_dict(ext_dict, plot_dir + run_name + os.sep 
    # #                             + "ext_dict")
    # # data_sets = [[ext_dict[i]["extractor"]] 
    # #             for i,z in enumerate(z_list)
    # #             ]
    # # save_dict(data_sets, plot_dir + run_name + os.sep 
    # #                             + "data_sets")
    # # ##### HACK END  load  instead          

    # # Analysis           
    # data_sets = load_dict(plot_dir + run_name + os.sep 
    #                             + "data_sets")
    # # print("data_sets:",data_sets)

    data_dir = (r"C:\\Users\\Christian\\Documents\\StudiumPhD\\python\\"
                + "Keysight-DMM-34461A\\analysis\\output\\flow_on_off\\")
    ext_dict = load_dict(data_dir + run_name + os.sep 
                                + "ext_dict")
    result_dict = make_result_dict(ext_dict)
    print(result_dict)
    print(result_dict.items())

    # Plot directly from saved ext object:
    # Chose 18th position = z = 1.5mm  for paper
    for i in range(37):
        ext = ext_dict[i]["extractor"]
        z = ext_dict[i]["z"]
        plot_path = (plot_dir + run_name + os.sep 
                     + f"{i}_{z:.2f}mm" + os.sep
                     + "paper" )
        os.makedirs(plot_path)
        ext.plot_all_ABA_fit_paper(
                            plot_path = plot_path,
                                )

#     #HACK Dump ext_dict to json. (Should be done as part of "extract" function
#     # in flow_on_off_cycle_analysis)
#      #connvert to json compatible dict
#     result_dict_unsorted_lists = {}
#     for key,val in result_dict.items():
#         # every numpy array must become a list
#         result_dict_unsorted_lists[key] = (
#             result_dict[key].tolist())
#     # Dump to json
#     with open((plot_dir + run_name + os.sep + "extractor_dict" + ".json"), 
#                 'w', encoding='utf-8') as f:
#         json.dump(result_dict_unsorted_lists, 
#                     f, ensure_ascii=False, indent=4)


#     power_plot(result_dict,
#                 x_lst= z_list,
#                 x_label = r"z-position [mm]",
#                 y_label = r"Detected Power [µW]",
#                 i_split = 50,
#                 plot_path= plot_dir + run_name + os.sep 
#                                 + "Detected_P_over_z_unsorted",
#                     flow_list = [1],
#                 )
    
#     # Plot individual  scans   
#     results_array = np.array([[result_dict["p_arr"][0:8], 
#                       result_dict["p_err_arr"][0:8]],
#                         [result_dict["p_arr"][8:16], 
#                       result_dict["p_err_arr"][8:16]],
#                         [result_dict["p_arr"][16:24], 
#                       result_dict["p_err_arr"][16:24]],
#                         [result_dict["p_arr"][24:32], 
#                       result_dict["p_err_arr"][24:32]],
#                         [result_dict["p_arr"][32:39], 
#                       result_dict["p_err_arr"][32:39]]

#     ])

#     #     z_list = [-11, -6, -1.5,  1,  3.5,
# #                     7, 12,  20, 
# #                    -10, -5,
# #                     -1, 1.5, 4, 8, 13, 
# #                     18, 
# #                       -9, -4, -0.5, 2,
# #                     4.5, 9, 14, 19,
# #                    -8,
# #                     -3, 0, 2.5, 5, 10,
# #                     15, 17,
# #                    -7, -2, 0.5,
# #                     3, 6, 11, 16
# #                     ]

#     power_plot_mult(results_array,
#                 x_arr = [z_list[0:8], z_list[8:16],
#                           z_list[16:24],
#                          z_list[24:32],z_list[32:39]],
#                 x_label = r"z-position [mm]",
#                 y_label = r"Detected Power [µW]",
#                 #i_split = 17,
#                 plot_path= plot_dir + run_name + os.sep 
#                                 + "Detected_P_over_z_scan_split",
#                 label_list = ["Scan1", "Scan2", "Scan3", "Scan4", "Scan5"],

#                 )
    
#     z_list, result_dict = sort_by_z_list(z_list_unsorted, result_dict)
#     print("sorted_results:", {"z_list": np.array(z_list)}, result_dict)


#     power_plot(result_dict,
#                 x_lst= z_list,
#                 x_label = r"z-position [mm]",
#                 y_label = r"Detected Power [µW]",
#                 i_split = 50,
#                 plot_path= plot_dir + run_name + os.sep 
#                                 + "Detected_P_over_z",
#                     flow_list = [1],
#                 )

#      # Measured as 35.17mm in CAD drawing, not actual distance
#     dist = 35.17
#     alpha_list = [np.tan(z/dist) * 180/np.pi for  z in z_list]
    
#     power_plot(result_dict,
#                 x_lst= alpha_list,
#                 x_label = r"$\alpha$" + r"-position [deg]",
#                 y_label = r"Detected Power [µW]",
#                 i_split = 50,
#                 plot_path= plot_dir + run_name + os.sep 
#                                 + "Detected_P_over_angle",
#                     flow_list = [1],
#                 )
    
#     # Plot k (µW per ohm)
#     k_arr = result_dict["µW_per_ohm"]

#     x_lst = z_list
#     #plot
#     fig = plt.figure(0, figsize=(8,6.5))
#     ax1=plt.gca()
#     ax1.plot(x_lst, k_arr,
#              ".",
#     #label = f"{flow_list[i]} sccm, up",
#     #color = f"C{2*(i//2)}",
#     )
#     ax1.set_xlabel(r"z-position [mm]")
#     ax1.set_ylabel("k [µW/Ohm]")

#     plt.grid(True)
#     plt.legend(shadow=True)
#     plt.tight_layout()
#     format_im = 'png' #'pdf' or png
#     dpi = 600
#     plt.savefig(plot_dir + run_name + os.sep 
#                                 + "k_over_z"
#                     + '.{}'.format(format_im),
#                     format=format_im, dpi=dpi)
#     ax1.cla()


###############################
# # Save ext_dict in exportable  formmats based on "legacy load"
#     ext_path = plot_dir + run_name + os.sep + "ext_dict"
#     if os.path.isfile(ext_path + ".pkl"):
#         result_dict_unsorted = legacy_load_pkl(pkl_path = ext_path)
#         save_dict(result_dict_unsorted, plot_dir + run_name + os.sep
#                     + "result_dict")
#         #connvert to json compatible dict
#         result_dict_unsorted_lists = {}
#         for key,val in result_dict_unsorted.items():
#             # every numpy array must become a list
#             result_dict_unsorted_lists[key] = (
#                 result_dict_unsorted[key].tolist())
#         # Dump to json
#         with open((plot_dir + run_name + os.sep  + "extractor_dict" + ".json"), 
#                     'w', encoding='utf-8') as f:
#             json.dump(result_dict_unsorted_lists, 
#                         f, ensure_ascii=False, indent=4)