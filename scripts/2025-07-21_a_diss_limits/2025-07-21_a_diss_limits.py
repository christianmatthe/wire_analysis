import os
import numpy as np
import matplotlib.pyplot as plt

import wire_analysis as wa
from wire_analysis.utils import load_json_dict, load_extractor_dict_json
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient, TC_to_T_Hack, ac_from_Abg, a_diss_from_A)
from wire_analysis.beamshape import (calc_norm_factor)


#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)


# Run through all the analysis dicts in question
work_dir = "./run_dicts/"
out_dir = "./output/"
os.makedirs(out_dir, exist_ok=True)

# # Do all fits
# for filename in os.listdir(work_dir):
#     print("filename:", filename)
#     beamfit = wa.Beamfit(
#         run_dict_path = work_dir + filename)
#     # HACK to change out dir without editing the files
#     # TODO If you do this implement feature to retroactively change out_dir in
#     # run_dict
#     ############### TODO Implement as base function
#     name, file_extension = os.path.splitext(filename)
#     beamfit.run_dict["out_dir_base"] = os.path.abspath("./output/") + os.sep 
#     beamfit.out_dir = beamfit.run_dict["out_dir_base"] + name + os.sep
#     beamfit.save_json_run_dict()
#     os.makedirs(beamfit.out_dir, exist_ok=True)
#     ############### END TODO


#     beamfit.default_fit()
# ##### END Do all fits


# for filename in os.listdir(work_dir):
#     print("filename:", filename)
#     path = work_dir + filename
#     rd = load_json_dict(path)
#     print(rd["fit_result"]["l_eff"])


flow_lst = [1, 0.2, 0.05,
            ]

flow_lst_1T = [10,
               1, 0.2, 0.05,
               0.01, 0.002
               ]
flow_lst_1T_str = ["10",
               "1", "02", "005",
               "001", "0002"
               ]
flow_lst_1T_TC_str = ["720",
               "720", "715", "715",
               "715", "715"
               ]

#NOTE Some of zhese files are not protected, in the sense, that they could 
# contain fit results that ate not of the 715TC, 390TC, 0A (2200, 1250, 300K)
# format. In particular the 1 sccm file is froma source where it could be 
# overwritten with other temperature combinations
file_list = [
    (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\scripts\\"
    + r"2025-05-10_6-point_method_1sccm_dch\output\1sccm_720TC_penumbra\\"
    + r"1sccm_720TC_penumbra.json"),

    (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\scripts\\"
     + r"2025-05-10_0.2sccm_cut_data_dch\output\02sccm_715TC_penumbra\\"
     + r"02sccm_715TC_penumbra.json"),

    (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\scripts\\"
     + r"2025-05-10_0.05sccm_Cv _cut_data_dch\output\005sccm_715TC_penumbra\\"
     + r"005sccm_715TC_penumbra.json"),

    #  #HACK 0.01sccm not  done  with  the same technique

    # (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\scripts\\"
    #    + r"2025-07-11_1Temp_leffs\output\001sccm_715TC\001sccm_715TC.json" ),
    # (r"C:\Users\Christian\Documents\StudiumPhD\python\wire_analysis\scripts\\"
    #    + r"2025-07-11_1Temp_leffs\output\0002sccm_715TC\0002sccm_715TC.json" ),
    ]

file_list_1T = [(r"C:\Users\Christian\Documents\StudiumPhD\python\\"
                + r"wire_analysis\scripts\2025-07-11_1Temp_leffs\output\\"
                + f"{s}sccm_{flow_lst_1T_TC_str[i]}TC\\"
                + f"{s}sccm_{flow_lst_1T_TC_str[i]}TC.json" ) 
                for i,s in enumerate(flow_lst_1T_str)
                ]



ac_lst = []
ac_err_lst = []
ac_err_leff_lst = []
for i,TC in enumerate(flow_lst):
    path = file_list[i]
    rd = load_json_dict(path)


    ad = a_diss_from_A(
        A = rd["fit_result"]["A"],
        flow = flow_lst[i] * 4.478 * 10**17 #sccm
        )
    ad_err = [ad - a_diss_from_A(
                A = rd["fit_result"]["A"] - rd["fit_result_errors"]["A"],
                flow = flow_lst[i] * 4.478 * 10**17 #sccm
                ), 
              ad - a_diss_from_A(
                A = rd["fit_result"]["A"] + rd["fit_result_errors"]["A"],
                flow = flow_lst[i] * 4.478 * 10**17 #sccm
                )
            ]

    ac_lst.append(ad)
    ac_err_lst.append(ad_err)
#     ac_err_leff_lst.append(ac_err_leff)
print("ad_lst", ac_lst)
print("ad_errs", ac_err_lst)
print("errs_abs",  np.abs(np.array(ac_err_lst)).T)

ad_1T_lst = []
ad_1T_err_lst = []
for i,flow in enumerate(flow_lst_1T):
    path = file_list_1T[i]
    rd = load_json_dict(path)


    ad = a_diss_from_A(
        A = rd["fit_result"]["A"],
        flow = flow * 4.478 * 10**17 #sccm
        )
    ad_err = [ad - a_diss_from_A(
                A = rd["fit_result"]["A"] - rd["fit_result_errors"]["A"],
                flow = flow * 4.478 * 10**17 #sccm
                ), 
              ad - a_diss_from_A(
                A = rd["fit_result"]["A"] + rd["fit_result_errors"]["A"],
                flow = flow * 4.478 * 10**17 #sccm
                )
            ]

    ad_1T_lst.append(ad)
    ad_1T_err_lst.append(ad_err)
#     ac_err_leff_lst.append(ac_err_leff)
print("ad_lst_1T", ad_1T_lst)
print("ad_errs_1T", ad_1T_err_lst)




# TODO Max power to max power comparison at different  temperatures
# ############ Plot results with errors
# fig = plt.figure(0, figsize=(8,6.5))
# ax1= plt.gca()
# # ax1.errorbar(T_lst, ac_lst, yerr = comb_err, fmt = ".",
# #     label = r"p and $l_{\rm{eff}}$ errs",
# #     ms = 18, lw = 5,
# #     )

# # ax1.errorbar(flow_lst, ac_lst, yerr = ac_T_errs, xerr = T_errs, fmt = ".",
# #     label = r"including  2\% $\Delta T$",
# #     ms = 18, lw = 5,
# #     )
# ax1.errorbar(flow_lst, ac_lst, yerr = np.abs(np.array(ac_err_lst)).T, fmt = ".",
#     label = r"propagating just $A$ fit errors",
#     ms = 18, lw = 5,
#     )



# ax1.set_xlabel(r"Flow [sccm]")
# ax1.set_ylabel(r"$\alpha_{\rm dissoc}$")


# # h, l = ax1.get_legend_handles_labels()
# # select = [1,0]
# # ax1.legend([h[i] for i in select], [l[i] for i in select],
# #                #shadow = True,
# #                )

# plt.legend(shadow=True)
# plt.tight_layout()
# #plt.grid(True)

# ax1.set_xscale('log')
# ax1.grid(which = "major")
# ax1.grid(which="minor", color = "0.9")


# format_im = 'png' #'pdf' or png
# dpi = 300
# plt.savefig(out_dir + "alpha_diss_from_A_basic"
#             + '.{}'.format(format_im),
#             format=format_im, dpi=dpi)
# # plt.show()
# ax1.cla()
# fig.clf()
# plt.close()
# ###################


#Digitized data from Alecs plot
#https://3.basecamp.com/3700981/buckets/3107037/uploads/8400184943
flows_mainz = [0.002,1,20]
ad_H =  [58.139534883720934, 31.007751937984494 ,11.886304909560726]
ad_H_err = (np.array([67.95865633,36.69250646,20.93023256
]) - np.array(ad_H))

ad_H2 = [65.11627907, 41.08527132, 29.97416021]

tschersich_list = np.array(
    [[0.002136235,	65.11627907],
     [0.008851836,	46.51162791],
     [0.099400037,	34.36692506],
     [0.457797489,	12.14470284],
    ]).T
tschersich_flows = tschersich_list[0]
tschersich_ads = tschersich_list[1]
print("tschersich_flows , tschersich_ads")
print(tschersich_flows , tschersich_ads)


############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()



# Wire det 1T analysis
arrow_len = 8
color_1T = "C7"
ax1.errorbar(flow_lst_1T, np.array(ad_1T_lst)*100,
              yerr = 100*np.abs(np.array(ad_1T_err_lst)).T, fmt = ",",
    label = r"wire det. T analysis",
    ms = 10, lw = 16,
    # lolims = True,
    capsize = 8,
    color = color_1T
    )

ax1.errorbar(flow_lst_1T, np.array(ad_1T_lst)*100, 
             yerr = 100*np.abs(np.array(ad_1T_err_lst)).T + arrow_len,
               fmt = ",",
    ms = 16, lw = 6,
    lolims = True,
    capsize = 8,
    color = color_1T
    )
#######
########### For Wire det data
# ax1.errorbar(flow_lst, np.array(ac_lst)*100,
#               yerr = 100*np.abs(np.array(ac_err_lst)).T, fmt = ",",
#     label = r"propagating just $A$ fit errors",
#     ms = 10, lw = 3,
#     # lolims = True,
#     capsize = 10
#     )

ax1.errorbar(flow_lst, np.array(ac_lst)*100,
              yerr = 100*np.abs(np.array(ac_err_lst)).T, fmt = ",",
    label = r"propagating just $A$ fit errors",
    ms = 10, lw = 16,
    # lolims = True,
    capsize = 8
    )
# ax1.errorbar(flow_lst, np.array(ac_lst)*100, 
#              yerr = 100*np.abs(np.array(ac_err_lst)).T, fmt = ",",
#     ms = 18, lw = 5,
#     lolims = True,
#     #capsize = 6,
#     color = "C0"
#     )

ax1.errorbar(flow_lst, np.array(ac_lst)*100, 
             yerr = 100*np.abs(np.array(ac_err_lst)).T + arrow_len, fmt = ",",
    ms = 16, lw = 6,
    lolims = True,
    capsize = 8,
    color = "C0"
    )
##########

####### tschersich:
ax1.errorbar(tschersich_flows, tschersich_ads, yerr = 0, fmt = "v",
             color = "k",
    label = r"Tschersich",
    ms = 18, lw = 3,
    ls = "-",
    # lolims = True,
    # capsize = 6
    )
# MAss spec data:
ax1.errorbar(flows_mainz, ad_H, yerr = ad_H_err, fmt = "^",
             color = "C3", 
    label = r"Mainz H",
    ms = 18, lw = 3,
    ls = "-",
    # lolims = True,
    capsize = 6
    )

ax1.errorbar(flows_mainz, ad_H2, yerr = 0, fmt = "o",
            color = "C2", fillstyle= "none", mew = 3,
    label = r"Mainz H2",
    ms = 15, lw = 5,
    # lolims = True,
    # capsize = 6
    )
##########




ax1.set_xlabel(r"Hydrogen Flow [sccm]")
ax1.set_ylabel(r"$\alpha_{\rm dissoc} [\%]$")


# h, l = ax1.get_legend_handles_labels()
# select = [1,0]
# ax1.legend([h[i] for i in select], [l[i] for i in select],
#                #shadow = True,
#                )
from matplotlib.lines import Line2D
# lower_cap = Line2D([0], [0], color='C0', marker=r'$\uparrow$',
#                     linestyle='None',
#                    markersize=12)


from matplotlib.patches import (FancyArrowPatch, FancyArrow, Arrow)
from matplotlib.legend_handler import HandlerTuple
from matplotlib.legend_handler import HandlerPatch
# Create the upward arrow using FancyArrow
arrow = FancyArrow(0, 0, 0, 0.3,  # x, y, dx, dy (note: dy is positive for upward)
                   width=0.02, length_includes_head=True,
                   head_width=0.08, head_length=0.1,
                   color='green')

# Add legend
# plt.legend([legend_handle], ['Lower limit'],
#                       handler_map={
#                tuple: HandlerTuple(ndivide=None),
#                FancyArrow: HandlerFancyArrow()
#            }, loc='upper right')

# plt.legend([arrow], ['Lower limit'],
#         handler_map={FancyArrow : HandlerPatch(patch_func=make_legend_arrow),
#            }, loc='upper right')

from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches

def make_legend_arrow(legend, orig_handle,
                      xdescent, ydescent,
                      width, height, fontsize):

    head_width=10
    head_length=10 / np.sqrt(2)
    width = 5

    x = head_width/2
    dx = 0
    dy = 10
    # y = -( dy) /2
    y  = -(head_length) /2

    p = mpatches.FancyArrow(x, y, dx, dy,
                            head_width = head_width,
                            head_length= head_length,
                            width = width,
                           )
    
    return p

arrows = [mpatches.FancyArrow(0,0,0,0 # dummy placeholders
                              , fc='C0',
                               ec = "C0"),
]

plt.legend(handles=arrows, 
           labels=["lower limit"], 
           handler_map={mpatches.FancyArrow : HandlerPatch(
               patch_func=make_legend_arrow),
                    })

# plt.legend(shadow=True)
plt.tight_layout()
#plt.grid(True)

ax1.set_ylim([0,100])
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.set_xscale('log')
ax1.grid(which = "major", color = "0.5")
ax1.grid(which="minor", color = "0.9")


format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + "alpha_diss_from_A"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
###################






