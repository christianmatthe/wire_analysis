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
        'size'   : 20
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


flow_lst = [1, 0.2, 0.05
            ]
flow_arr = np.array(flow_lst)

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
flows_mainz = np.array([0.002,1,20])
ad_H =  np.array([58.139534883720934, 31.007751937984494 ,11.886304909560726])
ad_H_err = (np.array([67.95865633,36.69250646,20.93023256
]) - np.array(ad_H))

ad_H2 = np.array([64, 39, 28])
ad_H2_err = np.array([ 
            np.abs(np.array([0, 100]) - ad_H2[0]),
            np.abs(np.array([28.57142, 51.9480]) - ad_H2[1]),
            np.abs(np.array([8.571428, 66.4935]) - ad_H2[2]),
                    ]).T


#Results of "c*T" refitting done by alec
# tschersich_list = np.array(
#     [[0.002136235,	65.11627907],
#      [0.008851836,	46.51162791],
#      [0.099400037,	34.36692506],
#      [0.457797489,	12.14470284],
#     ]).T
# Read off from Tschersich 2000 fig 6 at 2200K
tschersich_list = np.array(
    [[0.002136235,	60],
     [0.008851836,	46.15],
     [0.099400037,	30],
     [0.457797489,	22.17],
    ]).T

tschersich_flows = tschersich_list[0]
tschersich_ads = tschersich_list[1]
print("tschersich_flows , tschersich_ads")
print(tschersich_flows , tschersich_ads)


from matplotlib.legend_handler import HandlerPatch
import matplotlib.patches as mpatches


from matplotlib.patches import (FancyArrowPatch, FancyArrow, Arrow)
from matplotlib.legend_handler import HandlerTuple
from matplotlib.legend_handler import HandlerPatch

include_1T = True
############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()

arrow_len = 8
if include_1T:
    # Wire det 1T analysis
    color_1T = "C8"
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
    label = r"Mass Spec., Tschersich",
    ms = 15, lw = 3,
    ls = "-",
    ecolor = "None",
    # lolims = True,
    # capsize = 6
    )
# MAss spec data:
ax1.errorbar(flows_mainz, ad_H, yerr = ad_H_err, fmt = "^",
             color = "C3", 
    label = r"Mass Spec., Mainz,  $H$ production",
    ms = 13, lw = 3,
    elinewidth = 2, capthick = 2,
    ls = "-",
    # lolims = True,
    capsize = 6
    )

ax1.errorbar(flows_mainz, ad_H2, yerr = ad_H2_err, fmt ="o",
            color = "C2", fillstyle= "none", mew = 3,
    label = r"Mass Spec., Mainz, $H_2$ reduction",
    ms = 15, lw = 3,
    # elinewidth = 0,
    # ecolor = "None",
    # lolims = True,
    capsize = 6
    )
##########




ax1.set_xlabel(r"Hydrogen Flow $\Phi_{\rm in}$ [sccm]")
ax1.set_ylabel(r"Dissociation Fraction $\alpha_{\rm dissoc} [\%]$")


# Create the upward arrow using FancyArrow
arrow = FancyArrow(0, 0, 0, 0.3,  # x, y, dx, dy (note: dy is positive for upward)
                   width=0.02, length_includes_head=True,
                   head_width=0.08, head_length=0.1,
                   color='green')


class MyArrowObject(object):
    def __init__(
        self,
        patch_color = "C0"):
        self.patch_color = patch_color
    pass

class MyArrowHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, 
                      handlebox):
        patch_color = orig_handle.patch_color


        x0, y0 = handlebox.xdescent, handlebox.ydescent
        hb_width, hb_height = handlebox.width, handlebox.height
        patch = mpatches.Rectangle([x0, y0], hb_width, hb_height, 
                                   facecolor= "None",
                                   edgecolor= "None", 
                                   hatch='xx', lw=3,
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        ###
        scale  = fontsize/16  # Original  was designed with fontize 16
        head_width=10  *scale
        head_length=10 / np.sqrt(2)  *scale
        width = 4  *scale

        # x = head_width/2
        x = hb_width/2  
        dx = 0
        dy = 10  *scale
        # y = -( dy) /2
        y  = -(head_length) /2  *scale

        #patch_color =  self.patch_color

        p = mpatches.FancyArrow(x, y, dx, dy,
                                head_width = head_width,
                                head_length= head_length,
                                width = width,
                                transform=handlebox.get_transform(),
                                fc= patch_color,
                                ec = patch_color,
                            )
        handlebox.add_artist(p)
        h_rect = 3  *scale
        p_rect = mpatches.Rectangle([x - head_width/2, y -1], 
                                    width = head_width, height = h_rect, 
                                    facecolor= patch_color,
                                transform=handlebox.get_transform()
                            )
        handlebox.add_artist(p_rect)
        return 

# plt.legend(handles = [MyArrowObject(patch_color  = "C0"),
#              MyArrowObject(patch_color  = color_1T)],
#             labels=["lower limit 3T", "lower limit 1T"],
#            handler_map={MyArrowObject: MyArrowHandler()})
handles, labels = plt.gca().get_legend_handles_labels()

if include_1T:
    handles = [MyArrowObject(patch_color  = "C0"),
                MyArrowObject(patch_color  = color_1T),
                handles[-3], handles[-2], handles[-1]]
    labels = ["Wire Det., lower limit, 3-T-Point", 
              "Wire Det., lower limit, 1-T-Point", 
            labels[-3], labels[-2], labels[-1]]
else:
    handles = [MyArrowObject(patch_color  = "C0"),

                handles[-3], handles[-2], handles[-1]]
    labels = ["Wire Det., lower limit, 3-T-Point", 

            labels[-3], labels[-2], labels[-1]]

plt.legend(handles = handles,
            labels= labels,
           handler_map={MyArrowObject: MyArrowHandler()})


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
if include_1T:
    plt.savefig(out_dir + "alpha_diss_from_A"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
else:
    plt.savefig(out_dir + "alpha_diss_from_A_base"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
############################################################################


# Plot total atom flux
sccm = 4.478 * 10**17

############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()

flow_lst_1T = np.array(flow_lst_1T)

# Wire det 1T analysis
arrow_len = 8
if include_1T:
    color_1T = "C8"
    ax1.errorbar(flow_lst_1T, np.array(ad_1T_lst) * 2*flow_lst_1T * sccm,
                yerr = np.abs(np.array(ad_1T_err_lst)).T * 2*flow_lst_1T * sccm,
                fmt = ",",
        label = r"wire det. T analysis",
        ms = 10, lw = 14,
        # lolims = True,
        capsize = 7,
        color = color_1T
        )

    ax1.errorbar(flow_lst_1T,  np.array(ad_1T_lst) * 2*flow_lst_1T * sccm, 
                yerr =  2*flow_lst_1T * sccm * (np.abs(np.array(ad_1T_err_lst)).T 
                + np.array(ad_1T_lst) * 0.2),
                # yerr =  2*flow_lst_1T * sccm * (1-np.array(ad_1T_lst)),
                fmt = ",",
        ms = 16, lw = 6,
        lolims = True,
        capsize = 7,
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

ax1.errorbar(flow_lst, np.array(ac_lst)* 2*flow_arr * sccm,
              yerr = 2*flow_arr * sccm * np.abs(np.array(ac_err_lst)).T,
                fmt = ",",
    label = r"propagating just $A$ fit errors",
    ms = 10, lw = 14,
    # lolims = True,
    capsize = 7
    )
# ax1.errorbar(flow_lst, np.array(ac_lst)*100, 
#              yerr = 100*np.abs(np.array(ac_err_lst)).T, fmt = ",",
#     ms = 18, lw = 5,
#     lolims = True,
#     #capsize = 6,
#     color = "C0"
#     )

if include_1T:
    ax1.errorbar(flow_lst,  np.array(ac_lst)* 2*flow_arr * sccm, 
            yerr = 2*flow_arr * sccm * (
                        np.abs(np.array(ac_err_lst)).T + np.array(ac_lst) * 0.2),
            fmt = ",",
        ms = 16, lw = 6,
        lolims = True,
        capsize = 7,
        color = "C0"
        )
else:
    ax1.errorbar(flow_lst,  np.array(ac_lst)* 2*flow_arr * sccm, 
            yerr = 2*flow_arr * sccm * (
                        np.abs(np.array(ac_err_lst)).T 
                        + np.array(ac_lst) * 0.3),
            fmt = ",",
        ms = 16, lw = 6,
        lolims = True,
        capsize = 7,
        color = "C0"
        )
##########

#Plot 100% dissoc
flows = np.linspace(0.001,30)
ax1.plot(flows, 2*flows * sccm, "--", color = "C7",
          label = "100\% Dissociation")
ax1.set_xlim([0.001,30])
########

####### tschersich:
ax1.errorbar(tschersich_flows, tschersich_ads *0.01*2* tschersich_flows * sccm,
              yerr = 0, fmt = "v",
             color = "k",
    label = r"Mass Spec., Tschersich",
    ms = 15, lw = 3,
    ls = "-",
    ecolor = "None",
    # lolims = True,
    # capsize = 6
    )
# MAss spec data:
ax1.errorbar(flows_mainz, ad_H *0.01*2* flows_mainz * sccm,
              yerr = ad_H_err *0.01*2* flows_mainz * sccm, fmt = "^",
             color = "C3", 
    label = r"Mass Spec., Mainz,  $H$ production",
    ms = 11, lw = 3,
    elinewidth = 2, capthick = 2,
    ls = "-",
    # lolims = True,
    capsize = 6
    )

yerrs = 0.01*2* flows_mainz * sccm* (ad_H2_err) 
ax1.errorbar(flows_mainz, ad_H2 *0.01*2* flows_mainz * sccm, 
             yerr = yerrs, fmt ="o",
            color = "C2", fillstyle= "none", mew = 3,
    label = r"Mass Spec., Mainz, $H_2$ reduction",
    ms = 13, lw = 3,
    # elinewidth = 0,
    # ecolor = "None",
    # lolims = True,
    capsize = 6
    )
#########


ax1.set_xlabel(r"Hydrogen Flow $\Phi_{\rm in}$ [sccm]")
ax1.set_ylabel(r"Atom flux $\Phi_{\rm at}$ [atoms / s]")



class MyArrowObject(object):
    def __init__(
        self,
        patch_color = "C0"):
        self.patch_color = patch_color
    pass

class MyArrowHandler(object):
    def legend_artist(self, legend, orig_handle, fontsize, 
                      handlebox):
        patch_color = orig_handle.patch_color


        x0, y0 = handlebox.xdescent, handlebox.ydescent
        hb_width, hb_height = handlebox.width, handlebox.height
        patch = mpatches.Rectangle([x0, y0], hb_width, hb_height, 
                                   facecolor= "None",
                                   edgecolor= "None", 
                                   hatch='xx', lw=3,
                                   transform=handlebox.get_transform())
        handlebox.add_artist(patch)
        ###
        scale  = fontsize/16  # Original  was designed with fontize 16
        head_width=10  *scale
        head_length=10 / np.sqrt(2)  *scale
        width = 4  *scale

        # x = head_width/2
        x = hb_width/2  
        dx = 0
        dy = 10  *scale
        # y = -( dy) /2
        y  = -(head_length) /2  *scale

        #patch_color =  self.patch_color

        p = mpatches.FancyArrow(x, y, dx, dy,
                                head_width = head_width,
                                head_length= head_length,
                                width = width,
                                transform=handlebox.get_transform(),
                                fc= patch_color,
                                ec = patch_color,
                            )
        handlebox.add_artist(p)
        h_rect = 3  *scale
        p_rect = mpatches.Rectangle([x - head_width/2, y -1], 
                                    width = head_width, height = h_rect, 
                                    facecolor= patch_color,
                                transform=handlebox.get_transform()
                            )
        handlebox.add_artist(p_rect)
        return 

# plt.legend(handles = [MyArrowObject(patch_color  = "C0"),
#              MyArrowObject(patch_color  = color_1T)],
#             labels=["lower limit 3T", "lower limit 1T"],
#            handler_map={MyArrowObject: MyArrowHandler()})
handles, labels = plt.gca().get_legend_handles_labels()


if include_1T:
    handles = [handles[0],
               MyArrowObject(patch_color  = "C0"),
                MyArrowObject(patch_color  = color_1T),
                handles[-3], handles[-2], handles[-1]]
    labels = [labels[0],"Wire Det., lower limit, 3-T-Point", 
              "Wire Det., lower limit, 1-T-Point", 
            labels[-3], labels[-2], labels[-1]]
else:
    handles = [handles[0],
               MyArrowObject(patch_color  = "C0"),

                handles[-3], handles[-2], handles[-1]]
    labels = [labels[0], "Wire Det., lower limit, 3-T-Point", 

            labels[-3], labels[-2], labels[-1]]


plt.legend(handles = handles,
            labels= labels,
           handler_map={MyArrowObject: MyArrowHandler()},
           fontsize = 16)


# plt.legend(shadow=True)
plt.tight_layout()
#plt.grid(True)

#ax1.set_ylim([0,100])
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
ax1.yaxis.set_minor_locator(AutoMinorLocator())

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(which = "major", color = "0.5")
ax1.grid(which="minor", color = "0.9")


format_im = 'png' #'pdf' or png
dpi = 300
if include_1T:
    plt.savefig(out_dir + "atom_flux_from_A"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi)
else:
    plt.savefig(out_dir + "atom_flux_from_A_base"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
###################






