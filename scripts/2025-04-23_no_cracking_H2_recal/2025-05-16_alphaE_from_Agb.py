import os
import numpy as np
import matplotlib.pyplot as plt

import wire_analysis as wa
from wire_analysis.utils import load_json_dict, load_extractor_dict_json
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient, TC_to_T_Hack, ac_from_Abg)
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

# TODO Time to implement the accommodation_coeffcient methods
# folowing "accomodation_coeffficient_H2_2023-12-06.ipynb"
# We jsut needed all that, so we could call up  fit results and the 
# extractor_dict by filename :D

# HACK The below is not a clean solution (hiddenin the fucntions it uses)
TC_lst = [200, 300, 390, 475]
T_lst = [TC_to_T_Hack(TC) for TC in TC_lst]
print(T_lst)
# The conversion to T is  done by eye based on a plot by Max. For annyone who 
# has nto noticed yet: that is a major  HACK
# T_lst = [750,1000,1250,1500]
ac_lst = []
ac_err_lst = []
ac_err_leff_lst = []
ac_alt_lst = []
ac_T_err_lst = []
for i,TC in enumerate(TC_lst):
    path = work_dir + f"1sccm_{TC}TC.json"
    rd = load_json_dict(path)


    ac = ac_from_Abg(
        Abg = rd["fit_result"]["A"],
        T = T_lst[i],
        flow = 1 * 4.478 * 10**17 #sccm
        )[0]
    ac_err = [ac - ac_from_Abg(
                Abg = rd["fit_result"]["A"] - rd["fit_result_errors"]["A"],
                T = T_lst[i],
                flow = 1 * 4.478 * 10**17 #sccm
                )[0], 
              ac - ac_from_Abg(
                Abg = rd["fit_result"]["A"] + rd["fit_result_errors"]["A"],
                T = T_lst[i],
                flow = 1 * 4.478 * 10**17 #sccm
                )[0]
            ]

    ac_Tp = ac_from_Abg(
        Abg = rd["fit_result"]["A"],
        T = T_lst[i]*1.02,
        flow = 1 * 4.478 * 10**17 #sccm
        )[0]
    ac_Tm = ac_from_Abg(
        Abg = rd["fit_result"]["A"],
        T = T_lst[i]*0.98,
        flow = 1 * 4.478 * 10**17 #sccm
        )[0]
    #WRONG
    # ac_T_err = [ac_Tp - ac_from_Abg(
    #             Abg = rd["fit_result"]["A"] - rd["fit_result_errors"]["A"],
    #             T = T_lst[i]*0.98,
    #             flow = 1 * 4.478 * 10**17 #sccm
    #             )[0], 
    #           ac_Tm - ac_from_Abg(
    #             Abg = rd["fit_result"]["A"] + rd["fit_result_errors"]["A"],
    #             T = T_lst[i]*1.02,
    #             flow = 1 * 4.478 * 10**17 #sccm
    #             )[0]
    #         ]
    ac_lst.append(ac)
    ac_err_lst.append(ac_err)
    ac_T_err_lst.append(ac_err)
#     ac_err_leff_lst.append(ac_err_leff)
print("ac_lst", ac_lst)
# comb_err = np.transpose(np.array([
#             [np.abs(lst1[i][0] - lst2[i]),lst1[i][1] + lst2[i]]
#              for i in range(len(lst2))]))
comb_err = [ac_err_lst[i][0] for i in range(len(ac_err_lst))]
print("comb_err", comb_err)
T_errs = [i*0.02 for i in T_lst]
ac_T_errs = [ac_err_lst[i][0] + (np.abs(ac -ac_Tp)) 
             for i in range(len(ac_err_lst))]
# TODO Max power to max power comparison at different  temperatures
############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()
# ax1.errorbar(T_lst, ac_lst, yerr = comb_err, fmt = ".",
#     label = r"p and $l_{\rm{eff}}$ errs",
#     ms = 18, lw = 5,
#     )

ax1.errorbar(T_lst, ac_lst, yerr = ac_T_errs, xerr = T_errs, fmt = ".",
    label = r"including  2\% $\Delta T$",
    ms = 18, lw = 5,
    )
ax1.errorbar(T_lst, ac_lst, yerr = comb_err, fmt = ".",
    label = r"propagating just $A_{bg}$ fit errors",
    ms = 18, lw = 5,
    )



ax1.set_xlabel(r"T [K]")
ax1.set_ylabel(r"$\alpha_E$")


h, l = ax1.get_legend_handles_labels()
select = [1,0]
ax1.legend([h[i] for i in select], [l[i] for i in select],
               #shadow = True,
               )

#plt.legend(shadow=True)
plt.tight_layout()
plt.grid(True)

format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + "alphaE_from_Abg"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
###################







