import os
import numpy as np
import matplotlib.pyplot as plt

import wire_analysis as wa
from wire_analysis.utils import load_json_dict, load_extractor_dict_json
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient, TC_to_T_Hack)


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
for i,TC in enumerate(TC_lst):
    path = work_dir + f"1sccm_{TC}TC.json"
    rd = load_json_dict(path)
    extractor_dict = load_extractor_dict_json(rd["extractor_dict_path"])
    # THis is about to be the HACK iest extraction of beam center and out of
    # beam ever
    # p_measured = (0.28 - (-0.08)) * 1e-6 # Maximum power 
    # vs "out of beam" power
    p_measured = (np.max(extractor_dict["p_arr"])
                  - np.min(extractor_dict["p_arr"])) * 1e-6
    # print(p_measured)
    ac = calc_accomodation_coefficient(
        p_measured = p_measured,
        T = T_lst[i],
        l_eff = rd["fit_result"]["l_eff"],
        theta_max = rd["fit_result"]["theta_max"],
        y0 = rd["fit_result"]["y0"],
    )[0]
    ac_err_leff = [ac - calc_accomodation_coefficient(
        p_measured = p_measured,
        T = T_lst[i],
        l_eff = rd["fit_result"]["l_eff"] - rd["fit_result_errors"]["l_eff"],
        theta_max = rd["fit_result"]["theta_max"],
        y0 = rd["fit_result"]["y0"],
    )[0], 
    ac - calc_accomodation_coefficient(
        p_measured = p_measured,
        T = T_lst[i],
        l_eff = rd["fit_result"]["l_eff"] + rd["fit_result_errors"]["l_eff"],
        theta_max = rd["fit_result"]["theta_max"],
        y0 = rd["fit_result"]["y0"],
    )[0]
    ]

    # # HACK
    # ac = calc_accomodation_coefficient(
    #     p_measured = p_measured,
    #     T = T_lst[i],
    #     l_eff = rd["fit_result"]["l_eff"],
    #     theta_max = rd["fit_result"]["theta_max"],
    #     y0 = 40.5,
    # )[0]
    # ac_err_leff = [ac - calc_accomodation_coefficient(
    #     p_measured = p_measured,
    #     T = T_lst[i],
    #     l_eff = rd["fit_result"]["l_eff"] - rd["fit_result_errors"]["l_eff"],
    #     theta_max = rd["fit_result"]["theta_max"],
    #     y0 = 40.5,
    # )[0], 
    # ac - calc_accomodation_coefficient(
    #     p_measured = p_measured,
    #     T = T_lst[i],
    #     l_eff = rd["fit_result"]["l_eff"] + rd["fit_result_errors"]["l_eff"],
    #     theta_max = rd["fit_result"]["theta_max"],
    #     y0 = 40.5,
    # )[0]
    # ]


    print(f"ac_err_leff", ac_err_leff)
    print(f"alphaE({T_lst[i]}K):", ac)
    # Propagate errors in p only
    p_max_err = extractor_dict["p_err_arr"][np.argmax(extractor_dict["p_arr"])]
    p_min_err = extractor_dict["p_err_arr"][np.argmin(extractor_dict["p_arr"])]
    p_err = np.sqrt(p_max_err**2 + p_min_err**2)* 1e-6
    p_rel_err = p_err/p_measured
    ac_err = p_rel_err * ac
    print(f"alphaE({T_lst[i]}K):", f"{ac:2.3f}" , "+-", f"{ac_err:2.3f}")
    ac_lst.append(ac)
    ac_err_lst.append(ac_err)
    ac_err_leff_lst.append(ac_err_leff)
print("ac_lst", ac_lst)
print("ac_err_leff_lst", ac_err_leff_lst)
print("ac_err_leff_lst[0]", ac_err_leff_lst[0])

lst1 = ac_err_leff_lst
lst2 = ac_err_lst
comb_err = np.transpose(np.array([
            [np.abs(lst1[i][0] - lst2[i]),lst1[i][1] + lst2[i]]
             for i in range(len(lst2))]))
print("comb_err", comb_err)
# TODO Max power to max power comparison at different  temperatures
############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()
ax1.errorbar(T_lst, ac_lst, yerr = comb_err, fmt = ".",
    label = r"p and $l_{\rm{eff}}$ errs",
    ms = 18, lw = 5,
    )
ax1.errorbar(T_lst, ac_lst, yerr = ac_err_lst, fmt = ".",
    label = r"p stat errors",
    ms = 18, lw = 5,
    )


ax1.set_xlabel(r"T [K]")
ax1.set_ylabel(r"$\alpha_E$")

plt.legend(shadow=True)
plt.tight_layout()
plt.grid(True)

format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + "alphaE_perr_l_eff_err"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
###################







