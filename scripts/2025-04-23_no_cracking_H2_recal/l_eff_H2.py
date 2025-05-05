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


TC_lst = [200, 300, 390, 475]
T_lst = [TC_to_T_Hack(TC) for TC in TC_lst]
print(T_lst)
# The conversion to T is  done by eye based on a plot by Max. For annyone who 
# has nto noticed yet: that is a major  HACK
# T_lst = [750,1000,1250,1500]
leff_lst = []
leff_err_lst = []
for i,TC in enumerate(TC_lst):
    path = work_dir + f"1sccm_{TC}TC.json"
    rd = load_json_dict(path)
    extractor_dict = load_extractor_dict_json(rd["extractor_dict_path"])

    leff = rd["fit_result"]["l_eff"] 
    leff_err = rd["fit_result_errors"]["l_eff"]

    leff_lst.append(leff)
    leff_err_lst.append(leff_err)

# TODO Max power to max power comparison at different  temperatures
############ Plot results with errors
fig = plt.figure(0, figsize=(8,6.5))
ax1= plt.gca()
ax1.errorbar(T_lst, leff_lst, yerr = leff_err_lst, fmt = ".",
    label = r"$l_{\rm eff}$",
    ms = 18, lw = 5,
    )


ax1.set_xlabel(r"T [K]")
ax1.set_ylabel(r"$l_{\rm eff}$")

plt.legend(shadow=True, loc = "upper left")
plt.tight_layout()
plt.grid(True)

format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + "l_eff"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
###################







