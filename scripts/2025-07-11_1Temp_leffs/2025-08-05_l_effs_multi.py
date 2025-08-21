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

# extract l_effs
l_eff_1T_lst = [load_json_dict(path)["fit_result"]["l_eff"]
                 for path in file_list_1T ]
l_eff_1T_err_lst = [load_json_dict(path)["fit_result_errors"]["l_eff"]
                 for path in file_list_1T ]
print("l_effs_1T", l_eff_1T_lst )

# ac_lst = []
# ac_err_lst = []
# ac_err_leff_lst = []
# for i,TC in enumerate(flow_lst):
#     path = file_list[i]
#     rd = load_json_dict(path)


#     ad = a_diss_from_A(
#         A = rd["fit_result"]["A"],
#         flow = flow_lst[i] * 4.478 * 10**17 #sccm
#         )
#     ad_err = [ad - a_diss_from_A(
#                 A = rd["fit_result"]["A"] - rd["fit_result_errors"]["A"],
#                 flow = flow_lst[i] * 4.478 * 10**17 #sccm
#                 ), 
#               ad - a_diss_from_A(
#                 A = rd["fit_result"]["A"] + rd["fit_result_errors"]["A"],
#                 flow = flow_lst[i] * 4.478 * 10**17 #sccm
#                 )
#             ]

#     ac_lst.append(ad)
#     ac_err_lst.append(ad_err)
# #     ac_err_leff_lst.append(ac_err_leff)
# print("ad_lst", ac_lst)
# print("ad_errs", ac_err_lst)
# print("errs_abs",  np.abs(np.array(ac_err_lst)).T)

# Based on C:\Users\Christian\Documents\StudiumPhD\python\Keysight-DMM-34461A\
# analysis\Beam_profile_multi_flow_2025-03-02.ipynb

font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 20
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.sans-serif": "Computer Modern",
    "font.weight":"bold"
})

# 2025-03-02 sim eta, penumbra model, 3 point method, cut data
flow_list_3point = [0.05,0.2, 1]
l_eff_list_3point = [7.55,6.75,4.20]
# l_eff_err_list_3point = [0.71,0.31,0.22]

l_eff_err_list_3point = [0.71,0.31,0.22]

#Tschersich 2000
# Data from Table 1: J. Appl. Phys., Vol. 87, No. 5, 1 March 2000
l_eff_tscher = [3.8, 7, 11.4, 13]
flow_tscher = [1.89e17, 4.15e16, 3.68e15, 8.89e14]
particles_per_second_sccm = 4.48e17
flow_sccm_tscher = [flow / particles_per_second_sccm for flow in flow_tscher]

include_1T = True
###########
fig = plt.figure(0, figsize=(8,6))
ax1=plt.gca()

ax1.errorbar(flow_list_3point, l_eff_list_3point,
             yerr = l_eff_err_list_3point,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C0",
         label = "This Work, 3-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

#HACK HArdcoded span for multi T
ax1.errorbar(flow_list_3point[2], l_eff_list_3point[2],
             yerr = [[l_eff_list_3point[2]- 3.48], 
                     [np.abs(l_eff_list_3point[2]- 4.91)]],
         marker = "", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C1",
         label = "syst. span, 3-T-P ",
         capthick=2,elinewidth=2,capsize=3
         )


ax1.plot(flow_sccm_tscher, l_eff_tscher, ".", color = "C2", markersize = 30, 
         label = "Tschersich (2000)", markeredgecolor = "k")

# 1T
if include_1T:
    ax1.errorbar(flow_lst_1T, l_eff_1T_lst,
             yerr = l_eff_1T_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C3",
         label = r"This Work, 1-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

ax1.set_xlabel(r"Flow [sccm]")
ax1.set_ylabel(r"$l_{\rm eff}$")

ax1.set_xscale('log')

#Introduce seconndary axis
def sccm_to_atoms(sccm):
    atoms_per_sccm = 4.477962e17
    return sccm *atoms_per_sccm

def atoms_to_sccm(atoms):
    atoms_per_sccm = 4.477962e17
    return atoms / atoms_per_sccm

secax2 = ax1.secondary_xaxis(1, 
                            functions=(sccm_to_atoms, atoms_to_sccm))
secax2.set_xlabel('Flow [Molecules per s]')

plt.legend(shadow=True)
plt.tight_layout()
#plt.grid(True)

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
ax1.yaxis.set_minor_locator(AutoMinorLocator())
ax1.grid(which = "major", color = "0.5")
ax1.grid(which="minor", color = "0.9")

format_im = 'png' #'pdf' or png
dpi = 300
if include_1T:
    plt.savefig(out_dir + "l_effs"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi , bbox_inches="tight")
else:
    plt.savefig(out_dir + "l_effs_base"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi, bbox_inches="tight")
# plt.show()
ax1.cla()
fig.clf()
plt.close()





