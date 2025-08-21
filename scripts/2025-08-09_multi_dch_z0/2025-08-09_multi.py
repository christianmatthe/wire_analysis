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

# # extract l_effs
# l_eff_1T_lst = [load_json_dict(path)["fit_result"]["l_eff"]
#                  for path in file_list_1T ]
# l_eff_1T_err_lst = [load_json_dict(path)["fit_result_errors"]["l_eff"]
#                  for path in file_list_1T ]
# print("l_effs_1T", l_eff_1T_lst )


d_ch_1T_lst = [load_json_dict(path)["fit_result"]["d_ch"]
                 for path in file_list_1T ]
d_ch_1T_err_lst = [load_json_dict(path)["fit_result_errors"]["d_ch"]
                 for path in file_list_1T ]
print("d_chs_1T", d_ch_1T_lst )
z0_1T_lst = [load_json_dict(path)["fit_result"]["z0"]
                 for path in file_list_1T ]
z0_1T_err_lst = [load_json_dict(path)["fit_result_errors"]["z0"]
                 for path in file_list_1T ]

#3-T-P:
flow_lst_3T = [#10,
               1,
                0.2, 0.05,
                #0.01, 0.002
               ]
flow_lst_3T_str = [#"10",
                "1",
                 "02", "005",
               #"001", "0002"
               ]
flow_lst_3T_TC_str = [#"720",
               "720", 
               "715", "715",
               #"715", "715"
               ]

#HACK on 1sccm not guaranteed to be 3 point
file_list_3T = [(r"C:\Users\Christian\Documents\StudiumPhD\python\\"
                + r"wire_analysis\scripts\\"
                + r"2025-05-10_6-point_method_1sccm_dch\run_dicts\\"
                + r"1sccm_720TC_penumbra.json"),
    
    (r"C:\Users\Christian\Documents\StudiumPhD\python\\"
                + r"wire_analysis\scripts\\"
                + r"2025-05-10_0.2sccm_cut_data_dch\output\\"
                + r"02sccm_715TC_penumbra\02sccm_715TC_penumbra.json"),
                (r"C:\Users\Christian\Documents\StudiumPhD\python\\"
                + r"wire_analysis\scripts\\"
                + r"2025-05-10_0.05sccm_Cv _cut_data_dch\output\\"
                + r"005sccm_715TC_penumbra\005sccm_715TC_penumbra.json"),

                ]

d_ch_3T_lst = [load_json_dict(path)["fit_result"]["d_ch"]
                 for path in file_list_3T ]
d_ch_3T_err_lst = [load_json_dict(path)["fit_result_errors"]["d_ch"]
                 for path in file_list_3T ]
print("d_chs_3T", d_ch_3T_lst )

z0_3T_lst = [load_json_dict(path)["fit_result"]["z0"]
                 for path in file_list_3T ]
z0_3T_err_lst = [load_json_dict(path)["fit_result_errors"]["z0"]
                 for path in file_list_3T ]

###### H2 only fits:
flow_lst_H2 = [1,1,1,1, 0.2
               ]
flow_lst_H2_str = [#"10",
                "1","1", "1", "1", "02"
               ]
out_dir_suffix = ["","","","", "",
                  "_02sccm"]
flow_lst_H2_TC_str = [#"720",
               "475", 
               "390", "300", "200",
               #"715", "715",
               "390"
               ]
file_list_H2 = [(r"C:\Users\Christian\Documents\StudiumPhD\python\\"
                 + r"wire_analysis\scripts\\"
                + r"2025-04-23_no_cracking_H2_recal\\"
                + f"output{out_dir_suffix[i]}\\"
                + f"{s}sccm_{flow_lst_H2_TC_str[i]}TC\\"
                + f"{s}sccm_{flow_lst_H2_TC_str[i]}TC.json" ) 
                for i,s in enumerate(flow_lst_H2_str)
                ]

d_ch_H2_lst = [load_json_dict(path)["fit_result"]["d_ch"]
                 for path in file_list_H2 ]
d_ch_H2_err_lst = [load_json_dict(path)["fit_result_errors"]["d_ch"]
                 for path in file_list_H2 ]
print("d_chs_H2", d_ch_H2_lst )

z0_H2_lst = [load_json_dict(path)["fit_result"]["z0"]
                 for path in file_list_H2 ]
z0_H2_err_lst = [load_json_dict(path)["fit_result_errors"]["z0"]
                 for path in file_list_H2 ]

#############################################################
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

include_1T = True
###########
fig = plt.figure(0, figsize=(8,6))
ax1=plt.gca()

ax1.errorbar(flow_lst_3T, d_ch_3T_lst,
             yerr = d_ch_3T_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C0",
         label = r"2200K, 3-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

ax1.errorbar(flow_lst_H2, d_ch_H2_lst,
             yerr = d_ch_H2_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C1",
         label = r"H2, multi temp",
         capthick=2,elinewidth=2,capsize=3
         )

# 1T
if include_1T:
    ax1.errorbar(flow_lst_1T, d_ch_1T_lst,
             yerr = d_ch_1T_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C3",
         label = r"2200K, 1-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

ax1.set_xlabel(r"Flow [sccm]")
ax1.set_ylabel(r"$d_{ch} [mm]$")

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
    plt.savefig(out_dir + "d_ch"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi , bbox_inches="tight")
else:
    plt.savefig(out_dir + "d_ch_base"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi, bbox_inches="tight")
# plt.show()
ax1.cla()
fig.clf()
plt.close()


#############
#z0 plot:
###########
fig = plt.figure(0, figsize=(8,6))
ax1=plt.gca()

ax1.errorbar(flow_lst_3T, z0_3T_lst,
             yerr = z0_3T_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C0",
         label = r"2200K, 3-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

ax1.errorbar(flow_lst_H2, z0_H2_lst,
             yerr = z0_H2_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C1",
         label = r"H2, multi temp",
         capthick=2,elinewidth=2,capsize=3
         )

# 1T
if include_1T:
    ax1.errorbar(flow_lst_1T, z0_1T_lst,
             yerr = z0_1T_err_lst,
         marker = ".", ms = 15,
         ls = "",
         #lw = 9, 
         color = "C3",
         label = r"2200K, 1-T-P",
         capthick=2,elinewidth=2,capsize=3
         )

ax1.set_xlabel(r"Flow [sccm]")
ax1.set_ylabel(r"$z_0$ [mm]")

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
    plt.savefig(out_dir + "z0"
                + '.{}'.format(format_im),
                format=format_im, dpi=dpi , bbox_inches="tight")
else:
    plt.savefig(out_dir + "z0_base"
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi, bbox_inches="tight")
# plt.show()
ax1.cla()
fig.clf()
plt.close()




