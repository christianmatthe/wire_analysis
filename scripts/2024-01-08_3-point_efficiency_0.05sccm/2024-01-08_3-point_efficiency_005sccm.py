import os
import numpy as np
import matplotlib.pyplot as plt

import wire_analysis as wa
from wire_analysis.utils import load_json_dict, load_extractor_dict_json
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient, TC_to_T_Hack, Cv)
from wire_analysis.flow_on_off_cycle_analysis import (
                                        sort_by_z_list)



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

# Plot all the 1sccm runs in one plot
# Need to add run dicts for 0A and 720TC
# Use 0A: 2023-12-22_1sccm_0A_z-scan_jf+hg_wire
# 720TC: 2023-04-21_1sccm_15A_TC_z-scan_jf_wire


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
#     os.makedirs(beamfit.out_dir, exist_ok=True)
#     beamfit.save_json_run_dict()
#     ############### END TODO


#     beamfit.default_plot_data()

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
    return plot_dict

# plot data into shared plot   
#######
# TC_lst = [200, 300, 390, 475]
# indicator_list = ["720TC", "475TC", "390TC", "300TC","200TC",  "0A"]
indicator_list = ["715TC", "0A"]
filename_list = ["005sccm_" + indicator + ".json" 
                 for indicator in indicator_list]
# for TC_lst 44 is a HACK for 0 A ("room temp" =  44 = 297.7K)
#TC_lst = [715, 475, 390, 300, 200, 44] 
TC_lst = [715, 44] 
T_lst = [TC_to_T_Hack(TC) for TC in TC_lst]
pd_dict = {}
# Start plot:
plotname = "multi_run_power"
fig = plt.figure(0, figsize=(8,6.5), dpi =300)
ax1=plt.gca()
x_label = r"$z_{pos}$ [mm]"
for i,filename in enumerate(filename_list):
    print("filename:", filename)
    pd = p_data_plot_dict(filename)
    pd_dict[indicator_list[i]] = pd
    pd_dict[indicator_list[i]]["index"] = i
    pd_dict[indicator_list[i]]["T"] = T_lst[i]
    ### ax1
    ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
                label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))


ax1.set_ylabel(r"power [µW]")
ax1.set_xlabel(x_label)

ax1.grid(True)
ax1.legend(shadow=True, fontsize = 13)
# ax1.tight_layout()


fig.tight_layout()

format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + plotname
            + '.{}'.format(format_im),
            format=format_im, dpi=dpi)
# plt.show()
ax1.cla()
fig.clf()
plt.close()
# ######
# # print(pd_dict)
# # Define prototype "3-point efficiency"
# def three_point_CEB(low_pd, mid_pd, high_pd, ac_H = 1, gamma_H = 1):
#     # align data on z-axis (of most sparese data set)
#     z_lst = []
#     p_excess_lst = []
#     ceb_lst = []
#     for i,z in enumerate(high_pd["z_arr"]):
#         if (z in low_pd["z_arr"]) and (z in mid_pd["z_arr"]):
#             i_low = np.where(low_pd["z_arr"] == z)
#             i_mid = np.where(mid_pd["z_arr"] == z)
#             z_lst.append(z)
#             #print(z_lst)
#             #print(i,i_low,i_mid)

#             # 2 point line fit
#             Ts = [low_pd["T"], mid_pd["T"]]
#             xs = Ts
#             ps = [low_pd["p_arr"][i_low], mid_pd["p_arr"][i_mid]]
#             ys = ps
#             m = (ys[1] - ys[0])/(xs[1] - xs[0])
#             b = ys[0] - xs[0] * m
#             p_fit = lambda T : m * T + b

#             #Calculate excess power at high T
#             p_excess = high_pd["p_arr"][i] - p_fit(high_pd["T"])
#             p_excess_lst.append(p_excess)
#             # Calculate cracking efficiency based on comparison to H2 heating
#             #Assumes:
#             # equal T, equal distribution, pretends wire is at T=low_pd["T"]
#             # LAter try:
#             # 1. BAsed on equal H2 and H1 distribution
#             # 2. Fitting based on separate Background(300K) + H2(1250K) + H1
#             # THis is the basic and bodgy method previously called
#             # CEB (cracking efficinecy bodge)
#             # Calculate p_H2
#             p_H2 = lambda T : p_fit(T) - p_fit(low_pd["T"])
            
#             # Calculate energy delivered per H2
#             # ac_H2 estimated from my own measurements
#             # TODO do clean calculation for future use
#             ac_H2 = 0.45
#             #literature ac_H2 is 0.38 in one paper (at one temp)
#             N_A = 6.02214076e23  # [1/mol]
#             T_low = low_pd["T"]
#             E_H2 = lambda T: ac_H2 * (Cv(T)*T - Cv(T_low)*T_low)/N_A
#             # calculate expected ratio
#             # Energy per particle
#             E_rec = 7.1511 * 10**-19 # Joules
#             # https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.013001
#             # equivalent to 4.4634 eV
#             # Technically thermal and recombination should probably have its 
#             # own accomodation coefficient TODO
#             # Boltzman constant kb in [J/K]
#             kb = 1.38064852 * 10**-23 
#             # Effective energyper H assuming for every 2 H an H2 of equal T 
#             # would also have hit the wire
#             # using ac_H = beta_H the"reaction energy accomodation coefficient"
#             #  as per 
#             # Formally ac_H (alpha_H), beta_H and gamma_H are all separate
#             #  parameters
#             #https://pubs.acs.org/doi/epdf/10.1021/j100801a014
#             beta_H = ac_H
#             E_H = lambda T : (gamma_H * ac_H * (E_rec / 2) 
#                             + ac_H * (3/2 *kb * (T-low_pd["T"]))
#                             - E_H2(T)/2
#                                      ) 
#             # E_rec + differenc ein thermal energy from H2 and H at same T

#             # print(low_pd["T"], high_pd["T"])
#             n_H = p_excess/E_H(high_pd["T"])
#             # Assumes mid T can be scaled to high T
#             n_H2 = p_H2(mid_pd["T"])/E_H2(mid_pd["T"])
#             # number of H divided by number of total H
#             ceb = (n_H) / ((n_H) + 2 * n_H2)
#             ceb_lst.append(ceb)

#     # print(z_lst, p_excess_lst)

#     return z_lst,ceb_lst, p_excess_lst,

# z_lst, ceb_lst, p_excess_lst = three_point_CEB(
#     pd_dict["0A"],pd_dict["390TC"],pd_dict["720TC"])
# z_lst, ceb_lst_ac_H_lit, p_excess_lst = three_point_CEB(
#     pd_dict["0A"],pd_dict["390TC"],pd_dict["720TC"],ac_H=0.65, gamma_H=0.03)
# #ac_H lit from https://doi.org/10.1039/TF9716702711
# # Need to do some work with  this source verifiying it means what I think TODO
# # also referencing https://pubs.acs.org/doi/epdf/10.1021/j100801a014
# # y',the recombination coefficient defined as the fraction of incident 
# # atoms which recombine,

# # # plot p_excess_lst
# # # Start plot:
# # plotname = "excess_power_basic"
# # fig = plt.figure(0, figsize=(8,6.5), dpi =300)
# # ax1=plt.gca()
# # x_label = r"$z_{pos}$ [mm]"
# # # ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
# # #                 label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))
# # ax1.scatter(z_lst, p_excess_lst, #fmt = ".",
# #                 #label = "excess power"
# #                 )


# # ax1.set_ylabel(r"excess power [µW]")
# # ax1.set_xlabel(x_label)

# # ax1.grid(True)
# # ax1.legend(shadow=True, fontsize = 13)
# # # ax1.tight_layout()


# # fig.tight_layout()

# # format_im = 'png' #'pdf' or png
# # dpi = 300
# # plt.savefig(out_dir + plotname
# #             + '.{}'.format(format_im),
# #             format=format_im, dpi=dpi)
# # # plt.show()
# # ax1.cla()
# # fig.clf()
# # plt.close()

# # plot ceb_lst
# # Start plot:
# plotname = "CEB_plot"
# fig = plt.figure(0, figsize=(8,6.5), dpi =300)
# ax1=plt.gca()
# x_label = r"$z_{pos}$ [mm]"
# # ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
# #                 label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))
# ax1.scatter(z_lst, ceb_lst, #fmt = ".",
#                 label = "ac_H = 1" + r",$\gamma = 1$"
#                 )
# ax1.scatter(z_lst, ceb_lst_ac_H_lit, #fmt = ".",
#                 label = "ac_H = 0.65" + r",$\gamma = 0.03$"
#                 )


# ax1.set_ylabel(r"dissociation Efficiency (estimate)")
# ax1.set_xlabel(x_label)

# ax1.grid(True)
# ax1.legend(shadow=True, fontsize = 13)
# # ax1.tight_layout()


# fig.tight_layout()

# format_im = 'png' #'pdf' or png
# dpi = 300
# plt.savefig(out_dir + plotname
#             + '.{}'.format(format_im),
#             format=format_im, dpi=dpi)
# # plt.show()
# ax1.cla()
# fig.clf()
# plt.close()


# # plot ceb_lst
# # Start plot:
# plotname = "CEB_plot_spread"
# fig = plt.figure(0, figsize=(8,6.5), dpi =300)
# ax1=plt.gca()
# x_label = r"$z_{pos}$ [mm]"
# # ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
# #                 label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))
# ax1.scatter(z_lst, ceb_lst, #fmt = ".",
#                 label = "ac_H = 1" + r",$\gamma = 1$"
#                 )
# ax1.scatter(z_lst, ceb_lst_ac_H_lit, #fmt = ".",
#                 label = "ac_H = 0.65" + r",$\gamma = 0.03$"
#                 )
# # maximize CEB: 0.65+-0.2, 0.03+-0.01
# ac_H = 0.45
# gamma_H = 0.01
# z_lst, ceb_lst_ac_H_lit, p_excess_lst = three_point_CEB(
#     pd_dict["0A"],pd_dict["390TC"],pd_dict["720TC"],ac_H=ac_H, gamma_H=gamma_H)
# ax1.scatter(z_lst, ceb_lst_ac_H_lit, #fmt = ".",
#                 label = (f"ac_H = {ac_H:.2f}" + r",$\gamma =$" 
#                          + f"{gamma_H:.2f}")
#                 )
# # minimize CEB: 0.65+-0.2, 0.03+-0.01
# ac_H = 0.85
# gamma_H = 0.05
# z_lst, ceb_lst_ac_H_lit, p_excess_lst = three_point_CEB(
#     pd_dict["0A"],pd_dict["390TC"],pd_dict["720TC"],ac_H=ac_H, gamma_H=gamma_H)
# ax1.scatter(z_lst, ceb_lst_ac_H_lit, #fmt = ".",
#                 label = (f"ac_H = {ac_H:.2f}" + r",$\gamma =$" 
#                          + f"{gamma_H:.2f}")
#                 )


# ax1.set_ylabel(r"dissociation Efficiency (estimate)")
# ax1.set_xlabel(x_label)

# ax1.grid(True)
# ax1.legend(shadow=True, fontsize = 13)
# # ax1.tight_layout()


# fig.tight_layout()

# format_im = 'png' #'pdf' or png
# dpi = 300
# plt.savefig(out_dir + plotname
#             + '.{}'.format(format_im),
#             format=format_im, dpi=dpi)
# # plt.show()
# ax1.cla()
# fig.clf()
# plt.close()











