# Based on 2024-01-02_3-point-method
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


for filename in os.listdir(work_dir):
    print("filename:", filename)
    beamfit = wa.Beamfit(
        run_dict_path = work_dir + filename)
    # HACK to change out dir without editing the files
    # TODO If you do this implement feature to retroactively change out_dir in
    # run_dict
    ############### TODO Implement as base function
    name, file_extension = os.path.splitext(filename)
    beamfit.run_dict["out_dir_base"] = os.path.abspath("./output/") + os.sep 
    beamfit.out_dir = beamfit.run_dict["out_dir_base"] + name + os.sep
    os.makedirs(beamfit.out_dir, exist_ok=True)
    beamfit.save_json_run_dict()
    ############### END TODO


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
indicator_list = ["720TC", "475TC", "390TC", "300TC","200TC",  "0A"]
filename_list = ["1sccm_" + indicator + ".json" 
                 for indicator in indicator_list]
# for TC_lst 44 is a HACK for 0 A ("room temp" =  44 = 297.7K)
TC_lst = [720, 475, 390, 300, 200, 44] 
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
    # ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
    #             label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K"))
    ax1.errorbar(pd["z_arr"], pd["p_arr"],yerr = pd["p_err_arr"], fmt = ".",
        label = (r"$P_{\rm{meas}}(\rm{T}$ $\approx$ "+f"{T_lst[i]:.0f}K)"), 
        markersize =10)


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

######
# print(pd_dict)
def multi_point_CEB(pd_dict,
                    H2_indicators = ["475TC", "390TC", "300TC","200TC",  "0A"],
                    H_indicators = ["720TC"],
                    ac_H = 1, gamma_H = 1):
    # align data on z-axis (of most sparese data set)
    z_lst = []
    p_excess_lst = []
    p_excess_err_lst = []
    ceb_lst = []
    ceb_err_lst = []
    high_pd = pd_dict[H_indicators[0]]
    #Ts H2
    Ts = [pd_dict[key]["T"] for key in H2_indicators]

    for i,z in enumerate(high_pd["z_arr"]):
        checksum = 0
        for key in H2_indicators:
            if z in pd_dict[key]["z_arr"]:
                checksum += 1
        if checksum == len(H2_indicators):
            i_dict = {key:np.where(pd_dict[key]["z_arr"] == z) 
                      for key in H2_indicators}
            # i_low = np.where(low_pd["z_arr"] == z)
            # i_mid = np.where(mid_pd["z_arr"] == z)
            z_lst.append(z)
            #print(z_lst)
            #print(i,i_low,i_mid)
            #print(i_dict)

            # multi point line fit

            xs = Ts
            ps = np.array([pd_dict[key]["p_arr"][i_dict[key]] 
                    for key in H2_indicators]).flatten()
            ps_err = np.array([pd_dict[key]["p_err_arr"][i_dict[key]] 
                    for key in H2_indicators]).flatten()
            ys = ps
            # print("Ts", Ts, np.shape(Ts))
            # print( "ps", ps, np.shape(ps))
            # print( "ps_err", ps_err,  np.shape(ps_err))
            #Replace with real fit
            fit_func = lambda T,a,b : a*Cv(T)*T + b

            popt, pcov = curve_fit(fit_func, Ts, ps
                                ,sigma = ps_err
                                , absolute_sigma = True
                                )
            p_fit = lambda T: popt[0]*Cv(T)*T + popt[1]

            a_err = np.sqrt(pcov[0,0])
            b_err = np.sqrt(pcov[1,1])
            # print("a, aerr", (popt[0],a_err))
            # print("b, berr", (popt[1],b_err))

            p_fit_err =  lambda T : np.sqrt((a_err *Cv(T) * T)**2 + b_err**2)

            #Calculate excess power at high T
            p_excess = high_pd["p_arr"][i] - p_fit(high_pd["T"])
            p_excess_lst.append(p_excess)
            p_excess_err = np.sqrt(high_pd["p_err_arr"][i]**2
                            + (p_fit_err(high_pd["T"]))**2
                            )
            p_excess_err_lst.append(p_excess_err)
            # Calculate cracking efficiency based on comparison to H2 heating
            #Assumes:
            # equal T, equal distribution, pretends wire is at T=low_pd["T"]
            # LAter try:
            # 1. BAsed on equal H2 and H1 distribution
            # 2. Fitting based on separate Background(300K) + H2(1250K) + H1
            # THis is the basic and bodgy method previously called
            # CEB (cracking efficinecy bodge)
            # # Calculate p_H2 # TODO Fix
            # p_H2 = lambda T : p_fit(T) - p_fit(low_pd["T"])
            # p_H2_err = lambda T : np.sqrt(p_fit_err(T)**2 
            #                               + p_fit_err(low_pd["T"])**2
            #                             )
            
            # # Calculate energy delivered per H2
            # # ac_H2 estimated from my own measurements
            # # TODO do clean calculation for future use
            # ac_H2 = 0.45
            # #literature ac_H2 is 0.38 in one paper (at one temp)
            # N_A = 6.02214076e23  # [1/mol]
            # T_low = low_pd["T"]
            # E_H2 = lambda T: ac_H2 * (Cv(T)*T - Cv(T_low)*T_low)/N_A
            # # calculate expected ratio
            # # Energy per particle
            # E_rec = 7.1511 * 10**-19 # Joules
            # # https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.013001
            # # equivalent to 4.4634 eV
            # # Technically thermal and recombination should probably have its 
            # # own accomodation coefficient TODO
            # # Boltzman constant kb in [J/K]
            # kb = 1.38064852 * 10**-23 
            # # Effective energyper H assuming for every 2 H an H2 of equal T 
            # # would also have hit the wire
            # # using ac_H = beta_H the"reaction energy accomodation coefficient"
            # #  as per 
            # # Formally ac_H (alpha_H), beta_H and gamma_H are all separate
            # #  parameters
            # #https://pubs.acs.org/doi/epdf/10.1021/j100801a014
            # beta_H = ac_H
            # E_H = lambda T : (gamma_H * ac_H * (E_rec / 2) 
            #                 + ac_H * (3/2 *kb * (T-low_pd["T"]))
            #                 - E_H2(T)/2
            #                          ) 
            # # E_rec + differenc ein thermal energy from H2 and H at same T

            # # print(low_pd["T"], high_pd["T"])
            # n_H = p_excess/E_H(high_pd["T"])
            # n_H_err = p_excess_err/E_H(high_pd["T"])
            # # Assumes mid T can be scaled to high T
            # n_H2 = p_H2(mid_pd["T"])/E_H2(mid_pd["T"])
            # n_H2_err = p_H2_err(mid_pd["T"])/E_H2(mid_pd["T"])
            # # number of H divided by number of total H
            # ceb = (n_H) / ((n_H) + 2 * n_H2)
            # ceb_lst.append(ceb)
            # ceb_err = np.sqrt(
            #             (n_H_err * (2 * n_H2)/(((n_H) + 2 * n_H2))**2)**2
            #             + (n_H2_err * (2 * n_H)/(((n_H) + 2 * n_H2))**2)**2
            #             )
            # ceb_err_lst.append(ceb_err[0])

    # print(z_lst, p_excess_lst)

    return z_lst,ceb_lst,ceb_err_lst, p_excess_lst,p_excess_err_lst,


def background_subtract_p(low_pd, high_pd):
    #very basic background subtractionn without temperature axis alignment
    # align data on z-axis (of most sparese data set)
    z_lst = []
    p_excess_lst = []
    p_excess_err_lst = []
    for i,z in enumerate(high_pd["z_arr"]):
        if (z in low_pd["z_arr"]):
            i_low = np.where(low_pd["z_arr"] == z)
            z_lst.append(z)

            #Calculate excess power at high T
            p_excess = high_pd["p_arr"][i] - low_pd["p_arr"][i_low]
            p_excess_lst.append(p_excess)
            p_excess_err = np.sqrt(high_pd["p_err_arr"][i]**2
                            + low_pd["p_err_arr"][i_low]**2
                            )
            p_excess_err_lst.append(p_excess_err)
    return z_lst, p_excess_lst,p_excess_err_lst

H2_indicators = ["475TC", 
                             "390TC", 
                             "300TC",
                              "200TC",  
                            "0A"
                            ] 

#"Subsequences of the iterable from shortest to longest. with minimum length 2"
s = H2_indicators
indicator_powerset = chain.from_iterable(combinations(s, r) 
                    for r in range(2, len(s)+1))
set_indicator_powerset = set(indicator_powerset)


# #### START Comment  out here
# indicator_powerset = chain.from_iterable(combinations(s, r) 
#                     for r in range(2, len(s)+1))
# fit_res_dict = {}
# for i, indicators in tqdm(enumerate(indicator_powerset)):
#     print(indicators)
#     z_lst, ceb_lst, ceb_err_lst, p_excess_lst,p_excess_err_lst = multi_point_CEB(
#         pd_dict, H2_indicators = indicators)
#     ###################### Previous is reuse of old code to get p_excess

#     ####################################
#     # Broken filename = "1sccm_720TC_penumbra_"+ "-".join(indicators) + ".json"
#     filename = "1sccm_720TC_penumbra.json"
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

#     # print("p_excess_lst", p_excess_lst)
#     # print("p_arr=np.array(p_excess_lst)", np.array(p_excess_lst).flatten())

#     print("custom penumbra fit wait ~??min")
#     (popt, pcov) = beamfit.custom_fit(z_arr=np.asarray(z_lst).flatten(),
#                             p_arr=np.array(p_excess_lst).flatten()
#                             , p_err_arr = np.array(p_excess_err_lst).flatten()
#                     ,plotname="custom_fit_plot"+ "_"+ "-".join(indicators)
#                     )
#     fit_res_dict[str(indicators)] = {}
#     fit_res_dict[str(indicators)]["popt"] = popt.tolist()
#     fit_res_dict[str(indicators)]["pcov"] = pcov.tolist()


# save_json_dict("all_fit_params.json", fit_res_dict)
# #### END Comment  out here

indicator_powerset = set_indicator_powerset
#Just load:
fit_res_dict = load_json_dict("all_fit_params.json")
# print(fit_res_dict)
# print(fit_res_dict["('475TC', '390TC')"]["pcov"][0][0])
# print([str(key) for key in indicator_powerset] )
l_eff_list = [fit_res_dict[str(key)]["popt"][0] for key in indicator_powerset]
print([fit_res_dict[str(key)]["pcov"] for key in indicator_powerset])
l_eff_errs = [np.sqrt(fit_res_dict[str(key)]["pcov"][0][0])
               for key in indicator_powerset]
# print("l_eff_list:", l_eff_list)
print("l_eff:", f"{np.mean(l_eff_list):.3f}", " +- "
      , f"{np.std(l_eff_list):.3f}" )

# print("l_eff_errs", l_eff_errs)
values = np.array(l_eff_list)
errors = np.array(l_eff_errs)
n = len(values)
weights = 1/ (errors**2)
# Unweighted
mean = np.sum(values) / n
std = np.sqrt(np.sum((values - mean)**2)/(n-1))
#Weighted
w_mean = np.sum(weights * values) / np.sum(weights)
# weighted standat deviation
w_std = np.sqrt(np.sum(weights * (values - w_mean)**2)
                / (((n-1)/n) * np.sum(weights)  )
)
# the std of mean will liekly be a gross underestimate because it will
# assume fit errors to be  genuine and independent
std_of_mean = np.sqrt(1/(np.sum(weights)))
print("mean= ", mean)
print("std= ", std )
print("w_mean = ", w_mean)
print("w_std = ", w_std)
print("std_of_mean =", std_of_mean)

print(sorted(l_eff_list)[0:-1])
print("cut span mean= ", np.sum(sorted(l_eff_list)[0:-1]) / n)
print(sorted(l_eff_list))


# # With just length 2
# indicator_powerset = chain.from_iterable(combinations(s, r) 
#                     for r in range(2, 2+1))
# set_indicator_powerset = set(indicator_powerset)
# indicator_powerset = set_indicator_powerset

# l_eff_list = [fit_res_dict[str(key)]["popt"][0] for key in indicator_powerset]
# print([fit_res_dict[str(key)]["pcov"] for key in indicator_powerset])
# l_eff_errs = [np.sqrt(fit_res_dict[str(key)]["pcov"][0][0])
#                for key in indicator_powerset]
# # print("l_eff_list:", l_eff_list)
# print("l_eff:", f"{np.mean(l_eff_list):.3f}", " +- "
#       , f"{np.std(l_eff_list):.3f}" )

# # print("l_eff_errs", l_eff_errs)
# values = np.array(l_eff_list)
# errors = np.array(l_eff_errs)
# n = len(values)
# weights = 1/ (errors**2)
# # Unweighted
# mean = np.sum(values) / n
# std = np.sqrt(np.sum((values - mean)**2)/(n-1))
# #Weighted
# w_mean = np.sum(weights * values) / np.sum(weights)
# # weighted standat deviation
# w_std = np.sqrt(np.sum(weights * (values - w_mean)**2)
#                 / (((n-1)/n) * np.sum(weights)  )
# )
# # the std of mean will liekly be a gross underestimate because it will
# # assume fit errors to be  genuine and independent
# std_of_mean = np.sqrt(1/(np.sum(weights)))
# print("mean= ", mean)
# print("std= ", std )
# print("w_mean = ", w_mean)
# print("w_std = ", w_std)
# print("std_of_mean =", std_of_mean)

# print(sorted(l_eff_list))

dat = values
# nBins = 20 
nBins = 12 
title='Histogram'
hist_range = []

fig = plt.figure(0, figsize=(8,6.5), dpi =300)
if hist_range == []:
    hist_range = np.array([np.min(dat), np.max(dat)])
bin_width =  (hist_range[0]-hist_range[-1])/nBins
#nBins=round(abs(hist_range[-1]-hist_range[0])/bin_width)+1

binList=np.linspace(hist_range[0],
                    hist_range[-1],num=nBins + 1)

#print(binList)
entries,bin_edges=np.histogram(dat,bins=binList,range=hist_range)

bin_centers = (bin_edges[:-1]+bin_edges[1:]) / 2
#bin_width=(abs(bin_edges[-1]-bin_edges[0]))/(nBins)

a=np.empty([1,3])
for i in range(0,nBins):
    a=np.append(a,[[i,bin_centers[i],entries[i]]],axis=0)
    

l_eff = fit_res_dict["('475TC', '390TC', '300TC', '200TC', '0A')"]["popt"][0]
err = np.sqrt(
    fit_res_dict["('475TC', '390TC', '300TC', '200TC', '0A')"]["pcov"][0][0])
plt.axvline(l_eff, 
        #     label=r"All Sets " + r"$l_{\rm eff} = $"
        # + f"{l_eff:.2f}", 
        color = "C0", alpha = 1, ls = "--", lw  = 3)   

plt.axvspan(l_eff - err, l_eff + err, label=r"All Sets, " + r"$l_{\rm eff} = $"
        + f"{l_eff:.2f}" + r"$\pm$" + f"{err:.2f}", 
        color = "C0", alpha = 0.5, ls = "--")   

# l_eff = fit_res_dict["('390TC', '0A')"]["popt"][0]
# err = np.sqrt(
#     fit_res_dict["('390TC', '0A')"]["pcov"][0][0])
l_eff = 4.20
err = 0.22

plt.axvline(l_eff, 
        #     label=r"720TC-390TC-0A " + r"$l_{\rm eff} = $"
        # + f"{l_eff:.2f}", 
        color = "C2", alpha = 1, ls = "--", lw = 3)    
plt.axvspan(l_eff - err, l_eff + err, 
            label=r"2211K-1277K-298K,  " + r"$l_{\rm eff} = $"
        + f"{l_eff:.2f}" + r"$\pm$" + f"{err:.2f}", 
        color = "C2", alpha = 0.5, ls = "--")

#plt.title(title)
plt.ylabel('Count')
plt.xlabel(r'$l_{\rm eff}$')
plt.bar(bin_centers, entries, align='center',width=bin_width,edgecolor='k',
        color = "C1", label  = "Binned Data")
plt.legend(fontsize = 13, shadow = True)
plt.tight_layout()

format_im = 'png' #'pdf' or png
dpi = 300
plt.savefig(out_dir + title
        + '.{}'.format(format_im),
        format=format_im, dpi=dpi)
plt.close()
# return bin_centers, entries, bin_width
