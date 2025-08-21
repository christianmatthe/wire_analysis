# Based on 2024-01-02_3-point-method
import os
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import curve_fit

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
work_dir = "./run_dicts_02sccm/"
out_dir = "./output_02sccm/"


for filename in os.listdir(work_dir):
    print("filename:", filename)
    beamfit = wa.Beamfit(
        run_dict_path = work_dir + filename)
    # HACK to change out dir without editing the files
    # TODO If you do this implement feature to retroactively change out_dir in
    # run_dict
    ############### TODO Implement as base function
    name, file_extension = os.path.splitext(filename)
    beamfit.run_dict["out_dir_base"] = os.path.abspath(out_dir) + os.sep 
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
indicator_list = [ "390TC",  "0A"]
filename_list = ["02sccm_" + indicator + ".json" 
                 for indicator in indicator_list]
# for TC_lst 44 is a HACK for 0 A ("room temp" =  44 = 297.7K)
TC_lst = [390,  44] 
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

######
# print(pd_dict)
def multi_point_CEB(pd_dict,
                    H2_indicators = ["390TC",   "0A"],
                    H_indicators = ["715TC"],
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

            p_fit_err =  lambda T : np.sqrt((a_err *Cv(T) * T)**2 + b_err**2)

            #Calculate excess power at high T
            p_excess = high_pd["p_arr"][i] - p_fit(high_pd["T"])
            p_excess_lst.append(p_excess)
            p_excess_err = np.sqrt(high_pd["p_err_arr"][i]**2
                            + (p_fit_err(high_pd["T"]))**2
                            )
            p_excess_err_lst.append(p_excess_err)

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
                            )[0]
            p_excess_err_lst.append(p_excess_err)
            
    return z_lst, p_excess_lst,p_excess_err_lst

# H2_indicators = ["475TC", 
#                              "390TC", 
#                              "300TC",
#                               "200TC"
#                             ] 
H2_indicators = ["390TC"
                ]

#Run for all intermediate Temperatures
for i,indicator in enumerate(H2_indicators):
    z_lst, p_excess_lst,p_excess_err_lst = (
        background_subtract_p(
        pd_dict["0A"], pd_dict[indicator]
        ))
    # print("z_lst, p_excess_lst,p_excess_err_lst")
    # print(z_lst, p_excess_lst,p_excess_err_lst)


    #HACK Chop out the "probleatic" section around 20 deg
    # mask = np.ones(len(z_lst), dtype=bool)
    # mask[[len(z_lst) - 6,len(z_lst) - 5, len(z_lst) - 4]] = False
    # print("mask", mask)
    # z_lst, p_excess_lst,p_excess_err_lst = (z_lst[mask],
    #                                          p_excess_lst[mask],
    #                                          p_excess_err_lst[mask])
    
    # # del_arr = [len(z_lst) - 6,len(z_lst) - 5, len(z_lst) - 4]
    # del_arr = [len(z_lst) - 7,len(z_lst) - 6,len(z_lst) - 5, len(z_lst) - 4,
    #            len(z_lst) - 3]
    # z_lst = np.delete(z_lst, del_arr,axis=0)
    # p_excess_lst= np.delete(p_excess_lst, del_arr,axis=0)
    # p_excess_err_lst = np.delete(p_excess_err_lst, del_arr,axis=0)

    # plot p_excess_lst
    # Start plot:
    plotname = "excess_power_basic"
    fig = plt.figure(0, figsize=(8,6.5), dpi =300)
    ax1=plt.gca()
    x_label = r"$z_{pos}$ [mm]"
    ax1.errorbar(z_lst, p_excess_lst,yerr = p_excess_err_lst, fmt = ".",
                    #label = (f"{indicator_list[i]} ~ {T_lst[i]:.1f}K")
                    )
    # ax1.scatter(z_lst, p_excess_lst, #fmt = ".",
    #                 #label = "excess power"
    #                 )


    ax1.set_ylabel(r"excess power [µW]")
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
    ###################### Previous is reuse of old code to get p_excess

    ####################################
    filename = f"02sccm_{indicator}.json"
    beamfit = wa.Beamfit(
        run_dict_path = work_dir + filename)
    # HACK to change out dir without editing the files
    # TODO If you do this implement feature to retroactively change out_dir in
    # run_dict
    ############### TODO Implement as base function
    name, file_extension = os.path.splitext(filename)
    beamfit.run_dict["out_dir_base"] = os.path.abspath(out_dir) + os.sep 
    beamfit.out_dir = beamfit.run_dict["out_dir_base"] + name + os.sep
    os.makedirs(beamfit.out_dir, exist_ok=True)
    beamfit.save_json_run_dict()
    ############### END TODO

    # print("dch and y0 fit wait ~??min")
    # beamfit.fit_d_ch(z_arr=np.asarray(z_lst).flatten(),
    #                         p_arr=np.array(p_excess_lst).flatten()
    #                         , p_err_arr = np.array(p_excess_err_lst).flatten()
    #                 ,plotname="y0_fit_plot",
    #                 fit_y0=True
    #                 )
    
    print("dch fit wait ~??min")
    beamfit.fit_d_ch(z_arr=np.asarray(z_lst).flatten(),
                            p_arr=np.array(p_excess_lst).flatten()
                            , p_err_arr = np.array(p_excess_err_lst).flatten()
                    ,plotname="dch_fit_plot"
                    )



    # print("custom penumbra fit wait ~??min")
    # beamfit.custom_fit(z_arr=np.asarray(z_lst).flatten(),
    #                         p_arr=np.array(p_excess_lst).flatten()
    #                         , p_err_arr = np.array(p_excess_err_lst).flatten()
    #                 ,plotname="custom_fit_plot"
    #                 )










