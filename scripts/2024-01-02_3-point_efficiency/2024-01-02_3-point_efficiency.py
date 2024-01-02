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


    beamfit.default_plot_data()
#


# for filename in os.listdir(work_dir):
#     print("filename:", filename)
#     path = work_dir + filename
#     rd = load_json_dict(path)
#     print(rd["fit_result"]["l_eff"])

# TODO Time to implement the accommodation_coeffcient methods
# folowing "accomodation_coeffficient_H2_2023-12-06.ipynb"
# We jsut needed all that, so we could call up  fit results and the 
# extractor_dict by filename :D

# # HACK The below is not a clean solution (hiddenin the fucntions it uses)
# TC_lst = [200, 300, 390, 475]
# T_lst = [TC_to_T_Hack(TC) for TC in TC_lst]
# print(T_lst)
# # The conversion to T is  done by eye based on a plot by Max. For annyone who 
# # has nto noticed yet: that is a major  HACK
# # T_lst = [750,1000,1250,1500]








