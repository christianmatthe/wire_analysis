import wire_analysis as wa
import os
import numpy as np
from wire_analysis.utils import load_json_dict, load_extractor_dict_json
from wire_analysis.accommodation_coefficient import (
    calc_accomodation_coefficient)

# Run through all the analysis dicts in question
work_dir = "./run_dicts/"

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


for filename in os.listdir(work_dir):
    print("filename:", filename)
    path = work_dir + filename
    rd = load_json_dict(path)
    print(rd["fit_result"]["l_eff"])

# TODO Time to implement the accommodation_coeffcient methods
# folowing "accomodation_coeffficient_H2_2023-12-06.ipynb"
# We jsut needed all that, so we could call up  fit results and the 
# extractor_dict by filename :D

# HACK The below is not a clean solution (hiddenin the fucntions it uses)
TC_lst = [200, 300, 390, 475]
# The conversion to T is  done by eye based on a plot by Max. For annyone who 
# has nto noticed yet: that is a major  HACK
T_lst = [750,1000,1250,1500]
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
    )
    print(f"alphaE({T_lst[i]}K):", ac)





