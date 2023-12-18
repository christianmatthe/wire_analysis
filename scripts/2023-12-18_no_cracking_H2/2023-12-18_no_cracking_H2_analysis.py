import wire_analysis as wa
import os
from wire_analysis.utils import load_json_dict

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