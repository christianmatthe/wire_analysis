import wire_analysis as wa
import os
beamfit = wa.Beamfit(
    run_dict_path = 
        ("C:/Users/Christian/Documents/StudiumPhD/python/wire_analysis"
                            + os.sep +"\\scripts\\output\\test.json"))
beamfit.default_fit()