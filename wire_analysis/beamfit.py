import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.interpolate import interp1d

import os

#######Imports from  other submodules
# Put  these into a  utils.py?
from .flow_on_off_cycle_analysis import (load_dict, make_result_dict, 
                                        sort_by_z_list)

########

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)
mpl.rcParams['figure.dpi'] = 400


# make a class for fitting l_eff Tschersich model to beam shape
# Started 2023-Dec-11
# based on "HABS_beam_profile_fitting_4.1.1.1_cos3_low_temp_1sccm.ipynb"
# requirements: (non-exhaustive)
#   - allow for modular eta_wire
#   - Properly save outpput for further use
#   - Take output of flow_on_off_cycle_analysis_2.py
#

class Beamfit():
    """
    "The beamfit class contains functions and parameters necessary to
    fit beam profiles with the Tschersich model. 
    (K. G. Tschersich; V. von Bonin (1998))
    "

    Parameters
    ----------
    y0_default :  `float`
        default distance from wire to source [mm]
    
    """
    def __init__(self,
                y0_default = 35.17, # mm (estimated from CAD)
                ):
        #define parameters
        self.y0_default =  y0_default # mm (estimated from CAD)



        #define constants
        degree = np.pi/180 # convert form rad to degree
        
        # run init fucntions
        self.init_eta_wire_sim()
            #Allow for later adjustment of default eta 
        self.eta_wire_default = self.eta_wire_sim
        return


    # Auxiliary functions required for profile shape according to
    # Tschcersich  https://aip.scitation.org/doi/pdf/10.1063/1.368619
    #  (K. G. Tschersich; V. von Bonin (1998))

    # remember to input values in radians
    def beta(self, theta, l_eff):

        output = np.zeros_like(theta) 
        cond = np.abs(l_eff * np.tan(theta)) < 1
        output[cond] =  np.arccos(l_eff * np.tan(theta[cond]))
        return output


    def U(self, theta, l_eff):
        # Move conditional to beta only
        return ((2*self.beta(theta, l_eff)-np.sin(2*self.beta(theta, l_eff)))
                /np.pi)

    def V(self, theta, l_eff):
        # Move conditional to beta only
        return np.sin(self.beta(theta, l_eff))**3

    def jd(self, theta, l_eff):
        return np.cos(theta) * self.U(theta, l_eff)

    def jw(self, theta, l_eff):
        result = (
        (4/(3*np.pi))*(1-1/(2*l_eff + 1)) * (1/l_eff) 
        * (np.cos(theta)**2 / np.sin(theta)) 
        * (1-self.V(theta, l_eff))
        + (1/(2*l_eff + 1))*np.cos(theta) * (1-self.U(theta, l_eff))
        )
        return result

    def j(self, theta, l_eff):
        return self.jd(theta, l_eff) + self.jw(theta, l_eff)


    def beam_profile(self, theta, l_eff, theta_max):
        # Theta needs to be connverted to array if it isn't already
        theta = np.array(np.abs(theta))
        cond = (theta < theta_max) & (theta != 0)
        result = np.piecewise(theta, 
            [theta == 0, cond, theta >= theta_max],
            [
                1,self.j(theta[cond], l_eff), 0   
            ]
                )
        return result
    
    # Introduce wire_sensitivity options

  
    def init_eta_wire_sim(self,):
        # Simulated (default as of 2023-12-11 (and for evaluuatins before))
        # The hacky way: Data from 1µW simulation from 2023-06-12
        signal_arr_norm_1 = [0.0216639499386916, 
                            0.06400510252472039, 0.1053120495242299, 
        0.1455846628698668, 0.22303108769945335, 0.260209306315826, 
        0.29636195737211446, 0.3314934324624902, 0.3987144288917815,
        0.4308166043914375, 0.46192281167714666, 0.5211797239090269,
            0.5493479545349781, 0.5765551496525174, 0.6028111285969292,
            0.652510228987278, 0.6759743261085402, 0.6985290436389583, 
            0.7409536239600266, 0.7608451610624817, 0.7798705689299693,
                0.7980405003681781, 0.8318557396310409, 0.8475214304508306,
                0.8623723390671673, 0.8896675428283924, 0.9021298623376048,
                    0.913813396881382, 0.9247262045852894, 0.9442697530284341,
        0.9529144383444618, 0.9608162550428692, 0.9744139623616798, 
        0.9801199485931034, 0.9851032386668377, 0.989367596091101, 
        0.9957519299061783, 0.9978767811933983, 0.9992924418760787, 
        1.0, 1.0,
        0.9992924418760787, 0.9978767811933983, 0.9957519299061783,
        0.989367596091101, 0.9851032386668377, 0.9801199485931034,
            0.9744139623616798, 0.9608162550428692, 0.9529144383444618, 
            0.9442697530284341, 0.9247262045852894, 0.913813396881382,
            0.9021298623376048, 0.8896675428283924, 0.8623723390671673,
                0.8475214304508306, 0.8318557396310409, 0.7980405003681781,
                0.7798705689299693, 0.7608451610624817, 0.7409536239600266,
        0.6985290436389583, 0.6759743261085402, 0.652510228987278, 
        0.6028111285969292, 0.5765551496525174, 0.5493479545349781,
        0.5211797239090269, 0.46192281167714666, 0.4308166043914375,
            0.3987144288917815, 0.3314934324624902, 0.29636195737211446,
            0.260209306315826, 0.22303108769945335, 0.1455846628698668,
                0.1053120495242299, 0.06400510252472039, 0.0216639499386916]

        x_pos_seg = np.array([-0.0099, -0.0097, -0.0095, -0.009300000000000001, 
                        -0.0089, -0.0087,
            -0.0085, -0.0083, -0.0079, -0.0077, -0.0075, 
            -0.0070999999999999995, 
            -0.0069, -0.006699999999999999, -0.006499999999999999, -0.0061,
            -0.0059, -0.005700000000000001, -0.0053, -0.0051, -0.0049, -0.0047, 
            -0.004299999999999999, -0.0041, -0.0039000000000000003, 
            -0.0034999999999999996, -0.0032999999999999995, 
            -0.0030999999999999995, 
            -0.0028999999999999994, -0.0024999999999999988, -0.0023, -0.0021, 
            -0.0017000000000000003, -0.0014999999999999992, 
            -0.0012999999999999989, 
            -0.0010999999999999998, -0.0006999999999999996, 
            -0.0005000000000000004, 
            -0.0003000000000000003, 0.00010000000000000231])
        # Mirror and convert to mm
        x_pos_seg = np.concatenate([x_pos_seg, -1 * x_pos_seg[::-1]]) * 1e3
        self.eta_wire_sim = interp1d(x_pos_seg, signal_arr_norm_1, 
                            fill_value = "extrapolate")
        return 


    # Define Integrated Power on wire with variable eta (wire sensitivity)
    # Adapted from HABS_bpf_4.1.3_2023-11-22_dual_composition_fit.ipynb
    #"P_int_eta_laser"
    # #### Reimplemented below as a special case of the 2 component function
    # def P_int(self,
    #           z_pos, l_eff, theta_max, z0, A, P_0,
    #           y0 = None, # Set to constant for simplicity 
    #           eta = None
    #                  ):
    #     if eta is None:
    #         eta = self.eta_wire_default
    #     if y0 is None:
    #         y0 = self.y0_default

    #     # Edit of P_int_fast_array with selectabel eta
    #     z_space = z_pos
    #     result = np.zeros_like(z_space)
    #     for i,z_pos in enumerate(z_space): 
    #         def theta(lw):
    #             return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) / y0)
    #         integrant = lambda lw: ((self.beam_profile(theta(lw), l_eff, 
    #                                                    theta_max) 
    #                                 * np.cos(theta(lw))**3) 
    #                                 # Multiply by eta_wire
    #                                 * eta(lw)
    #                                         )   
    #         result[i] = integrate.quad(integrant, -10, 10,
    #                             epsabs=1e-1, epsrel=1e-1)[0]
    #     return A * result + P_0
    

    # define dual fit 
    def P_int_2_component(self,
                        z_pos, l_eff, theta_max, z0, A, P_0,
                        l_eff_2,A_2,
                        y0 = None, # Set to constant for simplicity 
                        eta = None
                        ):
        if eta is None:
            eta = self.eta_wire_default
        if y0 is None:
            y0 = self.y0_default
        # Edit of P_int_fast_array with wire sensitivity from laser test
        z_space = z_pos
        result = np.zeros_like(z_space)
        for i,z_pos in enumerate(z_space): 
            def theta(lw):
                        return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) 
                                         / y0)
            integrant = lambda lw: (
                    (A * self.beam_profile(theta(lw), l_eff, theta_max) 
                    + A_2 * self.beam_profile(theta(lw), l_eff_2, theta_max))
                                    * np.cos(theta(lw))**3 
                                    # Multiply by eta_wire
                                    * eta(lw)
                                            )   
            result[i] = integrate.quad(integrant, -10, 10,
                                epsabs=1e-1, epsrel=1e-1)[0]
        return result + P_0   

    def P_int(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              y0 = None, # Set to constant for simplicity 
              eta = None
                     ):
        result = self.P_int_2_component(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              l_eff_2 = 0,A_2 = 0,
              y0 = y0, # Set to constant for simplicity 
              eta = eta)
        return result


    ############################
    # Plotting function



        


if __name__ == "__main__":
    # Test using 1 sccm,  1500K
    # run_name = "2023-09-15_1sccm_475TC_z-scan_jf+hg_wire"
    # Make fit run dictionary
    # TODO Should this not rather be a class of its oown, so it can suggest its
    #  options?
    run_dict = {}
    rd = run_dict
    run_dict["base_dir"] =  ("C:/Users/Christian/Documents/StudiumPhD/python/"
                            + "Keysight-DMM-34461A/analysis/")
    run_dict["plot_dir"] = (run_dict["base_dir"] + os.sep 
                            + "output/flow_on_off/")
    run_dict["sc_dir"] = (run_dict["base_dir"]
                            + os.sep + "../" 
                            + "SC_downloads/")
                        
    run_dict["run_name"] = "2023-09-15_1sccm_475TC_z-scan_jf+hg_wire"
    run_dict["data_name"] = run_dict["run_name"]


    # requires flow_on_off_cycle_analysis_2
    # Which in turn  requires others, so we need to update the entire 
    # series to put them in the package

    ext_dict = load_dict(rd["plot_dir"] + rd["run_name"] + os.sep 
                            + "ext_dict")

    run_dict[""]

