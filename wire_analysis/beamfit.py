import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.interpolate import interp1d

import os
import json

#######Imports from  other submodules
# Put  these into a  utils.py?
# THis is currently a bad way to do it
from .flow_on_off_cycle_analysis import (load_dict, make_result_dict, 
                                        sort_by_z_list)
from .utils import (load_json_dict, save_json_dict, load_extractor_dict_json)

########

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)
mpl.rcParams['figure.dpi'] = 400
#From Pascal
mpl.rcParams['text.usetex'] = True
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['font.family'] = 'cmu serif'


# make a class for fitting l_eff Tschersich model to beam shape
# Started 2023-Dec-11
# based on "HABS_beam_profile_fitting_4.1.1.1_cos3_low_temp_1sccm.ipynb"
# requirements: (non-exhaustive)
#   - allow for modular eta_wire
#   - Properly save outpput for further use
#   - Take output of flow_on_off_cycle_analysis_2.py
#
degree = np.pi/180 # convert form rad to degree

class Beamfit():
    """
    "The beamfit class contains functions and parameters necessary to
    fit beam profiles with the Tschersich model. 
    (K. G. Tschersich; V. von Bonin (1998))
    "

    Parameters
    ----------
    # y0_default :  `float`
    #     default distance from wire to source [mm]
    run_dict_path = `str`
        path string to the run_dict.json which contains the information how
        to run the analysis, and where thefit reslsts are saved.
        Contents of run_dict:

    
    """
    def __init__(self,
                run_dict_path = ""
                ):
        # Check that valid run_dict path was entered
        if os.path.isfile(run_dict_path):
            pass
        else:
            print("please enter a valid run_dict_path")
            exit()
        #define parameters
        self.run_dict_path = run_dict_path


        #define constants
        degree = np.pi/180 # convert form rad to degree
        self.degree = degree
        # run init fucntions
        self.init_eta_wire_sim()
            #Allow for later adjustment of default eta 
        self.eta_wire_default = self.eta_wire_sim
        self.init_penumbra()

        #load run_dict
        self.run_dict = load_json_dict(run_dict_path)
        rd = self.run_dict
        self.y0_default =  self.run_dict["fit_start"]["y0"] 
            # mm (estimated from CAD) not  fitted by default

        # Initialize output dir
        self.out_dir = rd["out_dir_base"] + rd["analysis_run_name"] + os.sep 
        os.makedirs(self.out_dir, exist_ok=True)
        # initialize run_dict json name and save it into the out_dir
        # Chop filename off the end of the run_dict path
        self.run_dict_name = os.path.split(self.run_dict_path)[1]

        return


    # def save_json_dict(self, dict_path, dict):
    #     with open((dict_path), 'w', encoding='utf-8') as f:
    #         json.dump(dict, f, ensure_ascii=False, indent=4)
    #     return

    # def load_json_dict(self, dict_path):
    #     with open(dict_path, 'r', encoding='utf-8') as f:
    #         dict_load = json.load(f)
    #     return dict_load

    #####HACK to incuce calc_norm
    def integrate_H_angles_1D(self, theta_lim = [0,np.pi/2], l_eff = 7.96, y0 = 35.17,
                        theta_max = 90 * degree, norm_factor = 1):
        # Integral over phi  from [0, 2 * np.pi]
        integrant = lambda theta: (norm_factor * 2 * np.pi
                            * self.beam_profile(theta, l_eff, theta_max) 
                                * np.sin(theta)
                            )
        result, err = integrate.quad(integrant,
                                    theta_lim[0], theta_lim[1] # theta_lims
                                    )
        return result
        
    def calc_norm_factor(self, l_eff
                        ):
        return 1/self.integrate_H_angles_1D(l_eff = l_eff)
    ##################

    def save_json_run_dict(self, dict_path = None, dict = None):
        # If no values are supplied save the current run_dict to the current 
        # run_dict_path
        if dict_path is None:
            dict_path = self.run_dict_path
        if dict is None:
            dict = self.run_dict
        save_json_dict(dict_path, dict)
        #also  save to out_dir, so anything that is saved elsewhere can be
        # identified with its run
        save_json_dict(self.out_dir + self.run_dict_name, dict)
        return

    def load_json_dict(self, dict_path):
        with open(dict_path, 'r', encoding='utf-8') as f:
            dict_load = json.load(f)
        return dict_load



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

    def beam_profile_jw(self, theta, l_eff, theta_max):
        # Theta needs to be connverted to array if it isn't already
        theta = np.array(np.abs(theta))
        cond = (theta < theta_max) & (theta != 0)
        result = np.piecewise(theta, 
            [theta == 0, cond, theta >= theta_max],
            [
                1,self.jw(theta[cond], l_eff), 0   
            ]
                )
        return result
    
    def beam_profile_eibl062(self, theta, l_eff, theta_max):
        # Theta needs to be connverted to array if it isn't already
        theta = np.array(np.abs(theta))
        cond = (theta < theta_max) & (theta != 0)
        result = np.piecewise(theta, 
            [theta == 0, cond, theta >= theta_max],
            [
                1,
                0.62*self.jd(theta[cond], l_eff) 
                + self.jw(theta[cond], l_eff),
                0   
            ]
                )
        return result
    
    def beam_profile_eibl062_inverse(self, theta, l_eff, theta_max):
        # Theta needs to be connverted to array if it isn't already
        theta = np.array(np.abs(theta))
        cond = (theta < theta_max) & (theta != 0)
        result = np.piecewise(theta, 
            [theta == 0, cond, theta >= theta_max],
            [
                1,
                self.jd(theta[cond], l_eff) 
                + 0.62*self.jw(theta[cond], l_eff),
                0   
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
    
    def init_penumbra(self,radius_cap = 0.5, 
                      distance_cap_housing = 7.1, radius_housing = 3.0):
        # Approximates visible capillary fraction as a straight line passing
        # over a circle
        # and linearly interpolates the angles within the penumbra

        # calculate fraction of the capillary that is visible under a
        # given angle
        # Presumes capillary to be centered on openning  in HABS housing
        # radius_cap = 0.5 # mm
        # distance_cap_housing = 7.1 # mm
        # radius_housing = 3.0 # mm

        # calculate edges of penumbra
        # inner (close) edge
        theta_inner = np.arctan((radius_housing - radius_cap)
                                /distance_cap_housing)
        theta_outer = np.arctan((radius_housing + radius_cap)
                                /distance_cap_housing)

        degree = np.pi/180 # convert form rad to degree
        # print("inner edge:", theta_inner, theta_inner / degree)
        # print("outer edge:", theta_outer, theta_outer / degree)#
        self.theta_inner = theta_inner
        self.theta_outer = theta_outer
        # speed up by eliminating  this calculation
        # self.theta_inner = 0.3385556949116842 # 19.3978 deg 
        # self.theta_outer = 0.45799795159722173 # 26.2413 deg
        return (theta_inner, theta_outer)


    # Define Integrated Power on wire with variable eta (wire sensitivity)
    # Adapted from HABS_bpf_4.1.3_2023-11-22_dual_composition_fit.ipynb
    #"P_int_eta_laser"
    # #### Reimplemented below as a special case of the 2 component function
    def P_int(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              y0 = None, # Set to constant for simplicity 
              eta = None
                     ):
        if eta is None:
            eta = self.eta_wire_default
        if y0 is None:
            y0 = self.y0_default

        # Edit of P_int_fast_array with selectabel eta
        z_space = z_pos
        result = np.zeros_like(z_space)
        for i,z_pos in enumerate(z_space): 
            def theta(lw):
                return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) / y0)
            integrant = lambda lw: ((self.beam_profile(theta(lw), l_eff, 
                                                       theta_max) 
                                    * np.cos(theta(lw))**3) 
                                    # Multiply by eta_wire
                                    * eta(lw)
                                            )   
            result[i] = integrate.quad(integrant, -10, 10,
                                epsabs=1e-1, epsrel=1e-1)[0]
        return A * result + P_0
    

    def visible_fraction(self, theta, 
                         theta_inner =  None, theta_outer = None):
        # HACK the current implementation is merely an approximation
        if theta_inner is None:
            theta_inner = self.theta_inner
            theta_outer = self.theta_outer
        if np.abs(theta) <= theta_inner:
            fraction = 1
        elif np.abs(theta) >= theta_outer:
            fraction = 0
        else:
            # Take shortcut initially
            # HACK Linearly interpolate between inner and outer angles
            int_var = (2 *  ((theta_outer - np.abs(theta)) 
                        / (theta_outer - theta_inner))) - 1
            # Integrated circle equation
            vis_frac = lambda x: (np.sqrt(1-x**2) * x 
                                    + np.arcsin(x)
                                    + np.pi/2)/np.pi
            fraction = vis_frac(int_var)
        return fraction

    def P_int_penumbra(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              y0 = None, # Set to constant for simplicity 
              eta = None
                     ):
        # Attempt to model that the Capillary gets obscured by its housing 
        # in parts and not all at once
        if eta is None:
            eta = self.eta_wire_default
        if y0 is None:
            y0 = self.y0_default

        # Edit of P_int_fast_array with selectabel eta
        z_space = z_pos
        result = np.zeros_like(z_space)
        for i,z_pos in enumerate(z_space): 
            def theta(lw):
                return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) / y0)
            integrant = lambda lw: ((self.beam_profile(theta(lw), l_eff, 
                                                       theta_max) 
                                    * np.cos(theta(lw))**3) 
                                    # apply "visible capillary fraction"
                                    * self.visible_fraction(theta(lw))
                                    # Multiply by eta_wire
                                    * eta(lw)
                                            )   
            result[i] = integrate.quad(integrant, -10, 10,
                                epsabs=1e-1, epsrel=1e-1)[0]
        return A * result + P_0
    
    def P_int_penumbra_dch(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              d_ch,
              y0 = None, # Set to constant for simplicity 
              eta = None
                     ):
        # Attempt to model that the Capillary gets obscured by its housing 
        # in parts and not all at once
        if eta is None:
            eta = self.eta_wire_default
        if y0 is None:
            y0 = self.y0_default
        
        distance_cap_housing = d_ch
        radius_cap = 0.5 # mm
        radius_housing = 3.0 # mm
        theta_inner = np.arctan((radius_housing - radius_cap)
                                /distance_cap_housing)
        theta_outer = np.arctan((radius_housing + radius_cap)
                                /distance_cap_housing)

        # Edit of P_int_fast_array with selectabel eta
        z_space = z_pos
        result = np.zeros_like(z_space)
        for i,z_pos in enumerate(z_space): 
            def theta(lw):
                return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) / y0)
            integrant = lambda lw: ((self.beam_profile(theta(lw), l_eff, 
                                                       theta_max) 
                                    * np.cos(theta(lw))**3) 
                                    # apply "visible capillary fraction"
                                    * self.visible_fraction(theta(lw),
                                        theta_inner = theta_inner, 
                                        theta_outer = theta_outer)
                                    # Multiply by eta_wire
                                    * eta(lw)
                                            )   
            result[i] = integrate.quad(integrant, -10, 10,
                                epsabs=1e-1, epsrel=1e-1)[0]
        return A * result + P_0
    
    def P_int_penumbra_3par(self,
              z_pos, l_eff, theta_max, z0, A, P_0,
              d_ch, r_h, r_c,
              y0 = None, # Set to constant for simplicity 
              eta = None
                     ):
        # Attempt to model that the Capillary gets obscured by its housing 
        # in parts and not all at once
        if eta is None:
            eta = self.eta_wire_default
        if y0 is None:
            y0 = self.y0_default
        
        distance_cap_housing = d_ch
        radius_cap = r_c # mm
        radius_housing = r_h # mm
        theta_inner = np.arctan((radius_housing - radius_cap)
                                /distance_cap_housing)
        theta_outer = np.arctan((radius_housing + radius_cap)
                                /distance_cap_housing)


        # Edit of P_int_fast_array with selectabel eta
        z_space = z_pos
        result = np.zeros_like(z_space)
        for i,z_pos in enumerate(z_space): 
            def theta(lw):
                return np.arctan(np.sqrt(lw ** 2 + (z_pos - z0) ** 2) / y0)
            integrant = lambda lw: ((self.beam_profile(theta(lw), l_eff, 
                                                       theta_max) 
                                    * np.cos(theta(lw))**3) 
                                    # apply "visible capillary fraction"
                                    * self.visible_fraction(theta(lw),
                                        theta_inner = theta_inner, 
                                        theta_outer = theta_outer)
                                    # Multiply by eta_wire
                                    * eta(lw)
                                            )   
            result[i] = integrate.quad(integrant, -10, 10,
                                epsabs=1e-1, epsrel=1e-1)[0]
            
        #WARNIN NOTE NOTE NOTE 
        #changing this  and only this functionto be normalized
        nf  = self.calc_norm_factor(l_eff = l_eff)
        d_wire = 5e-6 # wire thickness not included in P_fit (folded into A)
        # given in m
        #y0 = 35.17 # fit in mm conver to m -> times 1e-3
        #Integration was performed in mm -> 
        # need to include another factor of 1e-3 from dx[mm]=1e-3*dx[m]
        result = result /(1e-3 * (((y0*1e-3)**2)/(nf*d_wire)))
        return A * result + P_0

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
    # Redefine P_int as special case of 2 component model
    # TypeError: P_int_2_component() got multiple values for argument 'l_eff_2'
    # def P_int(self,
    #           z_pos, l_eff, theta_max, z0, A, P_0,
    #           y0 = None, # Set to constant for simplicity 
    #           eta = None
    #                  ):
    #     result = self.P_int_2_component(self,
    #           z_pos, l_eff, theta_max, z0, A, P_0,
    #           l_eff_2 = 10,A_2 = 0,
    #           y0 = y0, # Set to constant for simplicity 
    #           eta = eta)
    #     return result
    

    #################################################

    # save to run dict:
    def save_fit_results(self, popt_abs, pcov_abs,
                         fit_y0=False, fit_dch=False):
        rd = self.run_dict
        rd["fit_result"] = {}
        rd["fit_result_errors"] = {}
        ############
        rd["fit_result"]["l_eff"] = popt_abs[0]
        rd["fit_result"]["A"] = popt_abs[1]
        rd["fit_result"]["z0"] = popt_abs[2]
        rd["fit_result"]["P_0"] = popt_abs[3]
        # currently not fitted
        rd["fit_result"]["theta_max"] = rd["fit_start"]["theta_max"]
        rd["fit_result"]["y0"] = rd["fit_start"]["y0"]

        # Also save errors
        # cannot be numpy array, to make json serializable
        rd["fit_result_errors"]["l_eff"] = np.sqrt(pcov_abs[0,0])
        rd["fit_result_errors"]["A"] = np.sqrt(pcov_abs[1,1])
        rd["fit_result_errors"]["z0"] = np.sqrt(pcov_abs[2,2])
        rd["fit_result_errors"]["P_0"] = np.sqrt(pcov_abs[3,3])

        # currently not fitted
        rd["fit_result_errors"]["theta_max"] = None
        rd["fit_result_errors"]["y0"] = None

        # HACK
        if fit_y0:
            rd["fit_result"]["y0"] = popt_abs[4]
            rd["fit_result_errors"]["y_0"] = np.sqrt(pcov_abs[4,4])

        #HACK
        if fit_dch:
            rd["fit_result"]["d_ch"] = popt_abs[4]
            rd["fit_result_errors"]["d_ch"] = np.sqrt(pcov_abs[4,4])

        # Save to file
        self.save_json_run_dict(dict_path= self.run_dict_path, 
                            dict = self.run_dict)
        return

    def default_fit(self,
                     ):
        rd = self.run_dict
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

        [a,b,c] = rd["selection_indices"]
        #neglegt leading selected points
        P_arr = P_arr[a:b:c]
        P_err_arr = P_err_arr[a:b:c]
        z_arr = z_arr[a:b:c]

        # initiate fit parameters
        l_eff = rd["fit_start"]["l_eff"]
        theta_max = rd["fit_start"]["theta_max"]
        z0 = rd["fit_start"]["z0"]
        A = rd["fit_start"]["A"]
        P_0 = rd["fit_start"]["P_0"]
        A_bound = rd["fit_start"]["A_bound"]
        # ##### All parameters except theta_max
        P_int_fit = lambda z_space, l_eff, A , z0, P_0: self.P_int(
                        z_space, l_eff, theta_max, z0, A, P_0)

        # Use errors as absolute to get proper error estimation
        popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                            sigma = P_err_arr,
                            absolute_sigma= False,
                            p0 = [l_eff, A,  z0, P_0], 
                            bounds=([2, A_bound[0],  z0 - 1, P_0 - 0.1],
                                    [20,  A_bound[1],  z0 + 1, P_0 + 0.1])
                            )
        ######

        print(popt_abs, pcov_abs)
        for i, p in enumerate(popt_abs):
            print(f"parameter {i:.0f}: {p:.5f}"
                +f"+-{np.sqrt(pcov_abs[i,i]):.5f}")
        # save to file
        self.save_fit_results(popt_abs, pcov_abs)
            
        #### plot
        z_space = np.linspace(-11,20,num=100)

        P_space_eye= P_int_fit(z_space, *popt_abs)
        P_arr_eye = P_int_fit(z_arr, *popt_abs)

        # Angle plot
        self.plot_fit(P_arr, P_arr_eye, P_err_arr, z_arr
            , z_space,P_space_eye, scale_residuals=True, 
            plot_angles=True, z0=popt_abs[2],
            theta_max=theta_max/self.degree,
            l_eff_str =(f"{popt_abs[0]:.2f}"+ r"$\pm$"
                     + f"{np.sqrt(pcov_abs[0,0]):.2f}"),
                     plotname = "default_fit_plot")
        return

    def custom_data_fit(self,z_arr, p_arr, p_err_arr, 
                        plotname = "custom_data_fit_plot"
                     ):
        # Crude function for fitting a dataset not equal to the run dict data
        # file. I suppose that really should not be the intended use case. 
        # mmmmmmh
        rd = self.run_dict
        # ######################################
        # # Load from newly reformatted result dict files
        # extractor_dict_unsorted = load_extractor_dict_json(
        #                         rd["extractor_dict_path"])
        # HACK brew up an equivalently formated dict:
        sort_this_dict = {} 
        sort_this_dict["p_arr"] = p_arr 
        sort_this_dict["p_err_arr"] = p_err_arr 
        # ##############
        z_array_unsorted = z_arr

        # Sort these:
        z_arr, sorted_dict = sort_by_z_list(z_array_unsorted,
                                            sort_this_dict)

        P_arr = sorted_dict["p_arr"]
        P_err_arr = sorted_dict["p_err_arr"]

        [a,b,c] = rd["selection_indices"]
        #neglegt leading selected points
        P_arr = P_arr[a:b:c]
        P_err_arr = P_err_arr[a:b:c]
        z_arr = z_arr[a:b:c]

        # initiate fit parameters
        l_eff = rd["fit_start"]["l_eff"]
        theta_max = rd["fit_start"]["theta_max"]
        z0 = rd["fit_start"]["z0"]
        A = rd["fit_start"]["A"]
        P_0 = rd["fit_start"]["P_0"]
        A_bound = rd["fit_start"]["A_bound"]
        # ##### All parameters except theta_max
        P_int_fit = lambda z_space, l_eff, A , z0, P_0: self.P_int(
                        z_space, l_eff, theta_max, z0, A, P_0)

        # Use errors as absolute to get proper error estimation
        popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                            sigma = P_err_arr,
                            absolute_sigma= False,
                            p0 = [l_eff, A,  z0, P_0], 
                            bounds=([2, A_bound[0],  z0 - 1, P_0 - 0.1],
                                    [20,  A_bound[1],  z0 + 1, P_0 + 0.1])
                            )
        ######

        print(popt_abs, pcov_abs)
        for i, p in enumerate(popt_abs):
            print(f"parameter {i:.0f}: {p:.5f}"
                +f"+-{np.sqrt(pcov_abs[i,i]):.5f}")
        # save to file
        self.save_fit_results(popt_abs, pcov_abs)
            
        #### plot
        z_space = np.linspace(-11,20,num=100)

        P_space_eye= P_int_fit(z_space, *popt_abs)
        P_arr_eye = P_int_fit(z_arr, *popt_abs)

        # Angle plot
        self.plot_fit(P_arr, P_arr_eye, P_err_arr, z_arr
            , z_space,P_space_eye, scale_residuals=True, 
            plot_angles=True, z0=popt_abs[2],
            theta_max=theta_max/self.degree,
            l_eff_str =(f"{popt_abs[0]:.2f}"+ r"$\pm$"
                     + f"{np.sqrt(pcov_abs[0,0]):.2f}"),
                     plotname = plotname)
        return
    
    def custom_fit(self,z_arr, p_arr, p_err_arr, 
                        plotname = "custom_fit_plot"
                     ):
        # Crude function for fitting a dataset not equal to the run dict data
        # file. I suppose that really should not be the intended use case. 
        # mmmmmmh
        rd = self.run_dict
        # ######################################
        # # Load from newly reformatted result dict files
        # extractor_dict_unsorted = load_extractor_dict_json(
        #                         rd["extractor_dict_path"])
        # HACK brew up an equivalently formated dict:
        sort_this_dict = {} 
        sort_this_dict["p_arr"] = p_arr 
        sort_this_dict["p_err_arr"] = p_err_arr 
        # ##############
        z_array_unsorted = z_arr

        # Sort these:
        z_arr, sorted_dict = sort_by_z_list(z_array_unsorted,
                                            sort_this_dict)

        P_arr = sorted_dict["p_arr"]
        P_err_arr = sorted_dict["p_err_arr"]

        [a,b,c] = rd["selection_indices"]
        #neglegt leading selected points
        P_arr = P_arr[a:b:c]
        P_err_arr = P_err_arr[a:b:c]
        z_arr = z_arr[a:b:c]

        # initiate fit parameters
        l_eff = rd["fit_start"]["l_eff"]
        theta_max = rd["fit_start"]["theta_max"]
        z0 = rd["fit_start"]["z0"]
        A = rd["fit_start"]["A"]
        P_0 = rd["fit_start"]["P_0"]
        A_bound = rd["fit_start"]["A_bound"]
        # ##### All parameters except theta_max
        P_int_fit = lambda z_space, l_eff, A , z0, P_0: self.P_int_penumbra(
                        z_space, l_eff, theta_max, z0, A, P_0)

        # Use errors as absolute to get proper error estimation
        popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                            sigma = P_err_arr,
                            absolute_sigma= False,
                            p0 = [l_eff, A,  z0, P_0], 
                            bounds=([2, A_bound[0],  z0 - 1, P_0 - 0.1],
                                    [20,  A_bound[1],  z0 + 1, P_0 + 0.1])
                            )
        ######

        print(popt_abs, pcov_abs)
        for i, p in enumerate(popt_abs):
            print(f"parameter {i:.0f}: {p:.5f}"
                +f"+-{np.sqrt(pcov_abs[i,i]):.5f}")
        # save to file
        self.save_fit_results(popt_abs, pcov_abs)
            
        #### plot
        z_space = np.linspace(-11,20,num=100)

        P_space_eye= P_int_fit(z_space, *popt_abs)
        P_arr_eye = P_int_fit(z_arr, *popt_abs)

        # Angle plot
        self.plot_fit(P_arr, P_arr_eye, P_err_arr, z_arr
            , z_space,P_space_eye, scale_residuals=True, 
            plot_angles=True, z0=popt_abs[2],
            theta_max=theta_max/self.degree,
            l_eff_str =(f"{popt_abs[0]:.2f}"+ r"$\pm$"
                     + f"{np.sqrt(pcov_abs[0,0]):.2f}"),
                     plotname = plotname)
        return (popt_abs, pcov_abs)
    
    def fit_d_ch(self,z_arr, p_arr, p_err_arr, 
                        plotname = "fit_d_ch_plot",
                        fit_y0 = False
                     ):
        # Crude function for fitting a dataset not equal to the run dict data
        # file. I suppose that really should not be the intended use case. 
        # mmmmmmh
        rd = self.run_dict
        # ######################################
        # # Load from newly reformatted result dict files
        # extractor_dict_unsorted = load_extractor_dict_json(
        #                         rd["extractor_dict_path"])
        # HACK brew up an equivalently formated dict:
        sort_this_dict = {} 
        sort_this_dict["p_arr"] = p_arr 
        sort_this_dict["p_err_arr"] = p_err_arr 
        # ##############
        z_array_unsorted = z_arr

        # Sort these:
        z_arr, sorted_dict = sort_by_z_list(z_array_unsorted,
                                            sort_this_dict)

        P_arr = sorted_dict["p_arr"]
        P_err_arr = sorted_dict["p_err_arr"]

        [a,b,c] = rd["selection_indices"]
        #neglegt leading selected points
        P_arr = P_arr[a:b:c]
        P_err_arr = P_err_arr[a:b:c]
        z_arr = z_arr[a:b:c]

        # initiate fit parameters
        l_eff = rd["fit_start"]["l_eff"]
        theta_max = rd["fit_start"]["theta_max"]
        z0 = rd["fit_start"]["z0"]
        A = rd["fit_start"]["A"]
        P_0 = rd["fit_start"]["P_0"]
        d_ch_0 = 7.1

        A_bound = rd["fit_start"]["A_bound"]
        # ##### All parameters except theta_max
        # P_int_fit = lambda z_space, l_eff, A , z0, P_0, d_ch: (
        #         self.P_int_penumbra_dch(
        #             z_space, l_eff, theta_max, z0, A, P_0, d_ch))

        # HACK to include r_c as well
        r_h_0 = 3.0
        #r_c_0 = 0.5
        # Guess that r_c_0 is actually a little smaller due to slight central
        # concentration  of gas
        r_c_0 = 0.4
        # ##### All parameters except theta_max
        if fit_y0:
            fit_dch = False
            # P_int_fit = lambda z_space, l_eff, A , z0, P_0, d_ch, y0: (
            #     self.P_int_penumbra_3par(
            #         z_space, l_eff, theta_max, z0, A, P_0, d_ch, r_h_0, r_c_0,
            #         y0 = y0))
            P_int_fit = lambda z_space, l_eff, A , z0, P_0, y0: (
                self.P_int_penumbra_3par(
                    z_space, l_eff, theta_max, z0, A, P_0, d_ch_0, r_h_0, r_c_0,
                    y0 = y0))

            
            # HACK bounds
            # Use errors as absolute to get proper error estimation
            y_base = self.y0_default
            popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                                sigma = P_err_arr,
                                absolute_sigma= False,
                                p0 = [l_eff, A,  z0, P_0, y_base], 
                                bounds=([1, A_bound[0]*0.3,  z0 - 1, P_0 - 0.1,
                                         y_base -5],
                                        [20,  A_bound[1],  z0 + 1, P_0 + 0.1,
                                         y_base +20])
                                )
        else:
            fit_dch = True
            P_int_fit = lambda z_space, l_eff, A , z0, P_0, d_ch: (
                self.P_int_penumbra_3par(
                    z_space, l_eff, theta_max, z0, A, P_0, d_ch, r_h_0, r_c_0))

            # Use errors as absolute to get proper error estimation
            popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                                sigma = P_err_arr,
                                absolute_sigma= False,
                                p0 = [l_eff, A,  z0, P_0, d_ch_0], 
                                bounds=([1, A_bound[0]*0.3,  z0 - 1, P_0 - 0.1,
                                        d_ch_0 - 1.5],
                                        [20,  A_bound[1],  z0 + 1, P_0 + 0.1,
                                        d_ch_0 + 1.5])
                                )
        ######

        print(popt_abs, pcov_abs)
        for i, p in enumerate(popt_abs):
            print(f"parameter {i:.0f}: {p:.5f}"
                +f"+-{np.sqrt(pcov_abs[i,i]):.5f}")
        # save to file
        #HAck to save either y_0 or d_chfit parameter
        self.save_fit_results(popt_abs, pcov_abs, 
                            fit_y0=fit_y0, fit_dch=fit_dch)
            
        #### plot
        z_space = np.linspace(-11,20,num=100)

        P_space_eye= P_int_fit(z_space, *popt_abs)
        P_arr_eye = P_int_fit(z_arr, *popt_abs)

        # Option to add table of parameters
        if True:
            mpl.rc('text', usetex=True)
            # col_labels=["parameter",'fit value']
            # row_labels=[r'$l_{\rm eff}$',r'$A$',r'$z_0$', r"$P_0$",
            #             r"$d_{ch}$"]
           
            # table_vals=[[f"{popt_abs[i]:.2f}"+ r"$\pm$"
            #          + f"{np.sqrt(pcov_abs[i,i]):.2f}"
            #          if np.sqrt(pcov_abs[i,i])>0.01
            #          else f"{popt_abs[i]:.2f}"+ r"$\pm$"
            #          + f"{np.sqrt(pcov_abs[i,i]):.1e}"] 
            #          for i in range(len(popt_abs))]
            
            col_labels=["parameter",'value', 'error']
            if fit_y0:
                # row_labels=[r'$l_{\rm eff}$',r'$A$',r'$z_0$', r"$P_0$",
                #         r"$d_{ch}$", r"$y_{0}$"]
                row_labels=[r'$l_{\rm eff}$',r'$A$',r'$z_0$', r"$P_0$",
                            r"$y_{0}$"]
            else:
                row_labels=[r'$l_{\rm eff}$',r'$A$',r'$z_0$', r"$P_0$",
                            r"$d_{ch}$"]
                # unit_labels = [r'1',r'${\rm µW}/{\rm mm}^2$', r'mm', 
                #             r'$\mu {\rm W}$', r'mm']
            
            table_vals=[[f"{popt_abs[i]:.2f}"
                     , f"{np.sqrt(pcov_abs[i,i]):.2f}"]
                     if np.sqrt(pcov_abs[i,i])>0.01
                     else [f"{popt_abs[i]:.2f}"
                     , f"{np.sqrt(pcov_abs[i,i]):.0e}"]

                    if popt_abs[i]>0.095
                    else [f"{popt_abs[i]:.1e}"
                     , f"{np.sqrt(pcov_abs[i,i]):.0e}"]

                     for i in range(len(popt_abs))]
            print("table_vals: ", table_vals)
            table = r'''\begin{tabular}{ c''' + (len(col_labels)-1)*" | c" +"}"
            #add column headers
            table = table + " & ".join(col_labels) + r" \\ \hline"
            # Add rows:
            for row in range(len(row_labels)):
                table = table + row_labels[row] + " & "
                table = table + " & ".join(table_vals[row])
                #line break
                table = table + r" \\ \hline"
            #end table
            table = table + r"\end{tabular}"
            
            print(table)

            # plt.text(0,0,table,size=12,
            # horizontalalignment='left',
            # verticalalignment='bottom', transform=ax1.transAxes)

            # # the rectangle is where I want to place the table
            # the_table = plt.table(cellText=table_vals,
            #         colWidths = [0.1]*len(col_labels),
            #         rowLabels=row_labels,
            #         colLabels=col_labels,
            #         loc='center right')

        # Angle plot
        self.plot_fit(P_arr, P_arr_eye, P_err_arr, z_arr
            , z_space,P_space_eye, scale_residuals=True, 
            plot_angles=True, z0=popt_abs[2],
            theta_max=theta_max/self.degree,
            l_eff_str = None,
                     plotname = plotname,
                     table_string = table)
        return (popt_abs, pcov_abs)
    
    def fit_3par_penumbra(self,z_arr, p_arr, p_err_arr, 
                        plotname = "fit_3par_penumb_plot"
                     ):
        # Crude function for fitting a dataset not equal to the run dict data
        # file. I suppose that really should not be the intended use case. 
        # mmmmmmh
        rd = self.run_dict
        # ######################################
        # # Load from newly reformatted result dict files
        # extractor_dict_unsorted = load_extractor_dict_json(
        #                         rd["extractor_dict_path"])
        # HACK brew up an equivalently formated dict:
        sort_this_dict = {} 
        sort_this_dict["p_arr"] = p_arr 
        sort_this_dict["p_err_arr"] = p_err_arr 
        # ##############
        z_array_unsorted = z_arr

        # Sort these:
        z_arr, sorted_dict = sort_by_z_list(z_array_unsorted,
                                            sort_this_dict)

        P_arr = sorted_dict["p_arr"]
        P_err_arr = sorted_dict["p_err_arr"]

        [a,b,c] = rd["selection_indices"]
        #neglegt leading selected points
        P_arr = P_arr[a:b:c]
        P_err_arr = P_err_arr[a:b:c]
        z_arr = z_arr[a:b:c]

        # initiate fit parameters
        l_eff = rd["fit_start"]["l_eff"]
        theta_max = rd["fit_start"]["theta_max"]
        z0 = rd["fit_start"]["z0"]
        A = rd["fit_start"]["A"]
        P_0 = rd["fit_start"]["P_0"]
        d_ch_0 = 7.1
        r_h_0 = 3.0
        r_c_0 = 0.5

        A_bound = rd["fit_start"]["A_bound"]
        # ##### All parameters except theta_max
        P_int_fit = lambda z_space, l_eff, A , z0, P_0, d_ch, r_h, r_c: (
            self.P_int_penumbra_3par(
                z_space, l_eff, theta_max, z0, A, P_0, d_ch, r_h, r_c))

        # Use errors as absolute to get proper error estimation
        popt_abs, pcov_abs = curve_fit(P_int_fit, z_arr, P_arr,
                            sigma = P_err_arr,
                            absolute_sigma= False,
                            p0 = [l_eff, A,  z0, P_0, d_ch_0, r_h_0, r_c_0], 
                            bounds=([2, A_bound[0],  z0 - 1, P_0 - 0.1,
                                     d_ch_0 - 1, r_h_0 * 0.8, r_c_0 * 0.2],
                                    [20,  A_bound[1],  z0 + 1, P_0 + 0.1,
                                     d_ch_0 + 5, r_h_0 * 1.5, r_c_0 * 1.2])
                            )
        ######

        print(popt_abs, pcov_abs)
        for i, p in enumerate(popt_abs):
            print(f"parameter {i:.0f}: {p:.5f}"
                +f"+-{np.sqrt(pcov_abs[i,i]):.5f}")
        # save to file
        self.save_fit_results(popt_abs, pcov_abs)
            
        #### plot
        z_space = np.linspace(-11,20,num=100)

        P_space_eye= P_int_fit(z_space, *popt_abs)
        P_arr_eye = P_int_fit(z_arr, *popt_abs)

        # Angle plot
        self.plot_fit(P_arr, P_arr_eye, P_err_arr, z_arr
            , z_space,P_space_eye, scale_residuals=True, 
            plot_angles=True, z0=popt_abs[2],
            theta_max=theta_max/self.degree,
            l_eff_str =(f"{popt_abs[0]:.2f}"+ r"$\pm$"
                     + f"{np.sqrt(pcov_abs[0,0]):.2f}"),
                     plotname = plotname)
        return (popt_abs, pcov_abs)


    #define plottign function:
    # TODO Rework: its just a straight copy for now
    def plot_fit(self,
                 P_arr, P_arr_eye, P_err_arr, z_arr, z_space, P_space_eye, 
              scale_residuals = False, 
             plot_angles = False, z0 = 1.5, theta_max = 22.6, 
             l_eff_str = None,
             # Add possibility to compare to second fit
             P_space_compare = None, P_compare_label = "P_compare",
            # Name plot
            plotname = "default",
            table_string = None,
             ):
        nParams = 5
        dof = len(P_arr) - nParams
        #chi2 = np.sum((P_arr - P_arr_eye)**2 / P_err_arr**2)
        chi2_red = np.sum((P_arr - P_arr_eye)**2 / P_err_arr**2) / dof

        fig = plt.figure(0, figsize=(8,6.5), dpi =300)
        ax1=plt.gca()
        gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1]) 
        gs.update(#wspace=0.05
                hspace = 0.005
            )

        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])


        if plot_angles == True:
            y0 = self.y0_default
            x_label = r"Central angle [deg]"
            z_arr = np.arctan((z_arr - z0) / y0) / self.degree
            z_space = np.arctan((z_space - z0) / y0) / self.degree
            # HACK for paper plot
            # ax1.axvline(theta_max, label=r"$\theta_{max} = $"
            #             + f"{theta_max:.2f} [deg]", 
            #             color = "gray", alpha = 0.6, ls = "--")
            # ax2.axvline(theta_max, label=r"$\theta_{max}$", 
            #             color = "gray", alpha = 0.6, ls = "--")
        else:
            x_label = r"$z_{pos}$ [mm]"
        ### ax1
        # ax1.errorbar(z_arr, P_arr,yerr = P_err_arr, fmt = ".",
        #             label = (r"data, $\chi^2_{red}$"
        #                     +" = {:2.3f} ".format(chi2_red)))
        # HACK for paper plot
        ax1.errorbar(z_arr, P_arr,yerr = P_err_arr, fmt = ".",
            label = (r"data".format(chi2_red)))

        # Plot Fit
        #xdata=p_arr
        if l_eff_str == None:
            ax1.plot(z_space, P_space_eye, "r-", label = r"$P_{\rm fit}$")
        else:
            # ax1.plot(z_space, P_space_eye, "r-", label = r"$P_{fit}$" 
            #         + r", $l_{eff}=$" + l_eff_str)
            # HACK for paper plot:
            ax1.plot(z_space, P_space_eye, "r-", 
                     label = r"$P_{\rm{rec}}(l_{\rm{eff}}=$" + l_eff_str + ")"
                   )

            #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [$\Omega$/µW],
            #  R$_0$={:5.3f} $\pm$ {:2.1e} [$\Omega$]'''.format(
            #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
            #           ))
                    
        if P_space_compare is None:
            pass
        else:
            ax1.plot(z_space, P_space_compare, "-",color = "C1", 
                    label = P_compare_label
                    )


        ax1.set_ylabel(r"Power [µW]")
        ax1.set_xlabel(x_label)

        ax1.grid(True)
        ax1.legend(shadow=True, fontsize = 13)
        # ax1.tight_layout()

        #ax 2
        if scale_residuals == False:
            ax2.errorbar(z_arr, P_arr-P_arr_eye, yerr = P_err_arr, fmt= ".",
                            label = "residuals")
            # Plot Fit
            #xdata=p_arr
            ax2.plot(z_space, [0.0 for z in z_space] , 'r-',
                #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW ],
                #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
                #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
                #           ))
                        )
            ax2.set_ylabel(r"Residuals [µW]")
        else:
            ax2.errorbar(z_arr, (P_arr-P_arr_eye)/P_err_arr, 
                        yerr = 1, fmt= ".",
                            label = "residuals")
            # Plot Fit
            #xdata=p_arr
            ax2.plot(z_space, [0.0 for z in z_space] , 'r-',
                #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [Ohm/µW ],
                #  R_0={:5.3f} $\pm$ {:2.1e} [Ohm]'''.format(
                #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
                #           ))
                        )
            ax2.set_ylabel(r"Residuals [$\sigma$]")

        ax2.set_xlabel(x_label)

        ax2.grid(True)

        #make custom pruning of uppper tick (do not plot ticks in upper 10%)
        #so that ax2 tick does nto interfere with  ax1 tick
        ax2.locator_params(axis="y", min_n_ticks = 3
                            )
        y_loc = ax2.yaxis.get_majorticklocs()
        x_loc = ax1.xaxis.get_majorticklocs()
        #print("y_loc: ", y_loc)
        #print("y_loc[1:-2]: ", y_loc[1:-2])
        #print("ylim: ", ax2.get_ylim())
        y2_min, y2_max = ax2.get_ylim()
        y_loc = [y for y in y_loc if y2_min < y < y2_max - (y2_max -
                                                             y2_min)*0.1]
        #print("y_loc: ", y_loc)
        ax2.set_yticks(y_loc)
        # set  x lims:
        x1_min, x1_max = ax1.get_xlim()
        ax2.set_xticks(x_loc)
        ax2.set_xlim(ax1.get_xlim())

        # Delete xticks on 1st axis
        ax1.set_xticklabels([])

        # Optional table
        if table_string == None:
            pass
        else:
            plt.text(0.2,0.02,table_string,size=12,
            horizontalalignment='left',
            verticalalignment='bottom', transform=ax1.transAxes,
            backgroundcolor = "w")


        # fig.tight_layout()
        fig.subplots_adjust(left=0.2)

        format_im = 'png' #'pdf' or png
        dpi = 300
        plt.savefig(self.out_dir + plotname
                    + '.{}'.format(format_im),
                    format=format_im, dpi=dpi)
        # plt.show()
        # Reintroduce diiagnostics plot with chi sqared:
        ax1.errorbar(z_arr, P_arr,yerr = P_err_arr, fmt = ".",
                     color = "C0",
                    label = (r"data, $\chi^2_{red}$"
                            +" = {:2.3f} ".format(chi2_red)))
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend([handles[i] for i in  [0,2]],[labels[i] for i in  [0,2]],
            shadow=True, fontsize = 13)
        plt.savefig(self.out_dir + plotname
            + '_chi2.{}'.format(format_im),
            format=format_im, dpi=dpi, bbox_inches="tight")

        ax1.cla()
        fig.clf()
        plt.close()
        # Make sure there is always a copy of the analysis_run_dict.json saved
        #  with the plot
        # Save to file (saves additional  copy  to out_dir)
        self.save_json_run_dict(dict_path= self.run_dict_path, 
                            dict = self.run_dict)

    def plot_data(self,
                 P_arr, P_err_arr, z_arr,
            # Name plot
            plotname = "default_data",
             ):

        fig = plt.figure(0, figsize=(8,6.5), dpi =300)
        ax1=plt.gca()
        x_label = r"$z_{pos}$ [mm]"
        ### ax1
        ax1.errorbar(z_arr, P_arr,yerr = P_err_arr, fmt = ".",
                    label = (r"data"))


        ax1.set_ylabel(r"Power [µW]")
        ax1.set_xlabel(x_label)

        ax1.grid(True)
        ax1.legend(shadow=True, fontsize = 13)
        # ax1.tight_layout()


        fig.tight_layout()

        format_im = 'png' #'pdf' or png
        dpi = 300
        plt.savefig(self.out_dir + plotname
                    + '.{}'.format(format_im),
                    format=format_im, dpi=dpi)
        # plt.show()
        ax1.cla()
        fig.clf()
        plt.close()
        # Make sure there is always a copy of the analysis_run_dict.json saved
        #  with the plot
        # Save to file (saves additional  copy  to out_dir)
        self.save_json_run_dict(dict_path= self.run_dict_path, 
                            dict = self.run_dict)
        return
    
    def default_plot_data(self,plotname= "data_plot"):
        rd = self.run_dict
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
        self.plot_data(P_arr = P_arr, 
                       P_err_arr = P_err_arr, z_arr = z_arr,
                       plotname = plotname)
        return

    # HACK HACK HACK
    # make old pkl files loadable by loading them with old legacy code and 
    # extracting the  result_dict
    def legacy_load_pkl(self, pkl_path):
        # import sys
        # import os
        # # Put legacy script  in path
        # SCRIPT_DIR = os.path.dirname(
        #     "C:\\Users\\Christian\\Documents\\StudiumPhD\\"
        #     + "python\\Keysight-DMM-34461A\\analysis\\")
        
        # sys.path.append(os.path.dirname(SCRIPT_DIR))
        # import flow_on_off_cycle_analysis_2 as fca2
        # result_dict_unsorted = fca2.legacy_load_pkl(pkl_path)
        # return result_dict_unsorted


        # # Option 2: exectue legacy code in place:
        # path = ("C:\\Users\\Christian\\Documents\\StudiumPhD\\"
        #     + "python\\Keysight-DMM-34461A\\analysis\\"
        #     + "flow_on_off_cycle_analysis_2.py")
        # exec(open(path).read())

        # Option 3 just read locallyreexported result_dict.pkl
        extractor_dict_unsorted = load_dict(pkl_path)


        return extractor_dict_unsorted 



        


if __name__ == "__main__":
    pass
    