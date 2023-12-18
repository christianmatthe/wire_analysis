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
    def save_fit_results(self, popt_abs, pcov_abs):
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
            x_label = r"central angle [deg]"
            z_arr = np.arctan((z_arr - z0) / y0) / self.degree
            z_space = np.arctan((z_space - z0) / y0) / self.degree
            ax1.axvline(theta_max, label=r"$\theta_{max} = $"
                        + f"{theta_max:.2f} [deg]", 
                        color = "gray", alpha = 0.6, ls = "--")
            ax2.axvline(theta_max, label=r"$\theta_{max}$", 
                        color = "gray", alpha = 0.6, ls = "--")
        else:
            x_label = r"$z_{pos}$ [mm]"
        ### ax1
        ax1.errorbar(z_arr, P_arr,yerr = P_err_arr, fmt = ".",
                    label = r"data, $\chi^2_{red}$"+" = {:2.3f} ".format(chi2_red))

        # Plot Fit
        #xdata=p_arr
        if l_eff_str == None:
            ax1.plot(z_space, P_space_eye, "r-", label = r"$P_{fit}$")
        else:
            ax1.plot(z_space, P_space_eye, "r-", label = r"$P_{fit}$" 
                    + r", $l_{eff}=$" + l_eff_str
            #  label=(('''fit: m={:5.4f} $\pm$ {:2.1e} [$\Omega$/µW],
            #  R$_0$={:5.3f} $\pm$ {:2.1e} [$\Omega$]'''.format(
            #      popt[0], np.sqrt(pcov[0,0]), popt[1], np.sqrt(pcov[1,1])) 
            #           ))
                    )
        if P_space_compare is None:
            pass
        else:
            ax1.plot(z_space, P_space_compare, "-",color = "C1", 
                    label = P_compare_label
                    )


        ax1.set_ylabel(r"power [µW]")
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


        fig.tight_layout()
        fig.subplots_adjust(left=0.2)

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
    