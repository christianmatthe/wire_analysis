import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy.interpolate import interp1d

#plot Options
import matplotlib as mpl
font = {#'family' : 'normal','weight' : 'bold',
        'size'   : 16
        #,'serif':['Helvetica']
        }
mpl.rc('font', **font)
mpl.rcParams['figure.dpi'] = 400


# make a class for fitting l_eff Tschersich model to beam shape
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
    "

    Parameters
    ----------
    test :  `bool`
        test parameter
    
    """
    def __init__(self,
                test = True,
                ):
        #define parameters
        self.test = test



        #define constants
        degree = np.pi/180 # convert form rad to degree
        return


    # Auxiliary functions required for profile shape according to
    # Tschcersich
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

if __name__ == "__main__":
    print(Beamfit().test)
