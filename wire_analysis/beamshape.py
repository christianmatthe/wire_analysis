# To contain fucntions for modeling the Tschersich beam shape
# Based on accomodation_coeffficient_H2_2023-12-06
import numpy as np
from scipy import integrate
from .beamfit import Beamfit
# class Beamshape():
#     def __init__(self,
#             ):
#             # initiate contstants
#             self.degree = np.pi/180
#             return
# Initially do not set this up as a class to "save time"

# Initially I will define the function  describing the emission of an atomic 
# Hydrogen beam form a hot capillary according to 
# https://aip.scitation.org/doi/pdf/10.1063/1.368619 
# (K. G. Tschersich; V. von Bonin (1998))
# remember to input values in radians
degree = np.pi/180
# remember to input values in radians
def beta(theta, l_eff):

    output = np.zeros_like(theta) 
    cond = np.abs(l_eff * np.tan(theta)) < 1
    output[cond] =  np.arccos(l_eff * np.tan(theta[cond]))
    return output


def U(theta, l_eff):
    # Move conditional to beta only
    return (2*beta(theta, l_eff)-np.sin(2*beta(theta, l_eff)))/np.pi

def V(theta, l_eff):
    # Move conditional to beta only
    return np.sin(beta(theta, l_eff))**3

def jd(theta, l_eff):
    return np.cos(theta) * U(theta, l_eff)

def jw(theta, l_eff):
    # result = (
    # (4/(3*np.pi))*(1-1/(2*l_eff + 1)) * (1/l_eff) 
    # * (np.cos(theta)**2 / np.sin(theta)) 
    # * (1-V(theta, l_eff))
    # + (1/(2*l_eff + 1))*np.cos(theta) * (1-U(theta, l_eff))
    # )
    jw_lambda = lambda theta, l_eff: (
        (4/(3*np.pi))*(1-1/(2*l_eff + 1)) * (1/l_eff) 
        * (np.cos(theta)**2 / np.sin(theta)) 
        * (1-V(theta, l_eff))
        + (1/(2*l_eff + 1))*np.cos(theta) * (1-U(theta, l_eff))
        )

    theta = np.array(np.abs(theta))
    cond = (theta != 0)
    result = np.piecewise(theta, 
        [theta == 0, cond],
        [
            0,jw_lambda(theta[cond], l_eff),
        ]
            )
    return result

def j(theta, l_eff):
    return jd(theta, l_eff) + jw(theta, l_eff)

#Chop at a maximum angle
def H_profile(theta, l_eff, theta_max):
    # Try including the condition again when inputing theta into j as Florian
    # showed me
    # cond = (theta != 0) & (theta < theta_max) & not (theta >= theta_max) 

    # Theta needs to be connverted to array if it isn't already
    ## also convert to positive anngles onyl while we are at it

    theta = np.array(np.abs(theta))
    cond = (theta < theta_max) & (theta != 0)
    result = np.piecewise(theta, 
        [theta == 0, cond, theta >= theta_max],
        [
            1,j(theta[cond], l_eff), 0   
        ]
            )
    return result

#### Normalize H_profile to 1 over solid angle
# H_profile is to be a pdf. The probability of finding a particle in the 
# entire half sphere of solid angle must be 1.
# This should translate to a probability of 1 of finding it in the projected 
# plane

# Define integrals over angle in order to efficiently integrate over the 
# Half-sphere in order to normalize with these
def integrate_H_angles(theta_lim = [0,np.pi/2], l_eff = 7.96, y0 = 35.17,
                       theta_max = 90 * degree, norm_factor = 1):
    integrant = lambda theta, phi: (norm_factor
                        * H_profile(theta, l_eff, theta_max) 
                             * np.sin(theta)
                        )
    result, err = integrate.dblquad(integrant,
                                0,  2* np.pi, # phi_limits
                                theta_lim[0], theta_lim[1] # theta_lims
                                )
    return result

def integrate_H_angles_1D(theta_lim = [0,np.pi/2], l_eff = 7.96, y0 = 35.17,
                       theta_max = 90 * degree, norm_factor = 1):
    # Integral over phi  from [0, 2 * np.pi]
    integrant = lambda theta: (norm_factor * 2 * np.pi
                        * H_profile(theta, l_eff, theta_max) 
                             * np.sin(theta)
                        )
    result, err = integrate.quad(integrant,
                                theta_lim[0], theta_lim[1] # theta_lims
                                )
    return result
    
def calc_norm_factor(l_eff
                     ):
    return 1/integrate_H_angles_1D(l_eff = l_eff)


def integrate_H_on_plane(x_lims = [-10,10], z_lims = [-2.5e-3,2.5e-3],
                        l_eff = 7.96, theta_max = 90 * degree, 
                        z0 = 0, y0  = 35.17,
                         err = 1e-2, norm_factor = None ):
    if norm_factor == None:
        # per default calculate noormalization factor
        # The function is much faster if norm_factor is provided
        norm_factor = calc_norm_factor(l_eff = l_eff)
          
    
    def theta(x,z):
            return np.arctan(np.sqrt(x ** 2 + (z - z0) ** 2) / y0)
    integrant = lambda x,z: ( norm_factor
                            * H_profile(theta(x,z), l_eff, theta_max)
                            * 1/(y0**2 * (1/np.cos(theta(x,z))**3))
                            # from solid angle to area on plane
                                )
    result, err = integrate.dblquad(integrant,
                                    z_lims[0], z_lims[1], # z boundaries
                                    x_lims[0], x_lims[1], # x boundaries
                                    epsabs=err,epsrel= err,
                                    )
    return result

def integrate_H_on_plane_1D(x_lims = [-10,10], z_lims = [-2.5e-3,2.5e-3],
                        l_eff = 7.96, theta_max = 90 * degree, 
                        z0 = 0, y0  = 35.17,
                        err = 1e-2, norm_factor = None):
    """
    This function serves to integrate the H distribution along a thin rectangle
    i.e. a projected wire. 
    DO NOT ENTER LARGE z_lims. Keep them well below mm

    The simplification to a 1D integral greatly speeds up the integration
    """
    if norm_factor == None:
        # per default calculate noormalization factor
        # The function is much faster if norm_factor is provided
        norm_factor = calc_norm_factor(l_eff = l_eff)
    
    z_center  = (z_lims[1] + z_lims[0])/2
    z_width = z_lims[1] - z_lims[0]
    def theta(x):
            return np.arctan(np.sqrt(x ** 2 + (z_center - z0) ** 2) / y0)
    integrant = lambda x: ( norm_factor
                           * z_width
                            * H_profile(theta(x), l_eff, theta_max)
                            * 1/(y0**2 * (1/np.cos(theta(x))**3))
                            # from solid angle to area on plane
                                )
    result, err = integrate.quad(integrant,
                                    x_lims[0], x_lims[1], # x boundaries
                                    epsabs=err,epsrel= err,
                                    )
    return result
    

def integrate_H_on_plane_1D_etaW(x_lims = [-10,10], z_lims = [-2.5e-3,2.5e-3],
                        l_eff = 7.96, theta_max = 90 * degree, 
                        z0 = 0, y0  = 35.17,
                        err = 1e-2, norm_factor = None):
    #HACK this path does not hold up on other systems
    eta_wire = Beamfit(run_dict_path=(  
        "C:\\Users\\Christian\\Documents\\StudiumPhD\\python\\wire_analysis"
        + "\\scripts\\2023-12-18_no_cracking_H2\\run_dicts\\1sccm_390TC.json"
                       )).eta_wire_default
    """
    This function serves to integrate the H distribution along a thin rectangle
    i.e. a projected wir, while accounting for wire sensitivity.
    DO NOT ENTER LARGE z_lims. Keep them well below mm

    The simplification to a 1D integral greatly speeds up the integration
    """
    if norm_factor == None:
        # per default calculate noormalization factor
        # The function is much faster if norm_factor is provided
        norm_factor = calc_norm_factor(l_eff = l_eff)
    
    z_center  = (z_lims[1] + z_lims[0])/2
    z_width = z_lims[1] - z_lims[0]
    def theta(x):
            return np.arctan(np.sqrt(x ** 2 + (z_center - z0) ** 2) / y0)
    integrant = lambda x: ( norm_factor
                           * z_width
                            * H_profile(theta(x), l_eff, theta_max)
                            * 1/(y0**2 * (1/np.cos(theta(x))**3))
                            # from solid angle to area on plane
                            * eta_wire(x) # weight by wire sensitivity
                                )
    result, err = integrate.quad(integrant,
                                    x_lims[0], x_lims[1], # x boundaries
                                    epsabs=err,epsrel= err,
                                    )
    return result