import numpy as np
from .beamshape import integrate_H_on_plane_1D_etaW, integrate_H_on_plane_1D

# Define Heat capacity
# molar gas constant R = k*N_Avogadro
R = 8.31446261815324 # [J/mol K]
N_A = 6.02214076e23  # [1/mol]
k_B = R/N_A
eV_per_Joule = 6.241509e18 # [eV/J]

def Cp(T):
    #HACK in case T was not alreaadyan array
    T = np.asarray(T)
    # Initialize output  array
    output = np.zeros_like(T) 
    # for T less 298K use inapproproate extension of fit beyond 298
    # See Davidson 1927 A note on the specific heat of the hydrogen molecule
    # for the physical description
    # https://royalsocietypublishing.org/doi/pdf/10.1098/rspa.1927.0105
    cond = T < 298
    A = 33.066178
    B = -11.363417
    C = 11.432816
    D = -2.772874
    E = -0.158558
    t = T[cond]/1000
    output[cond] = (np.zeros_like(T[cond]) 
                    + A + B*t + C*t**2 + D*t**3 + E/t**2)
    ### For temps 298 to 1000 K
    cond = T > 298
    t = T[cond]/1000
    res = A + B*t + C*t**2 + D*t**3 + E/t**2
    output[cond] = res

    ### For temps 1000 to 2500 K
    A = 18.563083
    B = 12.257357
    C = -2.859786
    D = 0.268238
    E = 1.977990
    cond = T > 1000
    t = T[cond]/1000
    res = A + B*t + C*t**2 + D*t**3 + E/t**2
    output[cond] = res

    ### For temps 2500 to 6000 K
    A = 43.413560
    B = -4.293079
    C = 1.272428
    D = -0.096876
    E = -20.533862
    cond = T > 2500
    t = T[cond]/1000
    res = A + B*t + C*t**2 + D*t**3 + E/t**2
    output[cond] = res
    return output



def Cv(T):
    return Cp(T) - R # Output in [J/(mol * K)]

degree = np.pi/180
# Calculate the power expected to hit the wire 
def predict_power_H2(T, x_lims = [-10, 10],
    z_lims = [-0.0025, 0.0025],
    l_eff = 7.96,
    theta_max =  90 * degree,
    z0 = 0,
    y0=  35.17,
    err = 0.01,
            ):
    sccm = 4.478 * 10**17 # particles per second
    flow = 1 * sccm
    integral = integrate_H_on_plane_1D_etaW(x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,)
    T_wire = np.array([350])  # Kelvins # SUPER ROUGH
    power_H2 = (Cv(T) * (T - T_wire)  # [j/mol]
            * (1/N_A) # convert from particle nmber to number of mol
            * flow * integral # Number of particles hitting wire
                 )
    return power_H2

# TODO Function for trasnforming from PID setpoint to (estimated) temperature

# Callculate accomodation coefficient
# TODO include compensation for calibration "abberation" due to electric 
# heating
def calc_accomodation_coefficient(p_measured, T, x_lims = [-10, 10],
    z_lims = [-0.0025, 0.0025],
    l_eff = 7.96,
    theta_max =  90 * degree,
    z0 = 0,
    y0=  35.17,
    err = 0.01,):


    # Effective total wire efficiency estimate
    # Before adjusting for beam distribution
    # Using the in vs outo-of-beam method (PROBLEMATIC HACK)
    # TODO account for th efac tthat "out-of-beam" isnt really
    effective_wire_eff  = p_measured/ (predict_power_H2(T,
                x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,) )




    eta_wire_H_weighted = (integrate_H_on_plane_1D_etaW(x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,) 
                /
                integrate_H_on_plane_1D(x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,) )
    # Calculate wieghted wire efficiency
    # Dividing by wieghted wire seinsitivity, accounts for the drop in signal 
    # compared to predicitonn expected simply becasue the wire is less 
    # sensitive at the edges
    wire_eff_eta_compensated  = effective_wire_eff / eta_wire_H_weighted
    return wire_eff_eta_compensated 
