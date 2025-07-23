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
    flow = 1 * 4.478 * 10**17 #sccm
            ):
    #sccm = 4.478 * 10**17 # particles per second
    #flow = 1 * sccm

    #TODO Check if the detour is necessary or w ecould just start with the
    #  weighted interal
    integral = integrate_H_on_plane_1D(x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,)
    T_wire = np.array([350])  # Kelvins # SUPER ROUGH
    power_H2 = ((Cv(T) * (T) - Cv(T_wire) * (T_wire))  # [j/mol]
            * (1/N_A) # convert from particle nmber to number of mol
            * flow * integral # Number of particles hitting wire
                 )
    return power_H2

def ac_from_Abg(Abg, T,
    # l_eff = 7.96,
    # theta_max =  90 * degree,
    # z0 = 0,
    # y0=  35.17,
    # err = 0.01,
    flow = 1 * 4.478 * 10**17 #sccm
    ):
    # nf  = calc_norm_factor(l_eff = rd["fit_result"]["l_eff"])
    # eta_norm = 1.5 # Approimate effective normallization required for sim norm
    # d_wire = 5e-6 # wire thickness not included in P_fit (folded into A)
    # # given in m
    # y0 = 35.17e-3 # in m
    # # Transform from mm^2 to m^2 -> /1e6
    # #Integration was performed in mm -> need to include another factor of 1e-3
    # #ac_alt_lst.append(ac_alt* ((y0**2)/(nf*eta_norm*d_wire)))
    # ac_alt_lst.append(ac_alt*1e-3 * ((y0**2)/(nf*d_wire)))
    # # based on which factors are included in 
    # # .beamshape integrateH_on_plane_1D_etaW
    # # but not in simplified  P_int_penumbra_3par 
    # # where these are instead subsumed in A
    # # It seems eta_w does not need to be included, becasue it is actually in 
    # # P_int_penumbra_3par 


    T_wire = np.array([350])  # Kelvins # SUPER ROUGH
    ac = Abg/ (flow * ((Cv(T) * (T) - Cv(T_wire) * (T_wire))  # [j/mol]
            * (1/N_A)) # convert from particle nmber to number of mol)
                )
    return ac

def a_diss_from_A(A
                    , flow = 1 * 4.478 * 10**17 #sccm
                    , eta_rec = 1  # set to 1 for lower limit result
                    ):
    E_rec =7.1511 * 10**-19  # Joules
        # https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.121.013001
        # equivalent to 4.4634 eV
    # print(A, flow, E_rec, eta_rec)
    # print(type(A), type(flow))
    a_diss = (A)/ (flow * E_rec * eta_rec)
    return a_diss


# TODO Function for trasnforming from PID setpoint to (estimated) temperature
def TC_to_T_Hack(TC_val):
    # Simply digitzed the PID Thermocouple to T plot in discourse at 2 points
    # https://discourse.project8.org/t/mainz-habs-power-supply-tdk-lambda/291
    Ts = [317.8571428571429, 2185.714285714286]
    TCs = [51.111111111111114, 711.1111111111111]
    m = (Ts[1] - Ts[0])/(TCs[1] - TCs[0])
    b = Ts[0] - TCs[0] * m
    T_val = m * TC_val + b
    return T_val

# Callculate accomodation coefficient
# TODO include compensation for calibration "abberation" due to electric 
# heating
def calc_accomodation_coefficient(p_measured, T, x_lims = [-10, 10],
    z_lims = [-0.0025, 0.0025],
    l_eff = 7.96,
    theta_max =  90 * degree,
    z0 = 0,
    y0=  35.17,
    err = 0.01,
    flow = 1 * 4.478 * 10**17 #sccm
    ):


    # Effective total wire efficiency estimate
    # Before adjusting for beam distribution
    # Using the in vs outo-of-beam method (PROBLEMATIC HACK)
    # TODO account for th efac tthat "out-of-beam" isnt really
    power_H2 = predict_power_H2(T,
                x_lims = x_lims,
                z_lims = z_lims,
                l_eff = l_eff,
                theta_max =  theta_max,
                z0 = z0,
                y0=  y0,
                err = err,
                flow = flow # in  number of particles 1sccm =4.478 * 10**17 
                )
    effective_wire_eff  = p_measured/ power_H2
    # print("P_meas", p_measured * 1e6, "µW")
    # print("P_predict", power_H2 * 1e6, "µW")
    # print("effective_wire_eff", effective_wire_eff)




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
