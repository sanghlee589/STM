import numpy as np
from numba import jit

@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def settling(Ds) :
    nu = 1.02*10**-6
    s = 2.65
    Rep = np.sqrt(9.81*(s-1)*Ds)*Ds/nu
    # Dietrich (1982) (manual 2-47a)
    b1 = 2.891394
    b2 = 0.95293
    b3 = 0.056835
    b4 = 0.002892
    b5 = 0.000245

    Rf_Dietrich = np.exp(-b1+b2*np.log(Rep)-b3*np.log(Rep)**2\
                         -b4*np.log(Rep)**3+b5*np.log(Rep)**4)
    Vs = Rf_Dietrich*np.sqrt(9.81*(s-1)*Ds) # m/s
  
    return Vs
    
@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def crit_shear_stress(Ds) :
    # particle and water density are set to 2.65 and 1, respectively.
    
    # kinematic viscosity of water (nu) (m2/s)
    # if temperature available, this can be estimated as follows.
    # nu = 1.79*10**-6/(1+0.03368*T+0.000221*T**2)
    nu = 1.02*10**-6 # at 20 oC
    
    # Particle Reynolds number
    R_ep = (1.65*9.81*Ds)**0.5*Ds/nu
    # dimensionlee critical shear stress
    tau_s_cr = 0.5*(0.22*R_ep**-0.6+0.06*np.exp(-17.77*R_ep**-0.6))
    # critical shear stress
    tau_cr = 1.65*9.81*Ds*tau_s_cr

    
    return tau_cr

@jit(nopython=True) # Set "nopython" mode for best performance, equivalent to @njit
def deposit(Vs, c, H, tau_b,tau_cr,threshold) : #settling v, concentration, water depth, bed and criticl shear stress
    
    dep = np.zeros(c.shape)
                
    for i in range(c.shape[0]) :
        for j in range(c.shape[1]) :
            if H[i,j] < threshold : # if no water on surface, all sediment deposited
                dep[i,j] = c[i,j]

    if np.max(tau_b) > tau_cr : # no deposition across the domain
        pass      
    else : # need to estimate deposition cell by cell
        for i in range(c.shape[0]) :
            for j in range(c.shape[1]) :
                if H[i,j] < threshold : # if no water on surface, all sediment deposited
                    pass              
                else :
                    dep[i,j] = Vs*c[i,j]/H[i,j]*(1-tau_b[i,j]/tau_cr)
                    if dep[i,j] > c[i,j]:
                        dep[i,j] = c[i,j]
    
    return dep
    

    