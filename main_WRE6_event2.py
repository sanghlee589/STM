# =============================================================================
'''
 Sediment transport in overland flow
 (Running 2D advection-dispersion equation with source and sink terms)
 @author: Sanghyun Lee

 - 3D Array configuration [time, x, y]
 c1, c2, c3, c4, c5 : 3D array for sediment concentration
 u and v : 3D array for flow fields in x and y directions
 D : 3D array for diffusion coefficients including 2 boundary layers


< Needed inputs for running WEPP-HE to obtain sediment source >
 - static input (2D array [x y]) : 
    vfs (fractin of very fine sand in the surface soil)
    sand (fraction of sand)
    silt (fraction of silt)
    clay (fraction of clay)
    orgmat (organic matter)
    slope (m/m)

 - dynamic input (3D array [time x y]) :
    H = surface runoff depth (m)
    Q_P = peak runoff rate (m/s)
    Q_T = effective runoff duration (s)
    P_I = effective rainfall intensity (m/s)
    P_T = effective rainfall duration (s)

 Adjusted parameters for interrill and rill erodibility, and critical shear stress
should be included as well (3D array).

< Needed inputs for advection-diffusion equaiton for sediment transport >
    H (m) as in the source model, 
    slope (m/m) as in the source model,
    flow fields in x and y (m/s) in 3D
    D (diffusion coefficients) will be calculated in the code.
'''
# =============================================================================

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
from mikeio import Dfs2 # to read MIKE SHE outputs
import timeit

start_t = timeit.default_timer()

# functions below are made for solving advection-diffusion equation
from bc import two_boundary 
from solver import D_C, diffusion, LW, advection_x, advection_y
from deposition import settling, crit_shear_stress, deposit
from source import WEPP_HE

path =r"F:\04_sediment_model_ch3\WRE watersheds\WRE6\sediment_transport_model"
MIKE_SHE = r"F:\04_sediment_model_ch3\WRE watersheds\WRE6\MIKE_SHE\FR6.she - Result Files"
os.chdir(path)

# time domain (for sediment model)
start = datetime(1994,4,6,0)
end = datetime(1994,5,21,0)
time_step = 'days' # days or hours

# Spatial domain
x = 80 # entire length in meter
y = 200 # entire length in meter
dx = dy = 20 # grid size in meter
nx = int( x / dx ) 
ny = int( y / dy ) 

if time_step == 'days':
    delta = timedelta(days=1)
    nt = (end-start).days + 1
    dt = 60*60*24 # in seconds
    
elif time_step == 'hours':
    delta = timedelta(hours=1)   
    nt = int((end-start).total_seconds()/60/60) + 1
    dt = 60*60 # in seconds

# Slope
S = np.zeros((nx,ny)) # slope
S[:,:] = 2.9 / 100 # slope

# =============================================================================
# for event-based simulations (for case start date does not match with the events)
# =============================================================================
sim_begin = datetime(1993,1,1)
shift = (start-sim_begin).days

# =============================================================================
# Stand-alone WEPP model inputs
# =============================================================================
# Particle fraction for topsoil
sand= np.round(np.loadtxt("SAND.txt") / 100 ,3)
silt = np.round(np.loadtxt("SILT.txt") / 100 ,3)
clay = np.round(np.loadtxt("CLAY.txt") / 100 ,3)
#orgmat = 1.724*orgc # if only organic carbon content is available (from WEPP document)
orgmat = np.round(np.loadtxt("orgmat.txt") / 100,3)
vfs = 0.2 # fractin of very fine sand in the surface soil 

# Time varying variables
RCF = np.loadtxt("RCF.txt") # residue cover fraction (0-1)
canhgt = np.loadtxt("canhgt.txt") # canopy height of plant (m)
RRC = np.loadtxt("RRC.txt") / 1000 # soil random roughness (m)

# KIADJ = np.loadtxt('ki_factor.txt') # adjusted factors for interrill erodibility
# KRADJ = np.loadtxt('kr_factor.txt') # adjusted factors for rill erodibility
# SHCRTADJ = np.loadtxt('crit_factor.txt') # adjusted factors for critical shear stress
KIADJ = 0.07 * np.ones(len(RCF))
KRADJ = 0.07 * np.ones(len(RCF))
SHCRTADJ = np.ones(len(RCF))
detention = 0.0002 # detention storage for deposition term [m]

# total rill frinction factor
FRCTRL = 1.11 + 4.5*RCF**1.55 + (canhgt/np.max(canhgt))

# adjust start date
RCF = RCF[shift:]
canhgt = canhgt[shift:]
KIADJ = KIADJ[shift:]
KRADJ = KRADJ[shift:]
SHCRTADJ = SHCRTADJ[shift:]
RRC = RRC[shift:]
FRCTRL = FRCTRL[shift:]

# =============================================================================
# Sediment transport model inputs
# =============================================================================
# Initialise input arrays
D = np.zeros((nt,nx+4,ny+4)) # diffusion coefficients (m2/s)
u = np.zeros((nt,nx,ny)) # velocities in x direction (m/s)
v = np.zeros((nt,nx,ny)) # velocities in y direction (m/s)
H = np.zeros((nt,nx,ny)) # water depth (m)

# =============================================================================
# Import hydrologic outputs from MIKE SHE (including outputs for overland flow)
# =============================================================================
dfs_Q = Dfs2(MIKE_SHE + "\FR6_overland.dfs2")
ds = dfs_Q.read()
H = ds.data[0] # depth of overland water <Water Depth> (meter)
u = ds.data[1] # overland flow in x-direction (m/sec)
v = ds.data[2] # overland flow in y-direction (m/sec)

# remove warm-up years (first 3 years (1990-1992))+shift and boundary data
H = H[1096+shift:,1:-1, 2:-1]
u = u[1096+shift:,1:-1, 2:-1]
v = v[1096+shift:,1:-1, 2:-1]

# =============================================================================
# Calculating sediment source using stand-alone WEPP
# =============================================================================
# source array by particle size
source_1 = np.zeros((nt,nx,ny)) # particle size: 0.002mm
source_2 = np.zeros((nt,nx,ny)) # particle size: 0.010mm
source_3 = np.zeros((nt,nx,ny)) # particle size: 0.030mm
source_4 = np.zeros((nt,nx,ny)) # particle size: 0.300mm
source_5 = np.zeros((nt,nx,ny)) # particle size: 0.200mm

# Initialise input arrays for stand-alone WEPP (should be imported as inputs)
Q_P = np.zeros((nt,nx,ny))

# hydrologic results (from hydrologic model's output)
P_I = np.loadtxt("rain_intensity.txt")  # effective rainfall intensity (m/s)
P_T = np.loadtxt("rain_duration.txt") # effective rainfall duration (s)
P_I = P_I[shift:]
P_T = P_T[shift:]
Q_T = P_T # effective runoff duration (s), assuming runoff duration=rainfall duration

# estimating sub-daily peak discharge (modified rational method equation)
# https://www.scielo.br/j/rbcs/a/crVVDQLbsdCqGwNfY3K9zvz/?format=pdf&lang=en
# Area: 0.02km2, time of concentration: 0.08 hr (based on TR-55)
# for i in range(nt):
#     Q_P[i,:,:] = H[i,:,:] * 1000 * 0.02 / 3.6 / 0.08 / 80 / 200 # m/s

# simul_time = start
# for k in range(nt):
#     for i in range(nx):
#         for j in range(ny):
#             sediment = WEPP_HE(dx, H[k,i,j], Q_P[k,i,j], Q_T[k], \
#                                 P_I[k], P_T[k], S[i,j],    \
#                                 vfs, sand[i,j], silt[i,j], clay[i,j],\
#                                 orgmat[i,j], RCF[k], canhgt[k], KIADJ[k], KRADJ[k], \
#                                 SHCRTADJ[k], RRC[k], FRCTRL[k])
#             source_1[k,i,j] = sediment[0] / 3600 / 24 # kg/m2/day to kg/m2/s
#             source_2[k,i,j] = sediment[1] / 3600 / 24  # kg/m2/day to kg/m2/s
#             source_3[k,i,j] = sediment[2] / 3600 / 24  # kg/m2/day to kg/m2/s
#             source_4[k,i,j] = sediment[3] / 3600 / 24  # kg/m2/day to kg/m2/s
#             source_5[k,i,j] = sediment[4] / 3600 / 24  # kg/m2/day to kg/m2/s
 
    
#     print('Sediment production simulation is done for '\
#           +simul_time.strftime("%Y-%m-%d-%H"))
#     simul_time += delta

# np.save('source_1_event2.npy', source_1)
# np.save('source_2_event2.npy', source_2)
# np.save('source_3_event2.npy', source_3)
# np.save('source_4_event2.npy', source_4)
# np.save('source_5_event2.npy', source_5)

source_1 = np.load('source_1_event2.npy')
source_2 = np.load('source_2_event2.npy')
source_3 = np.load('source_3_event2.npy')
source_4 = np.load('source_4_event2.npy')
source_5 = np.load('source_5_event2.npy')



# =============================================================================
# Calculating Diffusion coefficients
# =============================================================================
# bed shear stress for normal flow
tau_b = 1000*9.81*H[:nt,:,:]*S[:,:]
dep = np.zeros([nt,nx,ny])

# particle diameter 
# clay, silt, small and large aggregate, sand
# (0.002, 0.01, 0.03, 0.3, 0.2) mm
for Ds in [0.001/1000, 0.01/1000, 0.03/1000, 0.3/1000, 0.2/1000]: # this could be done in parallel
    s0 = np.zeros((nt,nx+4,ny+4)) # concentration (kg/m2)    
    Vs = settling(Ds) # setting velocity (m/s)
    
    # criticl shear stress
    tau_cr = crit_shear_stress(Ds)
    # Diffusion coefficient
    D[:,2:-2,2:-2] = D_C(tau_b, H[:nt,:,:], S) 

    # Stability criteria (Courant number)
    stability_criteria = 1/4
    
    # set the initial boundary condition (zero-gradient)
    s0 = two_boundary(s0,0)
    D = two_boundary(D, 0)

# =============================================================================
# Simulations 
# =============================================================================
    simul_time = start + delta
    n = 0
    original_dt = dt
    
    while simul_time <= end:
       
        # source (from stand-alone WEPP)
        if Ds == 0.001/1000:
            s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] + source_1[n,:,:]
        elif Ds == 0.01/1000:
            s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] + source_2[n,:,:]
        elif Ds == 0.03/1000:
            s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] + source_3[n,:,:]
        elif Ds == 0.3/1000:
            s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] + source_4[n,:,:]
        elif Ds == 0.2/1000:
            s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] + source_5[n,:,:]

        # sink (deposition)
        dep[n,:,:] = deposit(Vs, s0[n,2:-2, 2:-2], H[n,:,:], tau_b[n,:,:], tau_cr,detention)
        s0[n,2:-2, 2:-2] = s0[n,2:-2, 2:-2] - dep[n,:,:]
        # if components were below zero due to deposition, make it zero
        s0[n,2:-2, 2:-2][s0[n,2:-2, 2:-2] < 0] = 0
        
        # update boundary
        s0 = two_boundary(s0,n)
    
        dt = original_dt
        # Stability check (Courant number)
        check = max(np.max(D[n,:,:])*dt/dx**2 + np.max(D[n,:,:])*dt/dy**2, \
                    np.max(abs(u[n,:,:]))*dt/dx + np.max(abs(v[n,:,:]))*dt/dy)
        if check > stability_criteria :
            original_dt = dt
            while check > stability_criteria :
                dt = dt / 2
                check = max(np.max(D[n,:,:])*dt/dx**2 + np.max(D[n,:,:])*dt/dy**2, \
                            np.max(abs(u[n,:,:]))*dt/dx + np.max(abs(v[n,:,:]))*dt/dy)
            mu = int(original_dt / dt)
            print('    time step reduced to '+str(mu)+' times lesser than original setting')
            # additional runs using smaller dt until next time step
            s_mu = np.zeros((mu+1,nx+4,ny+4))
            s_mu[0,:,:] = s0[n,:,:]
            for w in range(mu) :
                s_mu[w,2:-2, 2:-2] = advection_x(s_mu, w, n, u, nx, dx, dt)
                s_mu[w+1,2:-2, 2:-2] = advection_y(s_mu, w, n, v, ny, dy, dt)
                s_mu[w+1,:,:] = diffusion(s_mu[w+1,:,:],D[n,:,:],dx,dy,dt)
                s_mu = two_boundary(s_mu,w+1)
            # update domain for the next time step                
            s0[n+1,:,:] = s_mu[w+1,:,:]    
                        
        else:
            # simulate with the regular time step
            s0[n,2:-2, 2:-2] = advection_x(s0, n, n, u, nx, dx, dt)
            s0[n+1,2:-2, 2:-2] = advection_y(s0, n, n, v, ny, dy, dt)
            s0[n+1,:,:] = diffusion(s0[n+1,:,:],D[n,:,:],dx,dy,dt)
            s0 = two_boundary(s0,n+1)
            
      
        D = two_boundary(D, n+1)
        
        
        print('Sediment transport simulation for D= '+str(Ds*1000)+\
              ' mm is done for '+simul_time.strftime("%Y-%m-%d-%H"))
        
        simul_time += delta
        n += 1 
    # sediment concentration results for each particle class
    if Ds == 0.001/1000:
        c1 = s0[:,2:-2, 2:-2] # sediment concentration for D = 0.001 mm (kg/m2/s)
        dep_1 = dep
    elif Ds == 0.01/1000:
        c2 = s0[:,2:-2, 2:-2] # sediment concentration for D = 0.01 mm (kg/m2/s)
        dep_2 = dep
    elif Ds == 0.03/1000:
        c3 = s0[:,2:-2, 2:-2] # sediment concentration for D = 0.03 mm (kg/m2/s)
        dep_3 = dep
    elif Ds == 0.3/1000:
        c4 = s0[:,2:-2, 2:-2] # sediment concentration for D = 0.3 mm (kg/m2/s)
        dep_4 = dep
    elif Ds == 0.2/1000:
        c5 = s0[:,2:-2, 2:-2] # sediment concentration for D = 0.2 mm (kg/m2/s)
        dep_5 = dep

# sediment load at the outlet
load_1 = c1 * dx * dy * 3600 * 24 # (kg)
load_2 = c2 * dx * dy * 3600 * 24 # (kg)
load_3 = c3 * dx * dy * 3600 * 24 # (kg)
load_4 = c4 * dx * dy * 3600 * 24 # (kg)
load_5 = c5 * dx * dy * 3600 * 24 # (kg)

stop_t = timeit.default_timer()
print('Time: ', stop_t - start_t)  

# calculate total load at the outlet
sed_simul = np.zeros(nt)
for i in range(nt):
    # watershed area = 16000 m2
    sed_simul[i] = load_1[i,:,0].sum() + load_2[i,:,0].sum() +\
        load_3[i,:,0].sum() +load_4[i,:,0].sum() +load_5[i,:,0].sum()

# sediment delivery ratio (SDR)
s_y = np.sum(sed_simul) # total sediment yield
g_e = ( np.sum(source_1) + np.sum(source_2) +np.sum(source_3)\
        +np.sum(source_4) +np.sum(source_5) ) * dx * dy * 3600 * 24 # gross erosion
SDR = s_y / g_e

# =============================================================================
# Comparison with sediment observations
# =============================================================================
os.chdir(path)
sed_observed = np.loadtxt('Sed_observed.txt') # 1993 to 1995
sed_observed = sed_observed[shift:(shift+nt)]
# cumulative sediment (mm)
sed_cum_simul = np.cumsum(sed_simul)
sed_cum_observed = np.cumsum(sed_observed)

# Nash-Sutcliffe coefficient
NSE = 1 - np.sum((sed_observed-sed_simul)**2) / np.sum((sed_observed-np.mean(sed_observed))**2)
print('NSE = '+str(np.round(NSE,2)))

# PBIAS
p_bias = np.sum(sed_simul-sed_observed) / np.sum(sed_observed) * 100
print('PBIAS ='+str(np.round(p_bias,2))+'%')
print('Sediment delivery ratio: '+str(np.round(SDR,2)))

# =============================================================================
# Plot the results
# =============================================================================
# sediment yield at outlet
dti = pd.date_range(start, periods=nt, freq="D")
sand_plt = np.zeros(nt)
lg_plt = np.zeros(nt)
sa_plt = np.zeros(nt)
silt_plt = np.zeros(nt)
clay_plt = np.zeros(nt)
# sand
for i in range(nt):
    sand_plt[i] = load_5[i,:,0].sum()
# large aggregate
for i in range(nt):
    lg_plt[i] = load_4[i,:,0].sum() +load_5[i,:,0].sum()
# small aggregate
for i in range(nt):
    sa_plt[i] = load_3[i,:,0].sum() +load_4[i,:,0].sum() +load_5[i,:,0].sum()
# silt
for i in range(nt):
    silt_plt[i] = load_2[i,:,0].sum() +load_3[i,:,0].sum() +load_4[i,:,0].sum() +load_5[i,:,0].sum()

# clay
for i in range(nt):
    clay_plt[i] = load_1[i,:,0].sum() + load_2[i,:,0].sum() +\
        load_3[i,:,0].sum() +load_4[i,:,0].sum() +load_5[i,:,0].sum()

dti = dti.strftime("%m/%d/%y")
interval = 10
plt.rcParams["figure.figsize"] = (10,8)
plt.rcParams.update({'font.size': 12})
plt.figure(2)
plt.subplot(311)
Q_observed = np.loadtxt('observed_Q_1993_1995.txt')
Q_observed = Q_observed[shift:shift+nt]
p = np.loadtxt('P_1993_1995.txt')
p = p[shift:shift+nt]
Q_simulated = np.zeros(sed_observed.shape)
q_temp = np.sqrt(u**2 + v**2)
# Convert to runoff depth (mm)
for i in range(nt):
    # watershed area = 16000 m2
    Q_simulated[i] = abs(q_temp[i,:,0].sum() * 3600 * 24 / 16) # m3/s to mm
# Nash-Sutcliffe coefficient for runoff
NSE_Q = 1 - np.sum((Q_observed-Q_simulated)**2) / np.sum((Q_observed-np.mean(Q_observed))**2)

plt.plot(dti,Q_observed, '-o', markersize=3,label='Measured Q')
plt.plot(dti,Q_simulated, label='Simulated Q')
plt.plot(dti,p,label='P',alpha = 0.3)
plt.ylabel('Runoff [mm/day]')
plt.xlabel('(a)')
ax = plt.gca()
ax.set_xticks(dti[::interval])
ax.set_xticklabels(dti[::interval])
ax.set_xticklabels([])
plt.text(dti[10],70, 'NSE = '+str(np.round(NSE_Q,2)),fontsize=12)
plt.legend()

plt.subplot(312)
plt.fill_between(dti,0, sand_plt, label='Sand')
plt.fill_between(dti,sand_plt, lg_plt, label='Large aggregate')
plt.fill_between(dti,lg_plt, sa_plt, label='Small aggregate')
plt.fill_between(dti,sa_plt, silt_plt, label='Silt')
plt.fill_between(dti,silt_plt, clay_plt, label='Clay')
plt.plot(dti,sed_observed, 'black',label='Measured')
#plt.plot(dti,sed_simul, label='simulated')
plt.ylabel('Sediment [kg/day]')
plt.xlabel('(b)')
ax = plt.gca()
ax.set_xticks(dti[::interval])
ax.set_xticklabels(dti[::interval])
ax.set_xticklabels([])
plt.text(dti[10],350, 'NSE = '+str(np.round(NSE,2)),fontsize=12)
plt.legend()

plt.subplot(313)
sand_plt_cum = np.cumsum(sand_plt)
lg_plt_cum = np.cumsum(lg_plt)
sa_plt_cum = np.cumsum(sa_plt)
silt_plt_cum = np.cumsum(silt_plt)
clay_plt_cum = np.cumsum(clay_plt)

plt.fill_between(dti,0, sand_plt_cum, label='Sand')
plt.fill_between(dti,sand_plt_cum, lg_plt_cum, label='Large aggregate')
plt.fill_between(dti,lg_plt_cum, sa_plt_cum, label='Small aggregate')
plt.fill_between(dti,sa_plt_cum, silt_plt_cum, label='Silt')
plt.fill_between(dti,silt_plt_cum, clay_plt_cum, label='Clay')
plt.plot(dti,sed_cum_observed, 'black',label='Measured')
plt.ylabel('Sediment [kg]')
plt.xlabel('(c)')
ax = plt.gca()
ax.set_xticks(dti[::interval])
ax.set_xticklabels(dti[::interval])
plt.text(dti[10],500, 'PBIAS = '+str(np.round(p_bias,2))+'%',fontsize=12)
plt.tight_layout()
plt.savefig('event2_WRE6.png',dpi=300)
plt.show()


# sediment production [kg]
s_p = (source_1 +source_2+source_3+source_4+source_5) * dx * dy * 3600 * 24 

print('')
print('fraction of clay            to total source: '\
      +str(np.round(np.sum(source_1)*24*3600*dx*dy / np.sum(s_p),2)))
print('fraction of silt            to total source: '\
      +str(np.round(np.sum(source_2)*24*3600*dx*dy / np.sum(s_p),2)))
print('fraction of small aggregate to total source: '\
      +str(np.round(np.sum(source_3)*24*3600*dx*dy / np.sum(s_p),2)))
print('fraction of large aggregate to total source: '\
      +str(np.round(np.sum(source_4)*24*3600*dx*dy / np.sum(s_p),2)))
print('fraction of sand            to total source: '\
      +str(np.round(np.sum(source_5)*24*3600*dx*dy / np.sum(s_p),2)))
print('')
print('fraction of clay            to total yield at outlet: '\
      +str(np.round(np.sum(load_1[:,:,0].sum()) / np.sum(sed_simul),2)))
print('fraction of silt            to total yield at outlet: '\
      +str(np.round(np.sum(load_2[:,:,0].sum()) / np.sum(sed_simul),2)))
print('fraction of small aggregate to total yield at outlet: '\
      +str(np.round(np.sum(load_3[:,:,0].sum()) / np.sum(sed_simul),2)))
print('fraction of large aggregate to total yield at outlet: '\
      +str(np.round(np.sum(load_4[:,:,0].sum()) / np.sum(sed_simul),2)))
print('fraction of sand            to total yield at outlet: '\
      +str(np.round(np.sum(load_5[:,:,0].sum()) / np.sum(sed_simul),2)))