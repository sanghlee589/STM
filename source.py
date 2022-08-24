# =============================================================================
# Running stand-alone WEPP model 
# The results provide source term (sediment production) in ADE
# Created by Sanghyun Lee
# =============================================================================

import numpy as np
import os
import re
import subprocess   

# path that contains all needed files to run stand-alone WEPP code
path = r"F:\04_sediment_model_ch3\WRE watersheds\WRE6\sediment_transport_model\Erosion"

def WEPP_HE(length,Q, Q_P, Q_T, P_I, P_T, slope, \
            vfs, sand, silt, clay, orgmat, RCF, canhgt, \
            KIADJ, KRADJ, SHCRTADJ, RRC, FRCTRL) :
    
# =============================================================================
# Assign input data
# =============================================================================
    width = length # for square grid cells (m)
    
# =============================================================================
# Make hydrologic pass file
# =============================================================================
    os.chdir(path)
    h = open("Hyd2er.PAS",'w')
    h.write("   "+"{:.4e}".format(Q)+"  "+"{:.4e}".format(Q_P)+"   "+str(Q_T)+\
            "      "+"{:.4e}".format(P_I)+"   "+str(P_T))
    h.close()
    
# =============================================================================
# Calculate input parameters
# =============================================================================
#    fraction of sediment in 5 sediment class
    NPART = 5 # number of sediment classes (fixed)
    Fcl = np.round(0.26 * clay  , 3) # primary clay
    if clay < 0.25 : # small aggregate
        Fsg = np.round(1.8 * clay, 3)
    elif clay > 0.50 :
        Fsg = np.round(0.6 * clay, 3)
    else :
        Fsg = np.round(0.45-0.6*(clay-0.25), 3)
    Fsi = silt - Fsg # primary silt
    if Fsi < 0 :
        Fsg = silt
        Fsi = 0
    Fsa = np.round(sand*(1-clay)**5, 3) # primary sand
    Flg = np.round(1-Fcl-Fsi-Fsg-Fsa, 3) # large aggregate
    if Flg < 0 :
        Flg = 0
        Fcl = np.round(Fcl / (Fcl+Fsi+Fsg+Fsa), 3)
        Fsg = np.round(Fsg / (Fcl+Fsi+Fsg+Fsa), 3)
        Fsi = np.round(Fsi / (Fcl+Fsi+Fsg+Fsa), 3)
        Fsa = np.round(Fsa / (Fcl+Fsi+Fsg+Fsa), 3)
    
    # particle composition for small and large aggregates
    clay_sg = np.round(clay / (silt+clay), 3)
    silt_sg = np.round(silt / (silt+clay), 3)
    clay_lg = np.round((clay-Fcl-(Fsg*clay_sg))/Flg, 3)
    if clay_lg < clay / 2 :
        Fsg_e = (0.3+0.5*(Fcl+Fsi+Fsg))*(clay+silt) / (1-0.5*(clay+silt))
        clay_lg = np.round((clay-Fcl-(Fsg_e*clay_sg))/Flg, 3)
    silt_lg = np.round(((silt-Fsi)-(Fsg*silt_sg))/Flg ,3)
    sand_lg = np.round((sand-Fsa)/Flg, 3)
    
# =============================================================================
# Code below should be active when adjusted values are not provided directly 
# time varying adjusted factors (0-1) should be provided (KIADJ, KRADJ, SHCRTADJ)   
# =============================================================================
    # baseline interrill erodibility (Kib)
    if sand >= 0.3 :
        vfs_e = vfs
        if vfs > 0.4:
            vfs_e = 0.4 # max vfs set as 0.40 in equation
        Kib = 2728000 + 192100 * vfs_e * 100
        clay_e = clay
    elif sand < 0.3 :
        clay_e = clay
        if clay < 0.10:
            clay_e = 0.10
        Kib = 6054000 - 55130 * clay_e * 100
    
    # baseline rill erodibility (Krb) and critical shear stress (tau_cd)
    if sand >= 0.3 :
        if vfs < 0.4 :
            vsf_e = 0.40
        if clay < 0.4 :
            clay_e = 0.40
        if orgmat < 0.0035 :
            orgmat = 0.0035
        Krb = 0.00197 + 0.0003*vfs_e * 100 + 0.03863*np.exp(-1.84*orgmat*100)
        tau_cd = 2.67 + 0.065*clay_e*100 - 0.058*vfs_e*100
    elif sand < 0.3 :
        if clay < 0.10 :
            clay_e = 0.10
        Krb = 0.0069 + 0.134*np.exp(-20*clay_e*100)
        tau_cd = 3.5
    KIADJ = Kib * KIADJ
    KRADJ = Krb * KRADJ
    SHCRTADJ = tau_cd * SHCRTADJ

    
    
# =============================================================================
# Make input file to run stand-alone WEPP model
# =============================================================================
    f = open("Erosion.INP","w")
    f.write("0    0\n") 
    f.write(str(NPART)+" "+str(width)+"\n") 
    f.write("2.65  2.65  1.80  1.60  2.65\n") # specific gravity of each sediment class
    f.write("\n")
    f.write("1.0 0.10   "+str(length)+"\n")
    f.write(str(Fcl)+"  "+str(Fsi)+"  "+str(Fsg)+"  "+str(Flg)+"  "+str(Fsa)+"\n") # fraction of sediment size class
    f.write("0.002  0.010  0.030  0.300  0.200\n") # diameter of each sediment clss (fixed)
    f.write("0.000  0.000  0.000  "+str(sand_lg)+"  1.000\n") # fraction of sand in each size class
    f.write("0.000  1.000  "+str(silt_sg)+"  "+str(silt_lg)+"  0.000\n") # fraction of silt in each sediment class
    f.write("1.000  0.000  "+str(clay_sg)+"  "+str(clay_lg)+"  0.000\n") # fraction of clay in each sediment class
    f.write("0.111  0.000  0.029  0.012  0.000\n") # fraction of OM in each sediment class
    # fraction of sand, silt, clay, OM in top soil
    f.write(str(sand)+"  "+str(silt)+"  "+str(clay)+"  "+str(orgmat)+"\n") 
    # time varying variables
    f.write(str(KIADJ)+"  "+str(KRADJ)+"  "+str(SHCRTADJ)+"  1.11  "\
            +str(FRCTRL)+"  "+str(RRC)+"\n") 
    f.write("2  "+str(length)+"\n")
    f.write("0.0  "+str(slope)+"  "+str(length)+"  "+str(slope)+"\n")
    f.close()
    
# =============================================================================
# Run stand-alone WEPP model
# =============================================================================
    DETACHED_PROCESS = 0x00000008
    subprocess.call(path+"/Erosion.exe", creationflags=DETACHED_PROCESS)
    
# =============================================================================
# Extract sediment yield for each sediment class from EROSION.OUT file
# (variables: total_yield and sediment)
# =============================================================================
    mylines = []
    with open('EROSION.OUT', 'rt') as myfile :
        for myline in myfile:
            mylines.append(myline)
    
        if mylines[33] == '\n' : # in case of no erosion 
            total_yield = 0
            sediment = np.zeros(NPART)
            
        else :
            A = re.findall("\d+\.\d+", mylines[102])
            total_yield = np.round(float(A[0]) / length , 3) # kg/m2/day    
            sediment = np.zeros(NPART)
            for i in range(NPART):
                B = re.findall("\d+\.\d+", mylines[114+i])
                # detached sediment fraction for each NPART
                fraction = B[-2] 
                # sediment leaving profile for each NPART
                sediment[i] = total_yield * float(fraction) # kg/m/day
    return sediment
#            if round(sum(sediment),0) == round(total_yield,0) :
#                print('Total sediment yield concentration: '+str(total_yield)+' kg/m2/day')
#                print('For each class: '+str(sediment))
#    #            del A, mylines, myline, i, fraction, Q
#            else:
#                print('NG')
            
        
        

    
