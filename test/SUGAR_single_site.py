#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Substrate Utilisation by Growth and Autotrophic Respiration (SUGAR) model
for a single site.

Instructions for use:

1. #-- Driving data:
If running for Caxiuana using JULES output, download jules.all.nc. 
In this script change the variable input_data_path to the path name of 
the directory in which you have saved jules.all.nc. If running with other data, 
you require GPP (kg m-2 s-1), canopy temperature (dC) and biomass (kg m-2) for 
running the model and NPP (kg m-2 s-1) to initialise the model if a value of alpha 
is not specified. If data is not in .nc file or variable names within nc are different,
#-- Driving data should be changed accordingly so that your input data is 
converted to numpy arrays.

2. #--Initialisation data
If values of phi and alpha are not specified by the user (see 3. #--Parameters)
they are estimated using a period of initialisation data. This should contain GPP, 
NPP and biomass. In the Caxiuana runs the first year of data from the driving 
dataset is used for this.

3. #--Parameters
All parameters can be changed. The default parameters below were used in the original study.
@author: simon jones
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from netCDF4 import Dataset


#-- Functions
def FQ(T):
	return q10**(0.1*(T-T_ref))

def f_R_maintenance(T,C_NSC,C_v):
	return R_m0*FQ(T)*C_v*C_NSC/(C_NSC+K_m*C_v)

def f_Growth(T,C_NSC,C_v):
	return G_0*FQ(T)*C_v*C_NSC/(C_NSC+K_m*C_v)

#-- Parameters
Y_g     = 0.75     # Growth efficiency coefficient (Thornley and Johnson, 1990)
a_km    = 0.5      # Relates half saturation constant (Km) to equilibrium NSC mass fraction (fNSC)
q10     = 2.0      # parameter in FQ(T) = exp(0.1*ln(q10)(T-Tref))
T_ref   = 25.0     # Reference temp in Q10 function ^ (degrees C)
f_NSC   = 0.08     # Equilibrium NSC mass fraction (default = 0.08 (Wurth et. al, 2005)
phi     = None     # Maximum specific rate of carbohydrate utilisation (s-1). If == None, it is calculated using initialisation data (default)
alpha   = None     # Ratio of growth to PCE. If == None, calculated using initialisation data (default)
     

#-- Switches
l_cv_const      = True # If True assume biomass is constant (default since SUGAR not coupled to DGVM but can be set to false and be prescribed Cv)


#-- Driving data
input_data_path = ''                                                     # file path for the driving GPP and temperature data
dataset_drive   = Dataset(input_data_path+'jules.all.nc')                # 'jules.all.nc' is the dataset that contains all output from the JULES simulations.
start_date      = pd.datetime(2001,1,1)                                  # start date of input data
end_date        = pd.datetime(2016,12,9,13,0)                            # end date of input data
dt              = 3600.0                                                 # timestep of input data in seconds.  
data_period     = pd.date_range(start_date,end_date,freq = str(dt)+'S')
P_len           = len(data_period)
gpp_gb          = dataset_drive.variables['gpp_gb'][:,0,0]               # Driving Gridbox GPP (kg m-2 s-1)
T1p5m_gb        = dataset_drive.variables['t1p5m_gb'][:,0,0]             # Driving Gridbox temperature (degrees C)
if l_cv_const == False:
	C_v_gb  = dataset_drive.variables['cv'][:,0,0]                   # Gridbox structural biomas (kg m-2)

#-- Initialisation
# only used if phi and alpha are not prescribed
C_v_init    = dataset_drive.variables['cv'][0:8760,0,0]
gpp_gb_init = dataset_drive.variables['gpp_gb'][0:8760,0,0]
npp_gb_init = dataset_drive.variables['npp_gb'][0:8760,0,0]
if l_cv_const == True:
	C_v_gb = np.ones(P_len)*np.mean(C_v_init)

if phi == None:
	tau = np.mean(C_v_init)/np.mean(gpp_gb_init)
	phi = (1.0+a_km)/(tau*np.mean(FQ(T1p5m_gb)))

if alpha == None:
    alpha = np.mean(npp_gb_init)/np.mean(gpp_gb_init)

K_m   = a_km*f_NSC
G_0   = alpha*phi
R_m0  = (1.0-alpha/Y_g)*phi


#-- Create output arrays
Rm_gb    = np.zeros(P_len)  # Gridbox maintenance respiration (kg m-2 s-1)
Rg_gb    = np.zeros(P_len)  # Gridbox growth respiration (kg m-2 s-1)
G_gb     = np.zeros(P_len)  # Gridbox growth (kg m-2 s-1)
U_gb     = np.zeros(P_len)  # Gridbox total carbohydrate utilisation (kg m-2 s-1)
C_NSC_gb = np.zeros(P_len)  # Gridbox NSC content (kg m-2)

#-- Initial conditions
C_NSC_gb[0] = f_NSC*C_v_gb[0]
Rm_gb[0]    = f_R_maintenance(T1p5m_gb[0],C_NSC_gb[0],C_v_gb[0])
G_gb[0]     = f_Growth(T1p5m_gb[0],C_NSC_gb[0],C_v_gb[0])
Rg_gb[0]    = ((1.0-Y_g)/Y_g)*G_gb[0]
U_gb[0]     = Rm_gb[0]+Rg_gb[0]+G_gb[0]

#-- Start cycle
for i in range(1,P_len):
	Rm_gb[i]    = f_R_maintenance(T1p5m_gb[i],C_NSC_gb[i-1],C_v_gb[i-1])
	G_gb[i]     = f_Growth(T1p5m_gb[i],C_NSC_gb[i-1],C_v_gb[i-1])
	Rg_gb[i]    = ((1.0-Y_g)/Y_g)*G_gb[i]
	U_gb[i]     = Rm_gb[i]+Rg_gb[i]+G_gb[i]
	C_NSC_gb[i] = max(0,C_NSC_gb[i-1]+(gpp_gb[i]-U_gb[i])*dt)

#-- Plotting with Pandas
Rm_gb_plot    = pd.Series(Rm_gb,index = data_period) 
Rg_gb_plot    = pd.Series(Rg_gb,index = data_period)
G_gb_plot     = pd.Series(G_gb,index = data_period)
U_gb_plot     = pd.Series(U_gb,index = data_period)
C_NSC_gb_plot = pd.Series(C_NSC_gb,index = data_period)
gpp_gb_plot  = pd.Series(gpp_gb,index = data_period)

plt.figure()
gpp_gb_plot.resample('M').mean().plot(color ='blue',label = 'Monthly mean GPP')
gpp_gb_plot.resample('D').mean().plot(color = 'blue',alpha = 0.15,label = 'Daily mean GPP')
U_gb_plot.resample('M').mean().plot(color ='orange',label = 'Monthly mean PCE')
U_gb_plot.resample('D').mean().plot(color ='orange',alpha = 0.25,label = 'Monthly mean PCE')
plt.title('GPP vs PCE',size = 14)
plt.ylabel('Flux (kg m$^{-2}$ s$^{-1}$)',size = 14)
plt.xlabel('Date',size = 14)
plt.legend()
plt.show()

plt.figure()
C_NSC_gb_plot.resample('D').mean().plot()
plt.title('NSC content',size = 14)
plt.ylabel('NSC content (kg m$^{-2}$)',size = 14)
plt.xlabel('Date',size = 14)
plt.show()

