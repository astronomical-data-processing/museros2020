# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 17:41:43 2020

@author: LJCHEN
"""

import numpy as np
import math as mp
import os
from muser.data_models.parameters import muser_path, muser_data_path
from astropy.time import Time
def sind(x):
    return np.sin(x*np.pi/180)
def cosd(x):
    return np.cos(x*np.pi/180)
def tand(x):
    return np.tan(x*np.pi/180)

def get_sun(obs_time):
    Y= obs_time.datetime.year
    M = obs_time.datetime.month
    D = obs_time.datetime.day
    h = obs_time.datetime.hour
    m = obs_time.datetime.minute
    s = obs_time.datetime.second
    ms = obs_time.datetime.microsecond
    # Y,M,D,h,m,s,ms = eval(input("UT Time YMDhmsmS:"))
    ephem_data_path = muser_path('configurations')

    Day = D + (ms/3600000000 + s/3600 + m/60 + h)/24
    DeltaT = np.loadtxt(os.path.join(ephem_data_path,'TData.txt'))
    T_Num = (Y-1973)*12+M
    Del_T = DeltaT[T_Num-1,3]
    if M==1 or M==2:
        Y = Y-1
        M = M+12
    B = 2-int(Y/100)+int(int(Y/100)/4)
    h0s  = -0.8333   #Sunrise elevation angle
    Lgt = 115.2505  # local longitude in degree
    Lat = 42.211833333  # local latitude in degree
    E = 1365  # local altitude in meter
    P = 820 # atmosphere pressure in mbar
    T = 11  # temperature in degree C
    JD = int(365.25*(Y+4716))+int(30.6001*(M+1))+Day+B-1524.5        # 3.1.1
    JDE = JD+Del_T/86400                                             # 3.1.2
    JC = (JD-2451545)/36525                                          # 3.1.3
    JCE = (JDE-2451545)/36525                                        # 3.1.3
    JME = JCE/10                                                     # 3.1.4
    L0T = np.loadtxt(os.path.join(ephem_data_path,'L0.txt'))                                  # 3.2.2
    L1T = np.loadtxt(os.path.join(ephem_data_path,'L1.txt'))
    L2T = np.loadtxt(os.path.join(ephem_data_path,'L2.txt'))
    L3T = np.loadtxt(os.path.join(ephem_data_path,'L3.txt'))
    L4T = np.loadtxt(os.path.join(ephem_data_path,'L4.txt'))
    L5T = np.loadtxt(os.path.join(ephem_data_path,'L5.txt'))
    L0 = np.sum(L0T[:,1]*np.cos(L0T[:,2]+L0T[:,3]*JME),axis=0)       # 3.2.2
    L1 = np.sum(L1T[:,1]*np.cos(L1T[:,2]+L1T[:,3]*JME),axis=0)       # 3.2.3
    L2 = np.sum(L2T[:,1]*np.cos(L2T[:,2]+L2T[:,3]*JME),axis=0)       # 3.2.3
    L3 = np.sum(L3T[:,1]*np.cos(L3T[:,2]+L3T[:,3]*JME),axis=0)       # 3.2.3
    L4 = np.sum(L4T[:,1]*np.cos(L4T[:,2]+L4T[:,3]*JME),axis=0)       # 3.2.3
    L5 = np.sum(L5T[1]*np.cos(L5T[2]+L5T[3]*JME),axis=0)             # 3.2.3
    L= (L0+L1*JME+L2*JME**2+L3*JME**3+L4*JME**4+L5*JME**5)/10**8     # 3.2.4
    Ld = (L*180/np.pi)%360                                           # 3.2.5
    B0T = np.loadtxt(os.path.join(ephem_data_path,'B0.txt'))
    B1T = np.loadtxt(os.path.join(ephem_data_path,'B1.txt'))                                      # 3.2.7
    B0 = np.sum(B0T[:,1]*np.cos(B0T[:,2]+B0T[:,3]*JME),axis=0)       # 3.2.7
    B1 = np.sum(B1T[:,1]*np.cos(B1T[:,2]+B1T[:,3]*JME),axis=0)       # 3.2.7
    Bs= (B0+B1*JME)/10**8                                            # 3.2.7
    Bd = Bs*180/np.pi                                                # 3.2.7
    R0T = np.loadtxt(os.path.join(ephem_data_path,'R0.txt'))                                       # 3.2.8
    R1T = np.loadtxt(os.path.join(ephem_data_path,'R1.txt'))
    R2T = np.loadtxt(os.path.join(ephem_data_path,'R2.txt'))
    R3T = np.loadtxt(os.path.join(ephem_data_path,'R3.txt'))
    R4T = np.loadtxt(os.path.join(ephem_data_path,'R4.txt'))
    R0 = np.sum(R0T[:,1]*np.cos(R0T[:,2]+R0T[:,3]*JME),axis=0)       # 3.2.8
    R1 = np.sum(R1T[:,1]*np.cos(R1T[:,2]+R1T[:,3]*JME),axis=0)       # 3.2.8
    R2 = np.sum(R2T[:,1]*np.cos(R2T[:,2]+R2T[:,3]*JME),axis=0)       # 3.2.8
    R3 = np.sum(R3T[:,1]*np.cos(R3T[:,2]+R3T[:,3]*JME),axis=0)       # 3.2.8
    R4 = np.sum(R4T[1]*np.cos(R4T[2]+R4T[3]*JME),axis=0)             # 3.2.8
    R= (R0+R1*JME+R2*JME**2+R3*JME**3+R4*JME**4)/10**8               # 3.2.8
    Theta = (Ld+180)%360
    Beta = -Bd                                              # 3.3.1
    X0 = 297.85036+445267.111480*JCE-0.0019142*JCE**2+JCE**3/189474  # 3.4.1
    X1 = 357.52772+35999.050340*JCE-0.0001603*JCE**2-JCE**3/300000   # 3.4.2
    X2 = 134.96298+477198.867398*JCE+0.0086972*JCE**2+JCE**3/56250   # 3.4.3
    X3 = 93.27191+483202.017538*JCE-0.0036825*JCE**2+JCE**3/327270   # 3.4.4
    X4 = 125.04452-1934.136261*JCE+0.0020708*JCE**2+JCE**3/450000    # 3.4.5
    Table43 = np.loadtxt(os.path.join(ephem_data_path,'TableA43.txt'))
    Del_Phi = sum((Table43[:,5]+Table43[:,6]*JCE)*sind(X0*Table43[:,0]+X1*Table43[:,1]+X2*Table43[:,2]+X3*Table43[:,3]+X4*Table43[:,4]))/36000000   # 3.4.7
    Del_Eps = sum((Table43[:,7]+Table43[:,8]*JCE)*cosd(X0*Table43[:,0]+X1*Table43[:,1]+X2*Table43[:,2]+X3*Table43[:,3]+X4*Table43[:,4]))/36000000   # 3.4.8
    U = JME/10
    Eps0=84381.448-4680.93*U-1.55*U**2+1999.25*U**3-51.38*U**4-249.67*U**5-39.05*U**6+7.12*U**7+27.87*U**8+5.79*U**9+2.45*U**10   # 3.5.1
    Eps = Eps0/3600+Del_Eps                                          # 3.5.2
    Del_Tao = -20.4898/(3600*R)                                      # 3.6
    Lamda = Theta+Del_Phi+Del_Tao                                    # 3.7
    V0 = 280.46061837+360.98564736629*(JD-2451545)+0.000387933*JC**2-JC**3/38710000   # 3.8.1
    V0 = V0%360
    V = V0+Del_Phi*cosd(Eps)
    Alpha = mp.atan2((sind(Lamda)*cosd(Eps)-tand(Beta)*sind(Eps)),cosd(Lamda))    # 3.9.1
    Alpha = Alpha *180/np.pi
    if Alpha<0:
        Alpha = Alpha+360
    Delta = mp.asin(sind(Beta)*cosd(Eps)+cosd(Beta)*sind(Eps)*sind(Lamda))    # 3.10
    Delta = Delta*180/np.pi
    H = V+Lgt-Alpha                                                  # 3.11
    H= H%360                                                         # 3.11
    Epslo = 8.794/(3600*R)                                           # 3.12.1
    u = mp.atan(0.99664719*tand(Lat))                                # 3.12.2
    x = np.cos(u)+E*cosd(Lat)/6378140                                # 3.12.3
    y = 0.99664719*np.sin(u)+E*sind(Lat)/6378140                     # 3.12.4
    Del_Alpha = mp.atan2((-x*sind(Epslo)*sind(H)),(cosd(Delta)-x*sind(Epslo)*cosd(H)))  # 3.12.5
    Del_Alpha = Del_Alpha*180/np.pi
    AlphaS = Alpha+Del_Alpha                                         # 3.12.6
    DeltaS = mp.atan2(((sind(Delta)-y*sind(Epslo))*cosd(Del_Alpha)),(cosd(Delta)-x*sind(Epslo)*cosd(H)))  # 3.12.7
    DeltaS = DeltaS*180/np.pi
    HS = H - Del_Alpha                                               # 3.13
    e0 = mp.asin(sind(Lat)*sind(DeltaS)+cosd(Lat)*cosd(DeltaS)*cosd(HS)) # 3.14.1
    e0 = e0*180/np.pi
    Delta_e = (P/1010)*(283/(273+T)*(1.02/(60*tand(e0+10.3/(e0+5.11))))) # 3.14.2
    e = e0+Delta_e                                                   # 3.14.3
    Thete_z = 90-e                                                   # 3.14.4
    Gamma = mp.atan2(sind(HS),(cosd(HS)*sind(Lat)-tand(DeltaS)*cosd(Lat))) # 3.15.1
    Gamma = Gamma*180/np.pi
    Phi = (Gamma+180)%360 # eastward from north                      # 3.15.2
    return Alpha,Delta, H,Thete_z, Phi
    # print("Right Ascension:",Alpha)
    # print("Declination:",Delta)
    # print("Hour Angle:",H)
    # print("Zenith Angle:",Thete_z)
    # print("Azimuth Angle:",Phi)