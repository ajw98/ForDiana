#!/usr/bin/env python

# This code is used to find the Chi interaction parameters for making phase diagrams


import numpy as np

R = 8.3145 # mL MPa/ mol K
T = 298 #K
#R^2_a = 4(ddp-dds)^2 + (dpp-dps)^2 + (dhp-dhs)^2

def R_sqr(dp,pp,hp,ds,ps,hs):
    r = 4 * (dp - ds)**2 + (pp - ps)**2 + (hp - hs)**2
    return r


#Cellulose Acetate
Pd = 14.9
Pp = 7.1
Ph = 11.1
#Glacial Acetic Acid
Sd = 14.5
Sp = 8
Sh = 13.5
#Water 
wd = 18.1
wp = 17.1
wh = 16.9
Rps = R_sqr(Pd,Pp,Ph,Sd,Sp,Sh) # Polymer and Solvent
Rpn = R_sqr(Pd,Pp,Ph,wd,wp,wh) # Polymer and Nonsolvent
#Rpn_ = R_sqr(Pd,Pp,Ph,wd,wp,wh)#Polymer and corrected water
#Rsn = R_sqr(Sd,Sp,Sh,Wd,Wp,Wh) # Solvent and Water 
Rsn = R_sqr(Sd,Sp,Sh,wd,wp,wh) # Solvent and corrected Water 
rho_s = 1.049 #g/mL
rho_p = 1.3 #g/mL
rho_w = 0.998 #g/mL

PMW = 162.15 #g/mol
SMW = 60.05 #g/mol
WMW = 18.015 #g/mol

PV = PMW/rho_p #mL/mol
SV = SMW/rho_s #mL/mol
WV = WMW/rho_w #mL/mol
V = (SV)
Chi13 = SV * Rpn / (4 * R * T) #Chi between P and w
Chi12 = SV * Rps / (4 * R * T) #Chi between P and PC 
Chi23 = SV * Rsn / (4 * R * T) # Chi between PC and W 
#Chi23_2b = V * Rsn_2b / (4 * R * T) # Chi between Pc and corrected w 
#print(PV,TPV,PCV,WV,V,'Vs')
#print('Polymer and 2-Pyr',Chi12)
print('Polymer and Solvent',Chi12)
print('Polymer and nonsolvent',Chi13)
#print('Polymer and corrected water',Chi13_)
#print('2-Pyr and water',Chi23)
#print('2-Pyr and corrected water',Chi23b)
print('Solvent and nonsolvent',Chi23)
#print('PC and corrected water',Chi23_2b)
