# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
# EOS
def CH4_EOS(Temp,P_PSI):
    T1=Temp+273.15 #temp (K)
    Tc=190.6 #critical temp (K)
    Pc=4.6e06 #critical pressure (Pa)
    P=P_PSI*6894.757 #pressure final (Pa)
    w=.008
    r=8.314 #gas constant (kg m2 s-2 mol-1 k-1)

    Tr=T1/Tc

    if w<=.49:
        k=.37464+1.54226*w-.26992*w**2
    else:
        k=.379664+w*(1.48503+w*(-.164423+.016666*w))
        
# Solve cubic EOS (z)

    A=.457236*((Tc/T1)**2)*(P/Pc)*((1+k*(1-np.sqrt(Tr)))**2)
    B=.07780*(Tc/T1)*(P/Pc)

    x=(-1+B)
    y=(A-3*B**2-2*B)
    n=(-A*B+B**2+B**3)

    Q=(3*y-x**2)/9
    R=(9*x*y-27*n-2*x**3)/54

    D=Q**3+R**2
    S=(R+np.sqrt(D))**(1/3)
    T=(R-np.sqrt(D))**(1/3)

    z1=(-1/3)*x+(S+T)
    z2=(-1/3)*-(1/2)*(S+T)+(1/2)*1j*3**.5*(S-T)
    z3=(-1/3)*x-(1/2)*(S+T)-(1/2)*1j*3**.5*(S-T)
    z=z1

#Calc molar density
    f=P/(z*r*T1) #density (mol/m3)
    return [f,z]