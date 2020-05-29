# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:14:57 2020

@author: Ihab Benyahia
"""


import numpy as np
import matplotlib.pyplot as plt
#%%

def compute_C1(taper):
    '''
    This formula was taken from: https://books.google.be/books?id=XtU4HVnWeZIC&pg=PA357&lpg=PA357&dq=C1+datcom+aspect+ratio&source=bl&ots=CgEhRyxJUO&sig=ACfU3U3J1WnRtmN-bZnme2-iXnQNTTdncw&hl=nl&sa=X&ved=2ahUKEwj1vKa3hqzpAhVNDuwKHf5vBg0Q6AEwAXoECAkQAQ#v=onepage&q=C1%20datcom%20aspect%20ratio&f=false
    '''
    return 0.5*np.sin(np.pi*(1-taper)**(1.5+0.8*(np.sin(np.pi*(1-taper)**2))))

def compute_ARcondition(C1,LE_sw):
    return 4/((C1+1)*np.cos(LE_sw))

def compute_CLa_theory(8*np.atan())
    
#=============================================================================
def compute_CD0_wing(Cf,tc_avg,x_tcmax,Swet,Sref,Rls):
    '''
    Cf from Figure 4.1.5.1-26/27 p.687-688
    Rls from Figure 4.1.5.1-28b p.689
    '''
    if x_tcmax >= 0.3:
        return Cf*(1+1.2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref
    else:
        return Cf*(1+2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref


def compute_CD0_body(Cf,lb,Sb,Ss,db):
    '''
    Cf from Figure 4.1.5.1-26/27 p.687-688
    Ss/Sb from Figure 2.3-2 and 2.3-3 p.179-180
    '''
    d = np.sqrt(Sb/0.7854)
    
    return Cf*(1+60/(lb/d)**3 + 0.0025*lb/d)*Ss/Sb

# Total CD0 is (CD0wing + CD0body)*Rwb, whereby Rwb is a correlation factor 
# from Figure 4.3.3.1-37 
def compute_CD_ind_wing(CL,AR,cla,v,w,theta,R,CLa):
    '''
    v from Figures 4.1.5.2-42a-i p.737-742
    w from Figures 4.1.5.2-48a-i p.743-747
    R from Figure 4.1.5.2-53 p.748

    '''
    e = 1.1*CLa/(AR*(R*CLa/AR + (1-R)*np.pi))
    return CL**2/(np.pi*AR*e) + CL*theta*cla*v + (theta*cla)**2*w

def compute_CD_ind_body(aoa,S0,Vb,k2k1,lb,eta,cdc = 1.2,d):
    '''
    k2k1 from Figure 4.2.1.1-20a p.774
    x0_lb from Figure 4.2.1.1-20b p.774
    eta from Figure 4.2.1.2-35a p.815
    cdc from Figure 4.2.1.2-35a p.815 and is 1.2 since we fly at low Mach number
    '''
    x0 = x0_lb*lb
    return 2*k2k1*S0/Vb**(2/3)*aoa**2 + eta*cdc*d/4*(lb-x0)*aoa**3
def compute_CD_ind_1(AR,LE_sw,t_MAC,lN_d,lA_d,taper,TR = 0,LER_MAC,theta,ymax_MAC,CLdes,Re,B):
    '''
    There are some constraints on when this method can be used
    B coefficients on p.1130
    '''
    return (B[0] + B[1]*1/AR + B[2]*AR + B[3]*np.sqrt(np.tan(LE_sw)) + B[4]*t_MAC +
            B[5]*lN_d + B[6]*lA_d + B[7]*taper + B[8]*taper**2 + B[9]*taper**3 +
            B[10]*TR + B[11]*LER_MAC + B[12]*theta + B[13]*ymax_MAC + B[14]*CLdes +
            B[15]*Re)

def compute_CD_ind_2(CD_ind_w,CD_ind_b,Sb,S):
    return CD_ind_w + CD_ind_b*Sb/S


def compute_CD(CD0,CDi):
    return CD0 + CDi
