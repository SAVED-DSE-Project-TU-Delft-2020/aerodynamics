# -*- coding: utf-8 -*-
"""
Created on Fri May 29 11:14:57 2020

@author: Ihab Benyahia
"""


import numpy as np
import matplotlib.pyplot as plt

# ---------------------------- General functions ---------------------------- #
def compute_sweep(LE_sw,taper,x_c,cr,b):
    return (np.tan(LE_sw*np.pi/180) - x_c*2*cr/b*(1-taper))*180/np.pi

def compute_cldes(m,rho,V,QC_sw,S):
    return (m*9.81)/(S*0.5*rho*(V*np.cos(QC_sw*np.pi/180))**2)

#===========AR condition======================================================
def compute_C1(taper):
    '''
    This formula was taken from: https://books.google.be/books?id=XtU4HVnWeZIC&pg=PA357&lpg=PA357&dq=C1+datcom+aspect+ratio&source=bl&ots=CgEhRyxJUO&sig=ACfU3U3J1WnRtmN-bZnme2-iXnQNTTdncw&hl=nl&sa=X&ved=2ahUKEwj1vKa3hqzpAhVNDuwKHf5vBg0Q6AEwAXoECAkQAQ#v=onepage&q=C1%20datcom%20aspect%20ratio&f=false
    '''
    return 0.5*np.sin(np.pi*(1-taper)**(1.5+0.8*(np.sin(np.pi*(1-taper)**2))))

def compute_ARcondition(C1,LE_sw,AR):
    '''
    Condition from p.568
    '''
    condition = 4/((C1+1)*np.cos(LE_sw*np.pi/180))
    if AR > condition:
        return "High AR"
    else:
        return "Low AR"

# ---------------------------------- AoA ----------------------------------- #
        
def compute_AOA_0_lift(clalpha,cldes,aoa_des):
    '''
    Formula from p.434
    No twist present
    '''
    return aoa_des - cldes/clalpha

def compute_AOA_0_lift_tw(aoa0_theta,theta,aoa_zero):
    '''
    Formula from p.434
    aoa0_theta from Figure 4.1.3.1-4 p.437-438
    '''
    return aoa_zero + aoa0_theta * theta
# ---------------------------------- Lift ---------------------------------- #
#===========High AR CLalpha===================================================
    
def compute_CLa(AR,HC_sw):
    
    return 2*np.pi*AR/(2+np.sqrt((AR/0.95)**2*(1+np.tan(HC_sw*np.pi/180)**2) + 4))

#===========Low AR CLalpha====================================================
    
# def compute_CLa_theory(AR,taper,LE_sw):
#     return 8*np.arctan2(np.pi*AR/(16+np.pi*AR/(1+2*taper*np.tan(LE_sw))))

# def compute_CLa_basic(CLa_th,CLa_ratio):
#     '''
#     CLa_ratio from Figure 4.1.3.2-50a p.???
#     '''
#     return CLa_th*CLa_ratio

# def compute_CLIII(x_c,LE_sw,taper):
#     '''
#     CLIII can be read of using the outcome of this function
#     '''
#     return x_c/((1+taper)*(1-x_c)*np.tan(LE_sw))

# def compute_CLa_limit(LE_sw):
#     return 0.0467*np.sin(np.pi/2 - LE_sw)**0.25 * 180/np.pi

# def compute_deltaCLa_II(CLa_limit,CLa_basic):
    
#     if CLa_limit - CLa_basic >= 0.0067*180/np.pi:
#         return 0.0067*180/np.pi
#     elif 0 < CLa_limit - CLa_basic < 0.0067*180/np.pi:
#         return CLa_limit - CLa_basic
#     else:
#         return 0

# def compute_CLa_II(CLa_basic,deltaCLa_II):
#     return CLa_basic + deltaCLa_II

# def compute_deltaCLa_III(CLa_limit,CLa_basic):
    
#     if CLa_limit - CLa_basic >= 0.0067*180/np.pi:
#         return 0.012*180/np.pi
#     elif 0 < CLa_limit - CLa_basic < 0.012*180/np.pi:
#         return CLa_limit - CLa_basic
#     else:
#         return 0
    
# def compute_CLa_III(CLa_basic,deltaCLa_III):
#     return CLa_basic + deltaCLa_III

#===========Stall characteristics for low AR==================================

def compute_CLmax_low(CLmax_base,delta_CLmax):
    '''
    Formula from p.570
    CLmax_base from Figure 4.1.3.4-23a/b p.587
    delta_CLmax from Figure 4.1.3.4-24a p.588
    '''
    return CLmax_base + delta_CLmax

def compute_aoa_stall_low(aoa_stall_base,delta_aoa):
    '''
    Formula from p.570
    aoa_stall_base from Figure 4.1.3.4-25a p.589
    delta_aoa from Figure 4.1.3.4-25b p.589
    '''
    return aoa_stall_base + delta_aoa

#===========Stall characteristics for high AR=================================

def compute_CLmax_high(CL_cl,clmax,delta_CLmax):
    '''
    Formula from p.568
    CL_cl from Figure 4.1.3.4-21a p.585
    delta_CLmax from Figure 4.1.3.4-22 p.586
    '''
    return CL_cl*clmax + delta_CLmax

def compute_aoa_stall_high(CLmax,CLalpha,aoa_0,delta_aoa):
    '''
    Formula from p.568
    delta_aoa from Figure 4.1.3.4-21b p.585
    '''
    return 180*CLmax/CLalpha/np.pi + aoa_0 + delta_aoa

#===========CL below stall====================================================

def compute_CN_prime_CLmax(CLmax,aoa_stall):
    '''
    Formula from p.507
    '''
    return CLmax/np.cos(aoa_stall*np.pi/180)
    
def compute_Jpar(C1,C2,LE_sw,AR):
    '''
    Formula from p.507
    Obtain delta_CNaa from this
    C2 from Figure 4.1.3.4-24b p.588
    '''
    return 0.3*(C1+1)*AR*np.cos(LE_sw*np.pi/180)*((C1+1)*(C2+1)-((C2+1)*AR*np.tan(LE_sw*np.pi/180)/7)**3)

def compute_linear(xb,xe,yb,ye,x):
    a = (ye - yb)/(xe - xb)
    y0 = yb - a*xb
    return a * x + y0
    
def compute_CNaa_ref(CN_prime,CLalpha,aoa_stall):
    '''
    Formula from p.509
    CN_prime @ CLmax
    '''
    return (CN_prime - CLalpha*0.5*np.sin(2*np.pi*aoa_stall/180))/(np.sin(np.pi*aoa_stall/180)*np.abs(np.sin(np.pi*aoa_stall/180)))

# def compute_aoa_camber(aoa,aoa_0,aoa_stall):
    
#     return aoa*(1+aoa_0/(90-aoa_stall)) - 90*aoa_0/(90-aoa_stall)

def compute_CNaa_below(CNaa_ref,delta_CNaa):
    '''
    Formula from p.507
    delta_CNaa from Figure 4.1.3.3-55a p.558
    '''
    return CNaa_ref + delta_CNaa

# def compute_CNaa_after(CNaa_ref,CNaa_90,aoa_stall,aoa,D,CLalpha,CLmax,CLmax_low):
#     '''
#     CNaa_90 from Figure 4.1.3.3-55b p.558
#     D from Figure 4.1.3.3-55a p.558
#     CLmax_low is the CLmax for low AR
#     aoa is for cambered wings
#     '''
#     return CNaa_ref + (CNaa_90 - CNaa_ref)*(1-np.tan(aoa_stall*np.pi/180)/np.tan(aoa*np.pi/180)) + D*CLalpha*CLmax**2/(2.3*CLmax_low**2)


def compute_CN_prime(CLalpha,aoa,CNaa):
    '''
    Formula from p.506
    '''
    return CLalpha*np.sin(2*np.pi*aoa/180)/2 + CNaa*np.sin(np.pi*aoa/180)**2

def compute_CL(CN_prime,aoa):
    '''
    Formula from p.506
    '''
    return CN_prime*np.cos(aoa*np.pi/180)

# ---------------------------------- Drag ---------------------------------- #
    
#===========CD0 of the wing===================================================
def compute_CD0_wing(Cf,tc_avg,x_tcmax,Swet,Sref,Rls):
    '''
    Formula from p.663
    Cf from Figure 4.1.5.1-26/27 p.687-688
    Rls from Figure 4.1.5.1-28b p.689
    '''
    if x_tcmax >= 0.3:
        return Cf*(1+1.2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref
    else:
        return Cf*(1+2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref

#===========CD0 of the body===================================================

def compute_CD0_body(Cf,lb,Ss_Sb,Sb):
    '''
    Formula from p.879
    Cf from Figure 4.1.5.1-26/27 p.687-688
    Ss_Sb from Figure 2.3-2 and 2.3-3 p.179-180
    '''
    d = np.sqrt(Sb/0.7854)
    
    return Cf*(1+60/(lb/d)**3 + 0.0025*lb/d)*Ss_Sb

#===========CD0 of the tail===================================================

def compute_CD0_tail(Cf,tc_avg,x_tcmax,Swet,Sref,Rls):
    '''
    Same method as for the wing
    Cf from Figure 4.1.5.1-26/27 p.687-688
    Rls from Figure 4.1.5.1-28b p.689
    
    '''
    if x_tcmax >= 0.3:
        return Cf*(1+1.2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref
    else:
        return Cf*(1+2*(tc_avg)+100*(tc_avg)**4)*Rls*Swet/Sref


#===========Lift induced drag of the wing=====================================

def compute_CD_ind_wing(CLalpha,AR,clalpha,v,w,theta,R,CL):
    '''
    Formula from p.699
    v from Figures 4.1.5.2-42a-i p.737-742
    w from Figures 4.1.5.2-48a-i p.743-747
    R from Figure 4.1.5.2-53 p.748

    '''
    e = 1.1*CLalpha/(AR*(R*CLalpha/AR + (1-R)*np.pi))
    return CL**2/(np.pi*AR*e) + CL*theta*clalpha*v + (theta*clalpha)**2*w,e

# #===========Lift induced drag of the body=====================================

# def compute_CD_ind_body(aoa,S0,Vb,x0_lb,k2k1,lb,eta,d,cdc = 1.2):
#     '''
#     k2k1 from Figure 4.2.1.1-20a p.774
#     x0_lb from Figure 4.2.1.1-20b p.774
#     eta from Figure 4.2.1.2-35a p.815
#     cdc from Figure 4.2.1.2-35a p.815 and is 1.2 since we fly at low Mach number
#     '''
#     x0 = x0_lb*lb
#     return 2*k2k1*S0/Vb**(2/3)*aoa**2 + eta*cdc*d/4*(lb-x0)*aoa**3

# def compute_CD_ind_1(AR,LE_sw,t_MAC,lN_d,lA_d,taper,LER_MAC,theta,ymax_MAC,CLdes,Re,B,TR = 0):
#     '''
#     There are some constraints on when this method can be used
#     B coefficients on p.1130
#     '''
#     return (B[0] + B[1]*1/AR + B[2]*AR + B[3]*np.sqrt(np.tan(LE_sw)) + B[4]*t_MAC +
#             B[5]*lN_d + B[6]*lA_d + B[7]*taper + B[8]*taper**2 + B[9]*taper**3 +
#             B[10]*TR + B[11]*LER_MAC + B[12]*theta + B[13]*ymax_MAC + B[14]*CLdes +
#             B[15]*Re)

# def compute_CD_ind_2(CD_ind_w,CD_ind_b,Sb,S):
#     return CD_ind_w + CD_ind_b*Sb/S


def compute_CD(CD0,CDi):
    return CD0 + CDi
