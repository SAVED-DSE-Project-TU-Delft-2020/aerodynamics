# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:19:42 2020

@author: Ihab Benyahia
"""


# =============================================================================
# Importing modules 
import numpy as np 
from isacalculator import compute_isa 
import matplotlib.pyplot as plt


# =============================================================================
# Functions 
def compute_reynolds(rho,V,L,mu):
    return (rho*V*L)/mu

def compute_AR(b,S):
    return b**2/S

def compute_sweep(pos,LE_sw,taper,b,c_root):
    return np.tan(LE_sw)-pos*2*c_root*(1-taper)/b

def compute_CLalpha(AR,eff,HC_sw):
    return 2*np.pi*AR/(2+np.sqrt(4+(AR/eff)**2*(1+np.tan(HC_sw)**2)))

def compute_CL(CLa,alpha,alpha0=-2*np.pi/180):
    return CLa*(alpha-alpha0)

def compute_e_straight(AR):
    return 1.78*(1-0.045*AR**0.68) - 0.64

def compute_e_sweep(AR,LE_sw):
    return 4.61*(1-0.045*AR**0.68)*(np.cos(LE_sw))**0.15 - 3.1

def compute_CD(e,AR,CL,CD0):
    return CD0 + CL**2/(np.pi*AR*e)

def compute_C1(taper):
    '''
    This formula was taken from: https://books.google.be/books?id=XtU4HVnWeZIC&pg=PA357&lpg=PA357&dq=C1+datcom+aspect+ratio&source=bl&ots=CgEhRyxJUO&sig=ACfU3U3J1WnRtmN-bZnme2-iXnQNTTdncw&hl=nl&sa=X&ved=2ahUKEwj1vKa3hqzpAhVNDuwKHf5vBg0Q6AEwAXoECAkQAQ#v=onepage&q=C1%20datcom%20aspect%20ratio&f=false
    '''
    return 0.5*np.sin(np.pi*(1-taper)**(1.5+0.8*(np.sin(np.pi*(1-taper)**2))))

def compute_ARcondition(C1,LE_sw):
    return 4/((C1+1)*np.cos(LE_sw))

# =============================================================================
#  Parameters for concept A and D
V_cruise = 28                                                           # [m/s]. Cruise speed
V_transition = 18                                                       # [m/s]. Transition speed 
h_cruise = 500                                                          # [m]. Cruise altitude 
h_transition = 20                                                       # [m]. Transition altitude 
p_cruise,rho_cruise,T_cruise = compute_isa(h_cruise)                    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature.
p_transition,rho_transition,T_transition = compute_isa(h_transition)    # [Pa, kg/m3, K]. Cruise pressure, density, and temperature.
mu_cruise = 17.73e-6                                                    # [Pa.s]. Dynamic viscocity of air at cruise
mu_transition= 17.88e-6                                                 # [Pa.s]. Dynamic viscocity of air at transition 
LE_sweep = 20 / 180 * np.pi                                             # [rad]. Leading edge sweep angle 
MTOM = 16.217                                                           # [kg]. Maximum Take-Off Weight
c_avg = 0.3                                                             # [m]. Average chord length
c_root = 0.76                                                           # [m]. Root chord length                                                         # [m]. Tip chord length
W_S = 121.256                                                           # [N/m2]. Wing loading
g = 9.80665                                                             # [m/s2]. Gravitational acceleration
S = MTOM*g/W_S                                                          # [m2]. Surface area
b = 3                                                                   # [m]. Wing span
c_tip = 2*S/b - c_root                                                  # [m]. Tip chord length
taper = c_tip/c_root                                                    # [-]. Taper ratio
CD0 = 0.003*2                                                           # [-]. Zero-lift drag coefficient
# =============================================================================

# Computed values
AR = compute_AR(b,S)                                                    # [-]. Aspect ratio
HC_sw = compute_sweep(0.5,LE_sweep,taper,b,c_root)                      # [rad]. Half-chord sweep
CLa = compute_CLalpha(AR, 0.95, HC_sw)                                  # [1/rad]. Lift curve slope
e_straight = compute_e_straight(AR)                                     # [-]. Oswald efficiency factor for straight wing 
e_sweep = compute_e_sweep(AR,LE_sweep)                                  # [-]. Oswald efficiency factor for swept wing
C1 = compute_C1(taper)                                                  # [-]. Taper ratio correction factor
condition = compute_ARcondition(C1, LE_sweep)                           # [-]. Condition to check for high or low AR
# =============================================================================

# CL and CD computation
alpha = np.arange(-10,10,0.01)*np.pi/180                                # [rad]. Angle of attack range
CL = compute_CL(CLa,alpha)                                              # [-]. Lift coefficient
CD = compute_CD(e_sweep,AR,CL,CD0)                                      # [-]. Drag coefficient

print("For Concept A and D we have:")
print("Tapered wing: CD0 =",CD0,"K =",1/(np.pi*AR*e_sweep))
print("Elliptical wing: CD0 =",CD0,"K =",1/(np.pi*AR))
print("Rectangular: CD0 =",CD0,"K =",1/(np.pi*AR*e_straight))
print("Maximum L/D =",np.max(CL/CD))

plt.figure(1)
plt.plot(CD,CL)
# =============================================================================
#  Parameters for concept C
taper = 1                                                               # [-]. Taper ratio
LE_sweep = 0                                                            # [rad]. Leading edge sweep angle 
MTOM = 16.217                                                           # [kg]. Maximum Take-Off Weight
W_S = 121.256                                                           # [N/m2]. Wing loading
g = 9.80665                                                             # [m/s2]. Gravitational acceleration
S = MTOM*g/W_S                                                          # [m2]. Surface area
b = 3                                                                   # [m]. Wing span
CD0 = 0.0045*5                                                          # [-]. Zero-lift drag coefficient
# =============================================================================

# Computed values
AR = compute_AR(b,S)                                                    # [-]. Aspect ratio
HC_sw = compute_sweep(0.5,LE_sweep,taper,b,c_root=0)                    # [rad]. Half-chord
CLa = compute_CLalpha(AR, 0.95, HC_sw)                                  # [1/rad]. Lift curve slope
e_straight = compute_e_straight(AR)                                     # [-]. Oswald efficiency factor for straight wing 
C1 = compute_C1(taper)                                                  # [-]. Taper ratio correction factor
condition = compute_ARcondition(C1, LE_sweep)                           # [-]. Condition to check for high or low AR
# =============================================================================

# CL and CD computation
alpha = np.arange(-10,10,0.01)*np.pi/180                                # [rad]. Angle of attack range
CL = compute_CL(CLa,alpha)                                              # [-]. Lift coefficient
CD = compute_CD(e_sweep,AR,CL,CD0)                                      # [-]. Drag coefficient

print("\nFor Concept C we have:")
print("Elliptical wing: CD0 =",CD0,"K =",1/(np.pi*AR))
print("Rectangular: CD0 =",CD0,"K =",1/(np.pi*AR*e_straight))
print("Maximum L/D =",np.max(CL/CD))
plt.figure(1)
plt.plot(CD,CL) 
