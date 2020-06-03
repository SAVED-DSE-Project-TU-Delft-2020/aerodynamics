""" 
This program calculates air pressure, air temperature and air density based on input altitude 

Author: Casper Kanaar 

"""
# =============================================================================
# Main

def compute_isa(h):
    g = 9.80665
    R = 287
    from math import e

    #Lists representing the base values for the lapse rate and the altitudes for each of the layers.

    h0_list = [0,11000,20000,32000,47000,51000,71000,86000]
    a_list = [-0.0065, 0, 0.0010, 0.0028, 0, -0.0028,-0.0020]

    #This part of the program loops through each of the h0_list values, in order to generate a list of base layer pressures and temperatures.
    #These lists will be used later to select base values from so the calculations can be done in one step.

    T0_list = []
    p0_list = []
    T0_list.append(288.15)
    p0_list.append(101325)

    for j in range(1,len(h0_list)):
        T0_list.append(T0_list[j-1] + a_list[j-1]*(h0_list[j]-h0_list[j-1]))

        if a_list[j-1] == 0:
            p0_list.append(p0_list[j-1] * (e ** -((g / (R * T0_list[j-1])) * (h0_list[j] - h0_list[j-1]))))
        else:
            p0_list.append(p0_list[j-1] * (((T0_list[j-1] + a_list[j-1] * (h0_list[j] - h0_list[j-1])) / T0_list[j-1]) ** (-g / (a_list[j-1] * R))))

    #T0_list = [288.15, 216.650, 216.65, 228.650, 270.65, 270.65, 214.65]
    #p0_list = [101325,22632.1,5474.89,868.019,110.906,66.9389,3.95642]
    #For reference ^^

    #This part of the program uses a while loop to loop through h0_list to find the base layer altitude which corresponds to the target altitude.
    #The index number is also used to select the base values from the other lists.

    hmax = h0_list[-1]
    h = min(h,hmax)
    i = 0
    while h > h0_list[i+1]:
        i = i + 1

    h0 = h0_list[i]
    a = a_list[i]
    T0 = T0_list[i]
    p0 = p0_list[i]

    #The final pressure, temperature and density are determined by a simple if statement, since all the other values have already been selected with the use of lists.

    a_list_gradient = [0,20000,32000,51000,71000]
    a_list_isothermal = [11000,47000]

    if h0 in a_list_gradient:
        p = p0 * (((T0 + a * (h - h0)) / T0) ** (-g / (a * R)))
        T = T0 + a*(h-h0)
        rho = p / (R*T)

    if h0 in a_list_isothermal:
        p = p0 * (e ** -((g / (R * T0)) * (h - h0)))
        T = T0
        rho = p / (R*T)

    #Final ouput

    return p,rho,T











































