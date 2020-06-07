# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 17:50:06 2020

@author: Casper
"""

import statsmodels.api as sm 
import numpy as np 


x = np.array([1,2,3,4,5]) 
y = x/2 

x1 = sm.add_constant(x)

results = sm.OLS(y,x1).fit()
