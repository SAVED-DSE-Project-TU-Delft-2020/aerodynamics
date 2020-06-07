# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 19:36:57 2020

@author: Casper
"""
import numpy as np 

a = np.array([[1],[1],[1]])
b = np.array([[2],[2],[2]])

assert np.allclose(a,b,rtol = 0.0001)