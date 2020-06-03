# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 19:36:57 2020

@author: Casper
"""
import numpy as np 

a = [[1,2,3],[3,4,5],[2,3,1]]
b = [[1],[1],[1]]

x = np.linalg.solve(a,b)
