#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 19 13:04:59 2020

@author: atkachev
"""

from magnetic import Ap, Bz
import numpy as np
from scipy import integrate
import time

R = np.linspace(0.001, 0.499, 10001)

t1 = time.time()
A1 = np.array([Ap(r,0,0.5,1000) for r in R])
t2 = time.time()
print(t2 - t1)

# f = lambda r: Bz(r,0,0.5,1000)*r
# t1 = time.time()
# A2 = np.array([integrate.quad(f,0,r)/r for r in R])[:,0]
# t2 = time.time()
# print(t2 - t1)

# err = A2 - A1

