#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 18 10:02:13 2024

@author: demroz
"""

import numpy as np
import pandas as pd
N = 1000
M = 2
x = np.random.uniform(size=(1000,M))
y = 0.5 + 0.1 * x[:,0] + 0.25 * x[:,1] + np.random.uniform(low=-0.1,high=0.1,size=(N,))

d = "/home/demroz/Documents/code/basic_math/tests/test_matrices/"

data = {}
data['x1'] = x[:,0]
data['x2'] = x[:,1]
data['y'] = y

df = pd.DataFrame(data)
df.to_csv(d+'regression_test_data.csv')

X = np.zeros((N,M+1))
X[:,0]=1
X[:,1:M+1] = x

coefficients, _, _, _ = np.linalg.lstsq(X, y, rcond=None)

import scipy as sp
slope, intercept, r_value, p_value, std_err = sp.stats.linregress(x, y)