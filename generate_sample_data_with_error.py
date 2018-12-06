# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:41:09 2018

@author: al064648
"""

import numpy as np
import random


fname_lin = 'C:\\Users\\al064648\\source\\repos\\Least_squares_method\\Least_squares_method\\linear_test.txt'
fname_quad = 'C:\\Users\\al064648\\source\\repos\\Least_squares_method\\Least_squares_method\\quadratic_test.txt'
fname_cub = 'C:\\Users\\al064648\\source\\repos\\Least_squares_method\\Least_squares_method\\cubic_test.txt'

N = 1000
coefs = [4., 3., 2., 1.]

x = np.linspace(-100., 100, num = N)

with open(fname_lin,'w') as f:
    f.write('linear: y=a*x+b, a={}, b={} \n'.format(coefs[1], coefs[0]))
    f.write('{} \n'.format(N))
    for i in range(N):
        y = 0.0
        for ord in range(2):
            y += coefs[ord] * x[i]**ord
        y += random.random()
        f.write('{} {} \n'.format(x[i],y))

with open(fname_quad,'w') as f:
    f.write('quadratic: y=a*x*x + b*x + c, a={}, b={}, c={} \n'.format(coefs[2], coefs[1], coefs[0]))
    f.write('{} \n'.format(N))
    for i in range(N):
        y = 0.0
        for ord in range(3):
            y += coefs[ord] * x[i]**ord
        y += random.random()
        f.write('{} {} \n'.format(x[i],y))
        
with open(fname_cub,'w') as f:
    f.write('cubic: y=a*x^3 + b*x^2 + c*x + d, a={}, b={}, c={}, d={} \n'.format(coefs[3], coefs[2], coefs[1], coefs[0]))
    f.write('{} \n'.format(N))
    for i in range(N):
        y = 0.0
        for ord in range(4):
            y += coefs[ord] * x[i]**ord
        y += random.random()
        f.write('{} {} \n'.format(x[i],y))

#for i in range(1,5):
#    for j in range(1,5):
#        print('A({},{}) = sum(x**({}))'.format(i,j,i+j-2))
#    print('B({}) = sum(y*x**({}))'.format(i,i-1))