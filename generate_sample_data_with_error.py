# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:41:09 2018

@author: al064648
"""

import numpy as np
import random

fname_lin = './linear_test.txt'
fname_quad = './quadratic_test.txt'
fname_cub = './cubic_test.txt'

N = 1000
coefs = [4., 3., 2., 1.]

x = np.linspace(-10., 10, num = N)

def write_block(ff, nsize, order, coefs, x):
    ff.write('{} \n'.format(nsize))
    for i in range(nsize):
        y = 0.0
        for ord in range(order+1):
            y += coefs[ord] * x[i]**ord
        y += random.random()
        ff.write('{} {} \n'.format(x[i],y))
    
with open(fname_lin,'w') as f:
    f.write('linear: y=a*x+b, a={}, b={} \n'.format(coefs[1], coefs[0]))
    write_block(f, N, 1, coefs, x)

with open(fname_quad,'w') as f:
    f.write('quadratic: y=a*x*x + b*x + c, a={}, b={}, c={} \n'.format(coefs[2], coefs[1], coefs[0]))
    write_block(f, N, 2, coefs, x)
        
with open(fname_cub,'w') as f:
    f.write('cubic: y=a*x^3 + b*x^2 + c*x + d, a={}, b={}, c={}, d={} \n'.format(coefs[3], coefs[2], coefs[1], coefs[0]))
    write_block(f, N, 3, coefs, x)
