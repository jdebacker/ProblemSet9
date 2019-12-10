# import packages
import numpy as np
import matplotlib.pyplot as plt
import SS
import TPI

# set model parameters
sigma = 1.5  # CRRA
beta = 0.8  # discount rate
alpha = 0.3   # capital share of output
delta = 0.1  # rate of depreciation
A = 1.0  # TFP
T = 20  # number of periods until SS

# parameter for convergence of GE loop
xi = 0.1

# Solving for SS
'''
Provided a guess of w = 0.05. Request to cross check if there is 
any better guess that can be used

We still need to add the TPI file and then update its execution here
'''

ss_params = (beta, sigma, alpha, A, delta, xi)
r_init = 1 / beta - 1
w_init = 0.5 # Not sure if there is a possible way to guess
r_ss, w_ss, b_sp1_ss, euler_errors_ss = SS.solve_ss(r_init, w_init, ss_params)
print('SS interest rate is ', r_ss)
print ('SS wage rate is ', w_ss)
print('Maximum Euler error in the SS is ',
      np.absolute(euler_errors_ss).max())


