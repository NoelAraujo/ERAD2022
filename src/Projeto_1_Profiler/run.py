"""
Gram-Schmidt process script

Created on Mon Feb 28 18:32:25 2022

@author: Edmilson Roque dos Santos
"""
import numpy as np
import matplotlib.pyplot as plt

import time
import GS_process as gs

# Set plotting parameters
params_plot = {'axes.labelsize': 16,
              'axes.titlesize': 18,
              'axes.linewidth': 1.0,
              'axes.xmargin':0, 
              'axes.ymargin': 0,
              'legend.fontsize': 18,
              'xtick.labelsize': 16,
              'ytick.labelsize': 16,
              'figure.figsize': (8, 6),
              'figure.titlesize': 18,
              'font.serif': 'Computer Modern Serif',
              'mathtext.fontset': 'cm',
              'axes.linewidth': 1.0
             }

plt.rcParams.update(params_plot)
plt.rc('text', usetex=True)


parameters = dict()
parameters['density'] = lambda x: 1
parameters['lower_bound'] = -1
parameters['upper_bound'] = 1
parameters['density_normalization'] = 1


for deg in range(1,15):
    start = time.process_time()
    orth_1 = gs.orthonormal(deg, parameters)    #orthonormal function
    print(deg, "   ",time.process_time() - start)

print_inner_product = True
if print_inner_product: #to check if the functions are orthonormal
    deg = 13         #degree
    
    start = time.process_time()
    orth_2 = gs.orthonormal(deg, parameters)    #orthonormal function
    print(time.process_time() - start)
    
    gs.l2_inner_product(orth_1, orth_2, parameters)
    



plot_orth = True
if plot_orth:
    x = np.arange(parameters['lower_bound'], parameters['upper_bound'], 0.01)
    plt.plot(x, orth_1(x))
    plt.plot(x, orth_2(x))
    plt.xlabel(r'$x$')
    plt.ylabel(r'$\varphi(x)$')