"""
Gram Schmidit process of polynomial up to degree d

Created on Mon Feb 28 18:29:33 2022

@author: Edmilson Roque dos Santos
"""


import numpy as np
from scipy.integrate import quad

    
def poly(order):
    '''
    Polynomial function to be used in the moment calculation

    Parameters
    ----------
    order : int
        degree of the polynomial

    Returns
    -------
    poly_order : function

    '''
    poly_order = lambda x: x**order
    return poly_order


def l2_inner_product(func1, func2, params):
    """
    Calculate inner product between two functions func1 and func2;
    
    \langle f_1, f_2 \rangle = \int f_1 f_2 d\mu. 

    Parameters
    ----------
    func1 : function
        
    func2 : function
        
    params : dictionary

    Returns
    -------
    inner product calculation

    """
    density = params['density']
    a = params['lower_bound']
    b = params['upper_bound']
    norm_factor = params['density_normalization']
    integrand = lambda x: func1(x)*func2(x)*density(x)/norm_factor
       
    proj, err = quad(integrand, a, b, epsabs = 1e-12, epsrel = 1e-12)
    
    return np.around(proj, 8)
    #return proj


def multiply(coefficient, func):
    return lambda x: coefficient*func(x)

def difference(func1, func2):
    return lambda x: func1(x) - func2(x)


def orthonormal(order, parameters):
       
    if(order == 0):
        norm = np.sqrt(l2_inner_product(poly(0), poly(0), parameters))
        
        return lambda x: poly(0)(x)/norm
    
    else:
        temp_func = lambda x: poly(order)(x)
        
        for deg in range(1, order + 1):
            ref_func = orthonormal(deg - 1, parameters)
            proj_deg = l2_inner_product(poly(order), ref_func, parameters)
            temp_func = difference(temp_func, multiply(proj_deg, ref_func))
            
        norm = np.sqrt(l2_inner_product(temp_func, temp_func, parameters))

        normed_func = lambda x: temp_func(x)/norm
            
        return normed_func