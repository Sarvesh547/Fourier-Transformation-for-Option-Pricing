#!/usr/bin/env python
# coding: utf-8

# In[7]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

def COSDensity(cf,x,N,a,b):
    i = np.complex(0.0,1.0) #define complex no. i=sqrt(-1)
    k = np.linspace(0,N-1,N) #define grid
    #No. of arrguments for charac fn
    u = np.zeros([1,N])
    u = k*np.pi / (b-a)
    
    #charac fn coeffients
    F_k    = 2.0 / (b - a) * np.real(cf(u) * np.exp(-i * u * a));
    F_k[0] = F_k[0] * 0.5; #first elements should multiply by half
    
    #final calculation
    f_X = np.matmul(F_k , np.cos(np.outer(u, x - a )))
    
    return f_X

def mainCalculation():
    i = np.complex(0.0, 1.0)
    #config on intergration domain
    a = -10.0
    b = 10.0
    
    #define range of no. explansion terms
    N = [16, 64, 128]
    
    #ND cfd and PDF
    mu = 0.5
    sigma = 0.2
    
     # Define characteristic function for the normal distribution
    cF = lambda u : np.exp(i * mu * u - 0.5 * sigma**2.0 * u**2.0);
    
    
    #define domain for density
    y = np.linspace(0.05,5,1000)
    
    plt.figure(1)
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("$f_Y(x)$")
    for n in N:
        f_Y = 1/y* COSDensity(cF,np.log(y),n,a,b) #calculation of density
        
        plt.plot(y,f_Y)
    plt.legend(["N=%.0f"%N[0],"N=%.0f"%N[1],"N=%.0f"%N[2]])
mainCalculation()

