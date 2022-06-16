#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as fft
import scipy.interpolate as interpolate


# In[23]:


def RecoveryDensity(cf,c,N= 8192):  # if 2**8 is empty then it will take 8192 i.e 8^13
    i = np.complex(0.0,1.0) #assigning -=sqrt(-1) we have only imaginery unit
    
    #specification of the grid for u 
    u_max = 20.0  # umax of our intergral it will choice when spining of charac fn
    du = u_max / N #grid size Δu
    u = np.linspace(0,N-1,N) * du  #divide points grid from min to max
    
    #specification of the grid for x
    b   = np.min(x)
    dx  = 2.0*np.pi / (N*du) #use relation of dx*du
    x_i = b + np.linspace(0,N-1,N) * dx #grid
    
    #chrac fn i.e we have phi as input chrac fn
    phi = np.exp(-i*b*u) * cf(u) #we have applied fourir trandformation
    #we applied first and last point to boundary
    gamma_1 = np.exp(-i*x_i*u[0])*cf(u[0]) 
    gamma_2 = np.exp(-i*x_i*u[-1])*cf(u[-1])
                     
    phi_boundary = 0.5 * (gamma_1 + gamma_2)   #it is boundarys of trapozal integration
    
    #FFT at this point we will recover density at grid points waht we calculates in x_i
    f_xi = du/np.pi * np.real(fft.fft(phi)- phi_boundary)
     
    #Interolation making sure input x is intergrate
    f_xiInterp = interpolate.interp1d(x_i,f_xi,kind='cubic') 
                     
    return f_xiInterp(x)  

def mainCalculation():
    i = np.complex(0.0,1.0)  # assigning i=sqrt(-1)
    
    #setting ND
    mu= 0.0
    sigma = 1.0
    
    #define charc fn for the ND
    cF = lambda u : np.exp(i * mu * u - 0.5 * sigma**2.0 * u**2.0);
        
    #define domain for density
    x = np.linspace(-8.0,8.0,100) #grid
    f_XExact = st.norm.pdf(x,mu,sigma)  #we have define pdf from libary
    
    
    f_XR = RecoverDensity(cF,x,2**8) #2**8 is no. of points i.e N
    
    plt.figure(1)
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("$f_X(x)$")
    plt.plot(x,f_XExact,'-r')
    plt.plot(x,f_XR,'--b')
    plt.legend(['Exact PDF','Approximated PDF'])
    
mainCalculation()


# In[ ]:





# In[ ]:




