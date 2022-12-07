#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 10:58:02 2022

@author: biancaboi
"""

import matplotlib.pyplot as plt
import scipy.integrate as spi
import numpy as np


# y = [V, I, F, Cn, E, Bn, P, Am]

#In this version we suppress tau_c znd tau_b which will be egal to 0

V0 = 1e4
T = 7e7
#teta_b = teta_c =0

def f(y, t, pv=210, delta_v=5, km=0.6, beta=5e-7, beta_prim=3e-8, rho=2.6, fi=0.33, delta_i=2, kn=2.5,
      ke=5e-5, pf=1e-5, delta_f=2, beta_cn=1, hc=1e4, pc=1.2, delta_e=0.57, beta_bn=0.03, hb=1e4,
      delta_p=0.5, pb=0.52, pm=8, delta_m=1.0075):
    
    z = np.zeros(np.size(y))
    
    z[0] = pv * y[1] - delta_v * y[0] - km * y[8] - beta * y[0] * T

    z[1] = beta_prim * y[0] * T - delta_i * y[1] - kn * y[1] * y[2] - ke * y[1] * y[5]
    
    z[2] = pf * y[1] - delta_f * y[2]
    
    z[3] = fi * y[2] * T - rho * y[3]
    
    z[4] = - beta_cn * (y[0] / (y[0] + hc)) * y[4]
    
    z[5] = beta_cn * (y[0]  / (y[0] + hc)) * y[4]  - delta_e * y[5]
    
    z[6] = - beta_bn * (y[0] / (y[0] + hb)) * y[6]
    
    z[7] = beta_bn * (y[0]  / (y[0] + hb)) * y[6]  - delta_p * y[7]
    
    z[8] = pm * y[7] - delta_m * y[8]

    return z
leg=['V', 'I', 'F', 'R', 'Cn', 'E', 'Bn', 'P', 'Am']

def afficher():
    y0 = np.array([V0, 0, 0, 0, 100, 0, 100, 0, 0])
    t0 = 0
    T = 100
    pas = 0.01
    tps = np.arange(t0, T, pas)

    y = spi.odeint(f, y0, tps)

    #for u in range(len(y0)):
    #plt.plot(tps, y[:, 0],label=leg[0])
    plt.plot(tps, y[:,4], label=leg[4])
    plt.plot(tps, y[:,6],label=leg[6]) #no change ... problem
    plt.legend()
    plt.show()


afficher()