# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def vanDerPol(dt,m,mu,v,tmax):
    t = 0
    x = 0.2
    
    vxn = v
    vv = [vxn]
    xn = [x]
    tt = [t]
    while t <= tmax:
        x = x + dt*vxn
        xn.append(x)
        vxn = dt*(mu*(1-x**2)*vxn-x) + vxn
        vv.append(vxn)
        t += dt
        tt.append(t)
    return (vv,xn,tt)

def vanDerPolRungeKutta(dt,m,mu,v,n_iter):
    t = 0
    x = 0.2
    
    vxn = v
    vv = [vxn]
    xn = [x]
    tt = [t]
    for i in range(n_iter):
        x = x + dt*vxn
        xn.append(x)
        k1 = dt*mu*(1 - x**2)*vxn - dt*x
        k2 = dt*mu*(1 - x**2)*(vxn + k1/2) - dt*x
        vxn = (1/2)*(k1+k2) + vxn
        vv.append(vxn)
        t += dt
        tt.append(t)
    return (vv,xn,tt)
        
if __name__ == '__main__':
    name = 'van_der_Pol'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    vv1, xn1, tt = vanDerPolRungeKutta(0.05,1.,.1,.10,2000)
    vv2, xn2, __ = vanDerPolRungeKutta(0.05,1.,.5,.10,2000)
    vv3, xn3, __ = vanDerPolRungeKutta(0.05,1.,1,.10,2000)
    vv4, xn4, __ = vanDerPolRungeKutta(0.05,1.,1.5,.10,2000)
    fig, axs = plt.subplots(2,2)
    fig.suptitle('Análise da Equação de Van Der Pol - Runge-Kutta')
    fig.set_size_inches(10, 7)
    axs[0][0].plot(tt,vv1,'k',linewidth=.7,label='$\\mu = 0.1$')
    axs[0][0].set(xlabel='t[s]',ylabel='velocidade[m/s]')
    axs[0][0].legend(loc='lower right')
    axs[0][1].plot(tt,vv2,'k',linewidth=.7,label='$\\mu = 0.5$')
    axs[0][1].set(xlabel='t[s]',ylabel='velocidade[m/s]')
    axs[0][1].legend(loc='lower right')
    axs[1][0].plot(tt,vv3,'k',linewidth=.7,label='$\\mu = 1.0$')
    axs[1][0].set(xlabel='t[s]',ylabel='velocidade[m/s]')
    axs[1][0].legend(loc='lower right')
    axs[1][1].plot(tt,vv4,'k',linewidth=.7,label='$\\mu = 1.5$')
    axs[1][1].set(xlabel='t[s]',ylabel='velocidade[m/s]')
    axs[1][1].legend(loc='lower right')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'analise_van-Der-Pol_RungeKutta.png'))