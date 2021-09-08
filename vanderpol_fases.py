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
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.title('Espaço de Fases da Equação de Van Der Pol - Runge-Kutta')
    plt.plot(xn1,vv1,'k',linewidth=.7,label='$\\mu = 0.1$')
    plt.xlabel='x[m]'
    plt.ylabel='velocidade[m/s]'
    plt.ylim(-4,4)
    plt.xlim(-4,4)
    plt.plot(xn2,vv2,'b',linewidth=.7,label='$\\mu = 0.5$')
    plt.plot(xn3,vv3,'r',linewidth=.7,label='$\\mu = 1.0$')
    plt.plot(xn4,vv4,'g',linewidth=.7,label='$\\mu = 1.5$')
    plt.legend(loc='lower right')
    ax.set_aspect('equal', adjustable='box')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'analise_van-Der-Pol_fases_Runge-Kutta.png'))