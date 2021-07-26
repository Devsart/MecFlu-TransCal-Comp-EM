# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 14:06:44 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def massamola_vertical(dt,k,m,v,tmax):
    t = 0
    y = 0
    g = 9.81
    vyn = v
    vv = [vyn]
    yn = [y]
    tt = [t]
    while t <= tmax:
        vyn = dt*(-k*y-m*g)/m + vyn
        vv.append(vyn)
        y = y + dt*vyn
        yn.append(y)
        t += dt
        tt.append(t)
    return (vv,yn,tt)
        
if __name__ == '__main__':
    name = 'massamola'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    vv, yn, tt = massamola_vertical(0.1,0.1,0.15,10,30)
    fig, axs = plt.subplots(2)
    fig.suptitle('Análise do Sistema Massa Mola Vertical')
    fig.set_size_inches(10, 10)
    axs[0].plot(tt[::],vv[::],'k.',linewidth=2,label = 'Solução Numérica')
    axs[0].set(xlabel='t[s]',ylabel='v[m/s]')
    axs[0].set_ylim([-20,20])
    axs[0].legend(loc='lower right')
    axs[0].grid()
    axs[1].plot(yn,vv,'k',linewidth=2)
    axs[1].set(xlabel='y[m]',ylabel='v[m/s]')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'analise_massamola_vertical.png'))
