# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def massamola(dt,k,m,v,tmax):
    t = 0
    x = 0
    
    vxn = v
    vv = [vxn]
    va = [vxn]
    xn = [x]
    xa = [x]
    tt = [t]
    w = np.sqrt((k/m))
    j = 0
    A = v/(-w*np.cos(w*t+j))
    B = v/(w**2*np.cos(w*t+j))
    while t <= tmax:
        x = x + dt*vxn
        xn.append(x)
        vxn = -dt*k*x/m + vxn
        vxa = -w*A*np.cos(w*t+j)
        vv.append(vxn)
        va.append(vxa)
        xxa = w*B*np.sin(w*t + j)
        xa.append(xxa)
        t += dt
        tt.append(t)
    mape_v = np.mean(np.abs(np.array(va).reshape(len(va),1) - np.array(vv).reshape(len(vv),1)) / np.abs(np.array(va).reshape(len(va),1) + 1e-20))
    return (vv,va,xn,xa,tt, mape_v)
        
if __name__ == '__main__':
    name = 'massamola'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    vv, va, xn, xa, tt, mape_v = massamola(0.01,0.1,1,10,100)
    fig, axs = plt.subplots(2)
    fig.suptitle('Análise do Sistema Massa Mola')
    fig.set_size_inches(10, 10)
    axs[0].plot(tt,va,'k',linewidth=2,label='Solução Analítica')
    axs[0].plot(tt[::200],vv[::200],'ko',linewidth=2,label = 'Solução Numérica')
    axs[0].set(xlabel='t[s]',ylabel='v[m/s]')
    axs[0].set_ylim([-20,20])
    axs[0].grid()
    axs[0].legend(title=f'Erro MAPE: {round(100 * mape_v, 2)} %',loc='upper right')
    axs[1].plot(xa,va,'k',linewidth=2)
    axs[1].set(xlabel='x[m]',ylabel='v[m/s]')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'analise_massamola.png'))