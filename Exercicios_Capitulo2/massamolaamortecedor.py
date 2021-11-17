# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def massamolaamortecedor(dt,m,w,B,v,tmax):
    t = 0
    x = 0
    
    vxn = v
    vv = [vxn]
    xn = [x]
    tt = [t]
    while t <= tmax:
        x = x + dt*vxn
        xn.append(x)
        vxn = dt*(-2*B*vxn-(w**2)*x) + vxn
        vv.append(vxn)
        t += dt
        tt.append(t)
    return (vv,xn,tt)
        
if __name__ == '__main__':
    name = 'massamola'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    __, xn1, tt = massamolaamortecedor(0.1,1.,.5,.0,10,35)
    __, xn2, __ = massamolaamortecedor(0.1,1.,.5,.1,10,35)
    __, xn3, __ = massamolaamortecedor(0.1,1.,.5,.75,10,35)
    __, xn4, __ = massamolaamortecedor(0.1,1.,.5,.5,10,35)
    fig, axs = plt.subplots(2,2)
    fig.suptitle('Análise do Sistema Massa Mola')
    fig.set_size_inches(15, 10)
    axs[0][0].plot(tt,xn1,'k.',linewidth=.5,label='$\\beta = 0$')
    axs[0][0].set(xlabel='t[s]',ylabel='posicao[m]')
    axs[0][0].set_ylim([-25,25])
    axs[0][0].legend(loc='lower right')
    axs[0][1].plot(tt,xn2,'k.',linewidth=.5,label='$\\beta < \\omega$')
    axs[0][1].set(xlabel='t[s]',ylabel='posicao[m]')
    axs[0][1].set_ylim([-25,25])
    axs[0][1].legend(loc='lower right')
    axs[1][0].plot(tt,xn3,'k.',linewidth=.5,label='$\\beta > \\omega$')
    axs[1][0].set(xlabel='t[s]',ylabel='posicao[m]')
    axs[1][0].set_ylim([-25,25])
    axs[1][0].legend(loc='lower right')
    axs[1][1].plot(tt,xn4,'k.',linewidth=.5,label='$\\beta = \\omega$')
    axs[1][1].set(xlabel='t[s]',ylabel='posicao[m]')
    axs[1][1].set_ylim([-25,25])
    axs[1][1].legend(loc='lower right')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'analise_massamola_amortecedor.png'))