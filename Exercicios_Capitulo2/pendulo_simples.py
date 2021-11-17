# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def pendulo_simples(dt,l,m,theta_i,v,tmax,g=9.81):
    t = 0    
    theta_i = (2*np.pi*theta_i)/360
    theta = theta_i
    v_theta_n = v
    vv = [v_theta_n]
    theta_n = [theta]
    theta_a = [theta]
    tt = [t]
    while t <= tmax:
        theta = theta + dt*v_theta_n
        theta_n.append(theta)
        v_theta_n = -dt*g*np.sin(theta)/l + v_theta_n
        vv.append(v_theta_n)
        theta_an = theta_i*np.cos(np.sqrt(g/l)*t)
        theta_a.append(theta_an)
        t += dt
        tt.append(t)
    return (vv,theta_n,theta_a,tt)
        
if __name__ == '__main__':
    name = 'pendulo simples'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    vv, xn, xa, tt = pendulo_simples(0.01,1.,1.,10,0.,10)
    fig, axs = plt.subplots(2)
    fig.suptitle('Análise do Pêndulo Simples')
    fig.set_size_inches(10, 10)
    axs[0].plot(tt,xa,'k',linewidth=2,label='Solução Analítica')
    axs[0].plot(tt[::20],xn[::20],'ko',linewidth=2,label = 'Solução Numérica')
    axs[0].set(xlabel='tempo [s]',ylabel='$\\theta$[rad]')
    axs[0].grid()
    axs[0].legend()
    axs[1].plot(xn,vv,'k',linewidth=2)
    axs[1].set(xlabel='$\\theta$[rad]',ylabel='velocidade [rad/s]')
    axs[1].set_ylim(-0.6,.6)
    axs[1].set_xlim(-0.2,.2)
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'pendulo_simples.png'))