# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def pendulo_amortecido_forcado(dt,theta_i,gamma,n_iter=2000,v=0.,g=9.81):
    t = 0    
    theta_i = (2*np.pi*theta_i)/360
    theta = theta_i
    v_theta_n = v
    w = 2*np.pi
    w0 = (3/2)*w
    B = w0/4.
    vv = [v_theta_n]
    theta_n = [theta]
    tt = [t]
    for i in range(n_iter):
        theta = theta + dt*v_theta_n
        theta_n.append(theta)
        v_theta_n = -dt*((w0**2)*np.sin(theta)+2*B*v_theta_n-gamma*(w0**2)*np.cos(w*t)) + v_theta_n
        vv.append(v_theta_n)
        t += dt
        tt.append(t)
    return (vv,theta_n,tt)
        
if __name__ == '__main__':
    name = 'pendulo amortecido forçado'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    vv1, xn1, tt1 = pendulo_amortecido_forcado(1e-2,0.,.2,6*(10**2))
    vv2, xn2, tt2 = pendulo_amortecido_forcado(1e-3,0.,.9,6*(10**3))
    vv3, xn3, tt3 = pendulo_amortecido_forcado(1e-4,0.,1.06,160000)
    vv4, xn4, tt4 = pendulo_amortecido_forcado(1e-5,0.,1.073,3*(10**6))
    fig, axs = plt.subplots(2,2)
    fig.suptitle('Espaço de Fases do Pêndulo Amortecido Forçado')
    fig.set_size_inches(10, 10)
    axs[0][0].plot(xn1[::],vv1[::],'k',linewidth=.7,label = '$\\gamma = 0.2$')
    axs[0][0].set(xlabel='$\\theta$[rad]',ylabel='velocidade[rad/s]')
    axs[0][0].grid()
    axs[0][0].legend()
    axs[0][1].plot(xn2,vv2,'k',linewidth=.7,label='$\\gamma = 0.9$')
    axs[0][1].set(xlabel='$\\theta$[rad]',ylabel='velocidade[rad/s]')
    axs[0][1].legend(loc='lower right')
    axs[1][0].plot(xn3,vv3,'k',linewidth=.7,label='$\\gamma = 1.06$')
    axs[1][0].set(xlabel='$\\theta$[rad]',ylabel='velocidade[rad/s]')
    axs[1][0].legend(loc='lower right')
    axs[1][1].plot(xn4,vv4,'k',linewidth=.7,label='$\\gamma = 1.073$')
    axs[1][1].set(xlabel='$\\theta$[rad]',ylabel='velocidade[rad/s]')
    axs[1][1].legend(loc='upper right')
    #axs[1].set_ylim([-10,10])
    plt.savefig(os.path.join('results', name, 'pendulo_amortecido_forcado_fases.png'))