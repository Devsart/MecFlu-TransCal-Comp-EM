# -*- coding: utf-8 -*-
"""
Created on Mon Jul 26 11:40:06 2021

@author: matheus.sartor
"""

import os
from mpl_toolkits import mplot3d
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def borboleta_lorenz(dt,sigma,ro,beta,n_iter):
    t = 0.
    x = 15.
    y = -25.
    z = -20.
    
    xn = [x]
    yn = [y]
    zn = [z]
    tt = [t]
    for i in range(n_iter):
        x_novo = x + dt*(sigma*(y - x))
        y_novo = y + dt*(x*(ro - z) - y)
        z_novo = z + dt*(x*y - beta*z)
        x = x_novo
        y = y_novo
        z = z_novo
        xn.append(x)
        yn.append(y)
        zn.append(z)
        t += dt
        tt.append(t)
    return (xn,yn,zn,tt)

def borboleta_lorenz_runge_kutta(dt,m,mu,v,n_iter):
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
    name = 'borboleta_lorenz'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    x, y, z, tt = borboleta_lorenz(1e-4,10,28,8/3,10**6)
    fig = plt.figure()
    fig.suptitle("Análise 3D Borboleta de Lorenz")
    ax = plt.axes(projection='3d')
    ax.plot3D(x, y, z, 'k', linewidth=.2)
    plt.savefig(os.path.join('results', name, 'analise_borboleta_lorenz.png'))