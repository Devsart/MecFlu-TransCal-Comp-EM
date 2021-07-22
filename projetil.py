# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 20:14:11 2021

@author: matheus.sartor
"""
import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np



def lancamento(x,y,theta,dt,xmax,m,beta,gamma,D,v,mode = "linear", g= 9.81):
    t=0
    b = beta*D
    c = gamma*D**2
    a = np.radians(theta)
    vx = v*np.cos(a)
    xx = [x]
    vy = v*np.sin(a)
    yy = [y]
    while x <= xmax:
        if(mode == "linear"):
            fdrag = -b*v
        elif(mode == "quadratico"):
            fdrag = -c*v**2
        elif(mode == "sem atrito"):
            fdrag = 0
        x = x + dt*vx
        xx.append(x)
        vx = dt*fdrag*np.cos(a)/m + vx
        y = y + dt*vy
        yy.append(y)
        vy = dt*(-g+fdrag*np.sin(a)/m) + vy
        if(y <= 0):
            vy = -vy
        v = np.sqrt(vx**2+vy**2)
        t += dt
        a = np.arctan(vy/vx)
    
    return (xx,yy)
        
if __name__ == '__main__':
    name = 'projetil'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    (xx_,yy_) = lancamento(0,.40,50,0.01,200,0.15,0.16,0.25,7e-2,30,"sem atrito")
    (xx_lin,yy_lin) = lancamento(0.,.40,50,0.01,200.,.15,.16,.25,7e-2,30)
    (xx_quad,yy_quad) = lancamento(0.,.40,50,0.01,200.,.15,.16,.25,7e-2,30,"quadratico")
    
    plt.title('Análise do projétil')
    plt.plot(xx_,yy_,linewidth=2, label = 'vácuo')
    plt.plot(xx_lin,yy_lin,'--',linewidth=2, label = 'atrito linear')
    plt.plot(xx_quad,yy_quad,'-.',linewidth=2, label = 'atrito quadrático')
    plt.ylabel('y[m]')
    plt.xlabel('x[m]')
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join('results', name, 'analise_projetil.png'))