# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 08:21:39 2021

@author: matheus.sartor
"""

# Primeiro exercício: calculo da curva de velocidade

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

name = 'carrinho'
Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
m = 1.0
b = 0.1
v0 = 10.0
dt = 0.01
t = 0
x=0
x0 = 0
vx = v0
vn = [v0]
va = [v0]
tt = [t]
xx=[x0]
xa=[x0]

while t <= 70:
    vx = (-dt*b/m + 1)*vx
    vn.append(vx)
    x = x + dt*vx
    xx.append(x)
    t += dt
    tt.append(t)
    van = v0*np.e**(-b*t/m)
    va.append(van)
    xan = (v0*m/b)*(1-np.e**(-b*t/m))
    xa.append(xan)

mape_u = np.mean(np.abs(np.array(va).reshape(len(va),1) - np.array(vn).reshape(len(vn),1)) / np.abs(np.array(va).reshape(len(va),1) + 1e-20))
mape_x = np.mean(np.abs(np.array(xa).reshape(len(xa),1) - np.array(xx).reshape(len(xx),1)) / np.abs(np.array(xa).reshape(len(xa),1) + 1e-20))
fig, axs = plt.subplots(2)
fig.suptitle('Análise do Movimento do Carrinho')
axs[0].plot(tt,va,linewidth=3, label = 'Solução Analítica')
axs[0].plot(tt[::200],vn[::200],'o',alpha = 0.5,linewidth=3, label = 'Solução Numérica')
axs[0].set(xlabel='t[s]',ylabel='v[m/s]')
axs[0].grid()
axs[0].legend(title=f'Erro MAPE: {round(100 * mape_u, 2)} %')
axs[1].legend(title=f'Erro MAPE: {round(100 * mape_x, 2)} %')
axs[1].plot(tt,xa,linewidth=3, label = 'Solução Analítica')
axs[1].plot(tt[::200],xx[::200],'o',alpha = 0.5,linewidth=3, label = 'Solução Numérica')
axs[1].set(xlabel='t[s]',ylabel='x[m]')
axs[1].grid()
plt.savefig(os.path.join('results', name, 'analise_carrinho.png'))



print(x)
plt.show()