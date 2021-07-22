# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 15:16:19 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

name = 'gota'
Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)

# Atrito linear de uma gota caindo

ro = 840.
g = 9.81
D = 1.5e-6
V = np.pi*(D**3)/6.
beta = 1.6e-4
dt = 1e-7
v = 0.
vn = [v]
va = [v]
b = beta*D
m = ro*V
t = 0
tt = [t]
tlim = 5.0e-5


while t <= tlim:
    v = dt*(g-b*v/m) + v
    t += dt
    vn.append(v)
    tt.append(t)
    van = va[0]*np.e**(-b*t/m)-(-m*g/b)*(1-np.e**(-b*t/m))
    va.append(van)

mape_u = np.mean(np.abs(np.array(va).reshape(len(va),1) - np.array(vn).reshape(len(vn),1)) / np.abs(np.array(va).reshape(len(va),1) + 1e-20))
# mape_x = np.mean(np.abs(np.array(xa).reshape(len(xa),1) - np.array(xx).reshape(len(xx),1)) / np.abs(np.array(xa).reshape(len(xa),1) + 1e-20))

plt.title('Análise da Velocidade da Gota')
plt.plot(tt,va,linewidth=3, label = 'Solução Analítica')
plt.plot(tt[::20],vn[::20],'o',alpha = 0.5,linewidth=3, label = 'Solução Numérica')
plt.xlabel('t[s]')
plt.ylabel('v[m/s]')
plt.grid()
plt.legend(title=f'Erro MAPE: {round(100 * mape_u, 2)} %')
# axs[1].legend(title=f'Erro MAPE: {round(100 * mape_x, 2)} %')
# axs[1].plot(tt,xa,linewidth=3, label = 'Solução Analítica')
# axs[1].plot(tt[::200],xx[::200],'o',alpha = 0.5,linewidth=3, label = 'Solução Numérica')
# axs[1].set(xlabel='t[s]',ylabel='x[m]')
# axs[1].grid()
plt.savefig(os.path.join('results', name, 'analise_gota.png'))
