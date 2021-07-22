# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 08:08:41 2021

@author: matheus.sartor
"""
import matplotlib.pyplot as plt
import numpy as np


def plotMalha1D(n_points=50,x_start=0,x_end= 10.,mode="linear"):
    x = np.linspace(x_start,x_end,n_points)
    y = np.zeros(x.shape)
    if(mode == "linear"):
        return plt.plot(x,y,'.')
    if(mode == "quad"):
        x = np.e**(-x)
        return plt.plot(x,y,'.')
    
# matriz de conectividade IEN
def IENcreate(ne):
    IEN = np.zeros((ne,2))
    for e in range(0,ne):
        IEN[e] = [e,e+1]
    print(IEN)
        


def plotMalha2D(Lx = 1, Ly = 1, n_pointsx=4,n_pointsy=4):
    x = []
    y = []
    for i in range(0,n_pointsy):
        for j in range(0,n_pointsx):
            x.append(j/(n_pointsx-1))
            y.append(i/(n_pointsy-1))
    plt.plot(x,y,'o')


if __name__ == '__main__':
    plotMalha2D();