# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 08:08:41 2021

@author: matheus.sartor
"""
import matplotlib.pyplot as plt
import numpy as np

def square(x):
    return 1-x**2
## Nao consegui fazer esse, ajude por favor
def gaussian(x,m,dp):
    return (((1/(np.sqrt(2*np.pi)*dp))*np.e**((-1/2)*((x-m)/dp)**2))/((1/(np.sqrt(2*np.pi)*dp))*np.e**((-1/2)*((m-m)/dp)**2)))
def exponencial(x):
    return 1-(np.e**x)/np.e


def plotMalha1D(n_points=100,x_start=0,x_end= 10.,mode="linear"):
    x = np.linspace(x_start,x_end,n_points)
    y = np.zeros(x.shape)
    xx = []
    for i in range(0,n_points):
        relative_linear_point = i/(n_points-1)
        if(mode == "quad"):
            relative_point = square(relative_linear_point)
        if(mode == "gaussian"):
            relative_point = gaussian(relative_linear_point, .5, .1)
        if(mode == "exponencial"):
            relative_point = exponencial(relative_linear_point)
        absolute_point = relative_point*x_end
        xx.append(absolute_point)
            
    plt.plot(xx,y,'.')
    if(mode == "linear"):
        return plt.plot(x,y,'.')
    #if(mode == "quad"):
     #   x = 1/(1/(np.sqrt(2*np.pi)*2))*np.e**((-1/2)*((x-(x_end/2))/2)**2)
      #  return plt.plot(x,y,'.')
    
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
    plotMalha1D(mode = "gaussian");