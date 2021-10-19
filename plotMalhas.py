# -*- coding: utf-8 -*-
"""
Created on Thu Jul 22 08:08:41 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

def square(x):
    return 1-x**2
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
        


def getXY(n_pointsx=4,n_pointsy=4,Lx = 1, Ly = 1, mode="square",function="None",coef=.2):
    x = []
    y = []
    for i in range(0,n_pointsy):
        for j in range(0,n_pointsx):
            x.append(j/(n_pointsx-1))
            y.append(i/(n_pointsy-1))
    IEN = IEN2D(n_pointsx,n_pointsy,mode)
    if(function=="sin"):
        y[-n_pointsx::] += y[-n_pointsx::]*np.sin(np.multiply(x[:n_pointsx:],2*np.pi))*coef
        y[-2*n_pointsx:-n_pointsx:] += y[-2*n_pointsx:-n_pointsx:]*np.sin(np.multiply(x[:n_pointsx:],2*np.pi))*coef
        y[-3*n_pointsx:-2*n_pointsx:] += y[-3*n_pointsx:-2*n_pointsx:]*np.sin(np.multiply(x[:n_pointsx:],2*np.pi))*coef
        y[-4*n_pointsx:-3*n_pointsx:] += y[-4*n_pointsx:-3*n_pointsx:]*np.sin(np.multiply(x[:n_pointsx:],2*np.pi))*coef
    x_array = np.asarray(x)
    x_array = (np.transpose(x_array.reshape(n_pointsx*n_pointsy)))
    y_array = np.asarray(y)
    y_array = (np.transpose(y_array.reshape(n_pointsx*n_pointsy)))
    IEN_array = np.asarray(IEN)
    return x_array,y_array,IEN_array

def plotMalhaQuadrada(X,Y,IEN):
    name = 'Malha Ordenada'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    ne = IEN.shape[0]
    crit_list = []
    for e in range(0,ne):
        [v1,v2,v3] = IEN[e]
        x1 = X[v1]; x2 = X[v2]; x3 = X[v3]
        y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3]
        a = np.sqrt((x3-x2)**2+(y3-y2)**2)
        b = np.sqrt((x3-x1)**2+(y3-y1)**2)
        c = np.sqrt((x2-x1)**2+(y2-y1)**2)
        s = (a + b + c)/2
        A = np.sqrt(s*(s-a)*(s-b)*(s-c)) # Formula de Heron
        ha = 2*A/a
        hb = 2*A/b
        hc = 2*A/c
        crit = ((ha*hb*hc)/(a*b*c))/((np.sqrt(3)/2)**3)
        crit_list.append((0,0,.6,crit))
    xy = np.stack((X, Y), axis=-1)
    verts = xy[IEN]
    ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts,edgecolors=('black',),
                                                 facecolor='None',
                                                 linewidths=(0.7,))
    pc.set_color(crit_list)
    ax.add_collection(pc)
    plt.plot(X,Y,'ko')
    plt.savefig(os.path.join('results', name, 'malha_senoidal_colorida.png'))

def IEN2D(nx,ny, mode="square"):
    if(mode=="square"):
        ne = (nx-1)*(ny-1)
        IEN = np.zeros((ne,4),dtype=int)
        i = 0
        for e in range(0,ne):
            if(e!=0 and ((e-1)%(nx-1)) == (nx - 2)):
                IEN[e] = [i+1,i+2,i+nx+2,i+nx+1]
                i += 2
            else:
                IEN[e] = [i,i+1,i+nx+1,i+nx]
                i += 1
    elif(mode=="triangle"):
        ne = (nx-1)*(ny-1)*2
        IEN = np.zeros((ne,3),dtype=int)
        i = 0
        for e in range(0,ne,2):
            if(e!=0 and (((e/2) -1)%((nx-1))) == (nx - 2)):
                IEN[e] = [i+1,i+nx+2,i+nx+1]
                IEN[e+1] = [i+1,i+nx+2,i+2]
                i += 2
            else:
                IEN[e] = [i,i+nx+1,i+nx]
                IEN[e+1] = [i,i+nx+1,i+1]
                i += 1
    return IEN
if __name__ == '__main__':
    X,Y,IEN = getXY(mode="triangle", function="sin")
    plotMalhaQuadrada(X,Y,IEN)