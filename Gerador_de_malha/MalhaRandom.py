# -*- coding: utf-8 -*-
"""
Created on Wed Jul 28 23:14:35 2021

@author: matheus.sartor
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import matplotlib

       
def getXYRandom(nx=6,ny=7,Lx = 1, Ly = 1):
    X = np.random.uniform(0.0,Lx,(nx*ny,1))
    X = (np.transpose(X.reshape(nx*ny)))

    Y = np.random.uniform(0.0,Ly,(nx*ny,1))
    Y = (np.transpose(Y.reshape(nx*ny)))
    Tri = matplotlib.tri.Triangulation(X,Y)
    return X,Y,Tri

def plotMalhaTri(X,Y,Tri):
    name = 'malhaTriangularRandom'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    IEN = Tri.triangles
    ne = IEN.shape[0]
    npoints = X.shape[0]
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
        crit_list.append(crit)
        
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    ax.set_aspect('equal')
    ax.set_title("Malha de Triângulos")

    ax = plt.triplot(X,Y,IEN,color='k',linewidth=0.5)
    ax = plt.plot(X,Y,'ko',markersize=2)
    for i in range(0,npoints):
        plt.text(X[i]+0.01,Y[i]+0.01,str(i),size='x-small',color='k')
    for e in range(0,ne):
        [v1,v2,v3] = IEN[e]
        xm = ( X[v1]+X[v2]+X[v3] )/3.0
        ym = ( Y[v1]+Y[v2]+Y[v3] )/3.0
        draw_circle = plt.Circle((xm+0.01, ym+0.01), 0.025, fill=False, color='gray')
        plt.text(xm,ym,str(e),size='xx-small',color='gray')
        plt.gcf().gca().add_artist(draw_circle)

    plt.gca().set_aspect('equal')
    plt.tripcolor(Tri, crit_list, edgecolors='k',cmap='Blues',vmin=0,vmax=1)
    plt.colorbar()
    plt.savefig(os.path.join('results', name, 'Malha_de_Triangulos.png'))


if __name__ == '__main__':
    X,Y,Tri = getXYRandom()
    plotMalhaTri(X,Y,Tri)