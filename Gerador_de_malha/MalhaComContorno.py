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
    IEN = Tri.triangles
    print(X,Y,IEN)
    return X,Y,IEN

def plotMalhaTri(X,Y,IEN):
    name = 'MalhaTriangularComContorno'
    Path(os.path.join('results', name)).mkdir(parents=True, exist_ok=True)
    borda = verificaBorda(X, Y, IEN)
    print(borda)
    ne = IEN.shape[0]
    npoints = X.shape[0]
    crit_list = []
    A_sum = 0
    for e in range(0,ne):
        [v1,v2,v3] = IEN[e]
        x1 = X[v1]; x2 = X[v2]; x3 = X[v3]
        y1 = Y[v1]; y2 = Y[v2]; y3 = Y[v3]
        a = np.sqrt((x3-x2)**2+(y3-y2)**2)
        b = np.sqrt((x3-x1)**2+(y3-y1)**2)
        c = np.sqrt((x2-x1)**2+(y2-y1)**2)
        s = (a + b + c)/2
        A = np.sqrt(s*(s-a)*(s-b)*(s-c)) # Formula de Heron
        A_sum += A
        ha = 2*A/a
        hb = 2*A/b
        hc = 2*A/c
        crit = ((ha*hb*hc)/(a*b*c))/((np.sqrt(3)/2)**3)
        crit_list.append(crit)
        
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 10)
    xy = np.stack((X, Y), axis=-1)
    verts = xy[borda]
    ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts,edgecolors=('black',),
                                               facecolor = 'none',
                                                 linewidths=(5))
    ax.add_collection(pc)
    ax.set_title("Malha de Triângulos Com Contorno")
    plt.plot(X,Y,'ko')
    ax = plt.triplot(X,Y,IEN,color='k',linewidth=0.5)
    ax = plt.plot(X,Y,'ko',markersize=2)
    for i in range(0,npoints):
       plt.text(X[i]+0.01,Y[i]+0.01,str(i),size = 'x-small',color='k')
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
    plt.savefig(os.path.join('results', name, 'Malha_de_Triangulos_com_Contorno.png'))

def verificaBorda(X,Y,IEN):
    IENBounds = []
    arestas_compartilhadas = []
    arestas_isoladas = []
    for i in range(0,len(IEN)):
        [v1,v2,v3] = IEN[i]
        a1 = [v1,v2]
        a2 = [v2,v3]
        a3 = [v3,v1]
        for j in range(0,len(IEN)):
            if(i != j):
                [w1,w2,w3] = IEN[j]
                a11 = [w1,w2]
                a22 = [w2,w3]
                a33 = [w3,w1]
                a11_ = [w2,w1]
                a22_ = [w3,w2]
                a33_ = [w1,w3]
                a = [a11,a22,a33,a11_,a22_,a33_]
                if(a1 in a):
                    arestas_compartilhadas.append(a1)
                if(a2 in a):
                    arestas_compartilhadas.append(a2)
                if(a3 in a):
                    arestas_compartilhadas.append(a3)
        if(a1 not in arestas_compartilhadas):
            arestas_isoladas.append(a1)
        if(a2 not in arestas_compartilhadas):
            arestas_isoladas.append(a2)
        if(a3 not in arestas_compartilhadas):
            arestas_isoladas.append(a3)
    IENBounds = organizadorIEN(arestas_isoladas)
    return IENBounds

def organizadorIEN(lista):
    ## Percorre a lista de arestas isoladas, organizando pelos vertices
    print(lista)
    for i in range(0,len(lista)-1):
        for j in range(0,len(lista)):
            if(lista[i][1] == lista[j][0]):
                lista_org = lista[i+1]
                lista[i+1] = lista[j]
                lista[j] = lista_org
    ## Separa os contornos diferentes (ex.: externo e interno(buracos))
    split = []
    n = 0
    for m in range(0,len(lista)):
        if(lista[m-1][1]!=lista[m][0]):
            split.append(lista[i:m])
            i = m
    split.append(lista[n:])
    print(split)
    ## Armazena numa IEN os diferentes contornos, 
    ## identificando tanto os externos quanto internos
    IENs = []
    for b in range(0,len(split)):
        contorno = split[b] 
        print(contorno)
        IEN = np.zeros((len(contorno),1),dtype=int)
        IEN[0] = contorno[0][0]
        for k in range(0,len(contorno)-1):
            IEN[k+1] = contorno[k][1]
        IEN = IEN.reshape((1,len(contorno)))
        IENs.append(IEN)
    return IENs
        

if __name__ == '__main__':
    X,Y,Tri = getXYRandom()
    plotMalhaTri(X,Y,Tri)