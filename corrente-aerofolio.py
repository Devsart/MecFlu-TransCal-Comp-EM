# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:38:53 2021

@author: matheus.sartor
"""

import plotMalhas as pm
import MalhaComContorno as contorno
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as mtri
import numpy as np
import os
from pathlib import Path
import glob
from PIL import Image
import meshio



name = 'funcao_corrente_mef_2d_triangular'
Path(os.path.join(name,'resultados')).mkdir(parents=True, exist_ok=True)

def montaKM(X,Y,IEN,regions):
    npoints = len(X)
    K = np.zeros( (npoints,npoints),dtype='float' )
    M = np.zeros( (npoints,npoints),dtype='float' )
    Gx = np.zeros( (npoints,npoints),dtype='float' )
    Gy = np.zeros( (npoints,npoints),dtype='float' )
    for elem in range(0,len(IEN)):
        [v1,v2,v3] = IEN[elem]
        xi,yi = X[v1],Y[v1]
        xj,yj = X[v2],Y[v2]
        xk,yk = X[v3],Y[v3]
        ai = xj*yk - xk*yj
        aj = xk*yi - xi*yk
        ak = xi*yj - xj*yj
        bi = yj - yk
        bj = yk - yi
        bk = yi - yj
        ci = xk - xj
        cj = xi - xk
        ck = xj - xi
        Area = (1.0/2.0)*np.linalg.det([[1.0,xi,yi],[1.0,xj,yj],[1.0,xk,yk]])
        B = (1.0/(2.0*Area)) * np.array([ [bi, bj, bk],[ci, cj, ck] ])
        # matriz da matriz B transposta
        BT = B.transpose()       
        # matriz de rigidez do elemento
        ke = Area*np.dot(BT,B)
        ke_x = (1.0/(4.0*Area))*np.array([[bi*bi,bi*bj,bi*bk],[bj*bi,bj*bj,bj*bk],[bk*bi,bk*bj,bk*bk]])
        ke_y = (1.0/(4.0*Area))*np.array([[ci*cj,ci*cj,ci*ck],[cj*ci,cj*cj,cj*ck],[ck*ci,ck*cj,ck*ck]])
        #ke = ke_x + ke_y
        ge_x = (1.0/6.0)*np.array([[bi,bj,bk],[bi,bj,bk],[bi,bj,bk]],dtype='float64')
        ge_y = (1.0/6.0)*np.array([[ci,cj,ck],[ci,cj,ck],[ci,cj,ck]],dtype='float64')
        ge = ge_x + ge_y
        me = (Area/12.0)*np.array([[2.0,1.0,1.0],[1.0,2.0,1.0],[1.0,1.0,2.0]])
        for i_loc in range(0,3):
            i_glb = IEN[elem,i_loc]
            for j_loc in range(0,3):
                j_glb = IEN[elem,j_loc]
                
                K[i_glb,j_glb] += ke[i_loc,j_loc]
                M[i_glb,j_glb] += me[i_loc,j_loc]
                Gx[i_glb,j_glb] += ge_x[i_loc,j_loc]
                Gy[i_glb,j_glb] += ge_y[i_loc,j_loc]
    return K,M,Gx,Gy

    
def solveCorrente():
    """
    Equation
    --------
    (M/dt+v*K+theta*V*G)@w = ((M/dt) - (1-theta)*v*G)@w
    

    Parameters
    ----------
    theta : float, optional
        theta = 0 - Explicit solver.
        theta = 1 - Implicit solver
        theta = 1/2 - Crank-Nicholson solver
        The default value is 0.

    Returns
    -------
    Matrix for transient Temperature solution.

    """
    nx = 10
    ny = 10
    c2 = 1
    c1 = 0
    v = .1
    e=1.
    msh = meshio.read('aerofolio.msh')
    X = msh.points[:,0]
    Y = msh.points[:,1]
    IEN = msh.cells['triangle']
    IENbound = msh.cells['line']
    IENboundTypeElem = list(msh.cell_data['line']['gmsh:physical'] - 1)
    boundNames = list(msh.field_data.keys())
    IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]
    npoints = len(X)
    ne = IEN.shape[0]
    
    # cria lista de nos do contorno
    cc = np.unique(IENbound.reshape(IENbound.size))
    ccName = [[] for i in range( npoints )]
    for elem in range(0,len(IENbound)):
    	ccName[ IENbound[elem][0] ] = IENboundElem[elem]
    	ccName[ IENbound[elem][1] ] = IENboundElem[elem]
    regions = msh.cell_data['triangle']['gmsh:geometrical']

    vx = np.zeros(npoints,dtype='float64')
    vy = np.zeros((npoints),dtype='float64')
    K,M,Gx,Gy = montaKM(X,Y,IEN,regions)
    index = 0
    Minv = np.linalg.inv(M)
    A = K.copy()
    b = np.zeros( (npoints),dtype='float' )
    bc = 0*np.ones( (npoints),dtype='float' )
    for i in cc:
        A[i,:] = 0.0
        A[i,i] = 1.0
        print(cc)
        if ccName[i] == 'top':
            bc[i] = 1.0
            b[i] = bc[i]
        if ccName[i] == 'bottom':
            b[i] = 0.0
        if ccName[i] == 'left':
            b[i] = Y[i]
        if ccName[i] == 'right':
            b[i] = Y[i]
        if ccName[i] == 'hole':
            b[i] = (max(Y)-min(Y))/2
    ## Solucao vorticidade
    Ainv = np.linalg.inv(A)
    # b -= bc
    psi = Ainv@b
    ## Campo de velocidades
    b_3 = np.dot(Gy,psi)
    M3 = M.copy()
    vx = np.linalg.solve(M3,b_3)
    
    b_4 = np.dot(Gx,psi)
    M4 = M.copy()
    vy = -np.linalg.solve(M4,b_4)
    ## resolvendo contornos
    Z1 = psi
    Z2 = vx
    print(Z2)
    Z3 = vy
    print(Z3)
    xy = np.stack((X, Y), axis=-1)
    verts = xy[IEN]
    verts2 = xy[IENbound]
    fig, ax = plt.subplots()
    triang = mtri.Triangulation(X,Y,IEN)
    surf = ax.tricontourf(triang,Z3,50,cmap=matplotlib.cm.jet)
    # surf2 = ax[0,1].tricontourf(X,Y,Z2,50,cmap=matplotlib.cm.jet, extent=(X.min(),
    #                    X.max(), Y.min(), Y.max()),vmin=0,vmax=1)
    # surf3 = ax[1,0].tricontourf(X,Y,Z3,50,cmap=matplotlib.cm.jet, extent=(X.min(),
    #                    X.max(), Y.min(), Y.max()),vmin=0,vmax=1)
    fig.set_size_inches(8, 4)
    ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts,edgecolors=('lightgray',),
                                                  facecolor='None',
                                                  linewidths=(0.2,))
    ax.add_collection(pc)
    pc2 = matplotlib.collections.PolyCollection(verts2,edgecolors=('black',),
                                                  facecolor='None',
                                                  linewidths=(0.7,))
    ax.add_collection(pc2)
    #ax[1,1].plot(X,Y,'k.')
    cbar = plt.colorbar(surf,shrink=1.0, aspect=20)
    cbar.set_label('Velocidade $v_{y}$[m/s]')
    plt.title(r'Solução da Equação: $v_{y} = M^{-1}(-G_{x}\psi)$')
    labx = np.linspace(X.min(),X.max(),nx)
    laby = np.linspace(Y.min(),Y.max(),ny)
    plt.xticks(labx)
    plt.yticks(laby)
    subname = 'resultados'
    plt.savefig(os.path.join(name, subname, f'velocidade_y_aerofolio.png'))
    plt.show()

solveCorrente()
# msh = meshio.read('aerofolio.msh')
# X = msh.points[:,0]
# Y = msh.points[:,1]
# IEN = msh.cells['triangle']
# IENbound = msh.cells['line']
# IENboundTypeElem = list(msh.cell_data['line']['gmsh:physical'] - 1)
# boundNames = list(msh.field_data.keys())
# IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]
# npoints = len(X)
# ne = IEN.shape[0]

# # cria lista de nos do contorno
# cc = np.unique(IENbound.reshape(IENbound.size))
# ccName = [[] for i in range( npoints )]
# for elem in range(0,len(IENbound)):
# 	ccName[ IENbound[elem][0] ] = IENboundElem[elem]
# 	ccName[ IENbound[elem][1] ] = IENboundElem[elem]
# regions = msh.cell_data['triangle']['gmsh:geometrical']
# contorno = msh.cells['line']
# parede_top = [1,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,2]
# parede_bot = [3,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,0]
# entrada = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
# saida = [67,68,69,70,71,72,73,74,75,76,77,78,79,80,81]

# print(ccName)  