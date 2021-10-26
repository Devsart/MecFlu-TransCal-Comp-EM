# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:38:53 2021

@author: matheus.sartor
"""

import plotMalhas as pm
import MalhaComContorno as contorno
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import os
from pathlib import Path
import glob
from PIL import Image
import meshio



name = 'funcao_corrente_mef_2d_triangular'
Path(os.path.join(name,'resultados_implicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_explicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_crank_nicholson')).mkdir(parents=True, exist_ok=True)

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

    
def solveCorrente(theta = 1,dt=0.5,lim_e=1e-5):
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
    msh = meshio.read('duto.msh')
    X = msh.points[:,0]
    Y = msh.points[:,1]
    IEN = msh.cells['triangle']
    npoints = len(X) 
    ne = IEN.shape[0]
    regions = msh.cell_data['triangle']['gmsh:geometrical']
    ##contorno_mae_v2 = [0,36,37,38,3,39,40,41,42,43,44,45,2,46,47,48,1,49,50,51,52,53,54,55]
    ##contorno_dell = [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,2,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,3,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,3,213,214,30,219,220,0,244,28,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,1]
    parede_top = [1,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,2]
    parede_bot = [3,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0]
    entrada = [4,5,6,7,8,9,10]
    saida = [34,35,36,37,38,39,40]
    contorno = [4,5,6,7,8,9,10,1,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,2,34,35,36,37,38,39,40,3,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,0]
    vx = np.zeros(npoints,dtype='float64')
    vy = np.zeros((npoints),dtype='float64')
    # psi = (((c2-c1)/(Y.max()-Y.min()))*Y + c1 - ((c2-c1)/(Y.max()-Y.min()))*Y.min())*np.ones((npoints),dtype='float64')
    # for i in parede_top:
    #     psi[i] = c2
    # for j in parede_bot:
    #     psi[j] = c1
    # for k in entrada:
    #     vx[k] = 1.0
    #     psi[k] = ((c2-c1)/(Y.max()-Y.min()))*Y[k] + c1 - ((c2-c1)/(Y.max()-Y.min()))*Y.min()
    # for l in saida:
    #     psi[l] = ((c2-c1)/(Y.max()-Y.min()))*Y[l] + c1 - ((c2-c1)/(Y.max()-Y.min()))*Y.min()
    K,M,Gx,Gy = montaKM(X,Y,IEN,regions)
    cc = IENbound = msh.cells['line']
    # A,Q,beta = Drichlet2D(npoints,K,M,ccL,ccR,ccTop,ccBot)
    index = 0
    Minv = np.linalg.inv(M)
    A = K.copy()
    b = np.zeros( (npoints),dtype='float' )
    bc = 0*np.ones( (npoints),dtype='float' )
    for i in parede_top:
        bc[i] = 1.0
    for j in parede_bot:
        bc[j] = 0.0
    for k in entrada:
        vx[k] = 1.0
        vy[k] = 0.0
        bc[k] = Y[k]
    for l in saida:
        bc[l] = Y[l]
    for j in contorno:
        # A[:,j] = 0
        A[j,:] = 0.0
        A[j,j] = 1.0
        b[j] = bc[j]
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
    verts2 = xy[cc]
    fig, ax = plt.subplots(2,2)
    surf = ax[0,0].tricontourf(X,Y,Z1,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                       X.max(), Y.min(), Y.max()),vmin=0,vmax=1)
    surf2 = ax[0,1].tricontourf(X,Y,Z2,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                       X.max(), Y.min(), Y.max()),vmin=0,vmax=1)
    surf3 = ax[1,0].tricontourf(X,Y,Z3,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                       X.max(), Y.min(), Y.max()),vmin=0,vmax=1)
    fig.set_size_inches(9, 4)
    ax=plt.gca()
    pc = matplotlib.collections.PolyCollection(verts,edgecolors=('lightgray',),
                                                  facecolor='None',
                                                  linewidths=(0.7,))
    ax.add_collection(pc)
    pc2 = matplotlib.collections.PolyCollection(verts2,edgecolors=('black',),
                                                  facecolor='None',
                                                  linewidths=(0.7,))
    ax.add_collection(pc2)
    #ax[1,1].plot(X,Y,'k.')
    cbar = ax[0,0].colorbar(surf,shrink=1.0, aspect=20)
    cbar.set_label('Temperatura [°C]')
    labx = np.linspace(X.min(),X.max(),nx)
    laby = np.linspace(Y.min(),Y.max(),ny)
    plt.xticks(labx)
    plt.yticks(laby)
    subname = '???'
    if(theta == 0):
        subname = 'resultados_explicito'
    if(theta == 1):
        subname = 'resultados_implicito'
    if(theta == 1/2):
        subname = 'resultados_crank_nicholson'
    plt.savefig(os.path.join(name, subname, f'{index:04}.png'))
    plt.show()
    index += 1

#solveCorrente(1,0.1)
msh = meshio.read('duto-furo.msh')
X = msh.points[:,0]
Y = msh.points[:,1]
IEN = msh.cells['triangle']
npoints = len(X) 
ne = IEN.shape[0]
regions = msh.cell_data['triangle']['gmsh:geometrical']
contorno = msh.cells['line']

print(msh.cells['line'])  