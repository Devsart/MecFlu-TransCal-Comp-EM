# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:38:53 2021

@author: matheus.sartor
"""

"""
Equacionamento: Problema térmico transiente 2D
dT/dt = alpha*Nabla**2(T) + Q / rho*cv

MdTi/dt = alpha*K*Ti + M*Qi/(rho*cv)

discretização temporal:

M ( Tfuturo - Tpresente )/ dt
colocar num loop

Sempre construir as matrizes M e K do elemento

Elementos finitos:
    Joga todo mundo pra esquerda, integra e iguala a zero
    T(t=0) : condição inicial
    Q = 100
    
    Integral(w*(dT/dt - alpha*Nabla2(T) - Q/ rhocv )dA) = 0
    discretiza:
        contorno:
            Integral ( w*alhpa*nablaT dTau )
        miolo:
            Integral ( nablaW * alpha * nablaT dA )
        Tempo:
            Integral ( w * dT/dt * dA )
        condução:
            Integral ( w * Q / rho*cv * dA )
    
    Eq:
        Tempo - contorno + miolo - condução = 0

T= soma Ni(x,y)*ti(x,y,t)
Q = soma Nk(x,y)*qk
w = soma Nj(x,y)*wj

substituindo para Galerkin Ni=Nj=Nk:
    
    1) soma(i,j)  wj * dti/dt * Integral (Nj * Ni) dA
    2) 0 -> no contorno w é 0 -> Drichlet
    3) soma(i,j) alpha*ti*wj * Integral (NablaNi*NablaNj) dA
    4) soma(j,k) wj*qk / (rho*cv) * Integral (Nj*Nk) dA
    
    M => Integral (N * N) dA
    K => Integral (NablaN * NablaN) dA
    
    
Forma matricial:
    M*dti/dt + alpha * K * ti - M * qk / (rho*cv) = 0
    
Discretização temporal:
    M (t_futuro - t_presente)/dt + alpha*K*ti - M*qk/rho*cv = 0
    
    ti explicito ou implicito?
    
    Pode-se usar theta


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


name = 'calor_transiente_mef_2d_triangular'
Path(os.path.join(name,'resultados_implicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_explicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_crank_nicholson')).mkdir(parents=True, exist_ok=True)

def montaKM(X,Y,IEN):
    npoints = X.shape[0]
    alpha = np.ones(len(IEN), dtype='float')
    K = np.zeros( (npoints,npoints),dtype='float' )
    M = np.zeros( (npoints,npoints),dtype='float' )
    for elem in range(1,len(IEN)):
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
        ck = xi - xk
        cj = xj - xi
        if(yi < 0.5 and yk < 0.5 and yj < 0.5):
            alpha[elem] = 0.1
        A = np.absolute((1/2)*np.linalg.det([[1,xi,yi],[1,xj,yj],[1,xk,yk]]))
        ke_x = (1/(4*A))*np.array([[bi**2,bi*bj,bi*bk],[bj*bi,bj**2,bj*bk],[bk*bi,bk*bj,bk**2]])
        ke_y = (1/(4*A))*np.array([[ci**2,ci*cj,ci*ck],[cj*ci,cj**2,cj*ck],[ck*ci,ck*cj,ck**2]])
        ke = ke_x + ke_y         
        me = (A/12)*np.array([[2,1,1],[1,2,1],[1,1,2]])
        for i_loc in range(0,3):
            i_glb = IEN[elem,i_loc]
            for j_loc in range(0,3):
                j_glb = IEN[elem,j_loc]
                
                K[i_glb,j_glb] += alpha[elem]*ke[i_loc,j_loc]
                M[i_glb,j_glb] += me[i_loc,j_loc]
    return K,M,npoints

def contornoPlaca(nx, ny):
    ccL = []
    ccR = []
    ccTop = []
    ccBot = []
    cc = []
    inner = []
    for i in range(1,ny-1):
        ccL.append(nx*i)    
        ccR.append(nx*(i+1) - 1)
    for j in range(nx):
        ccBot.append(j)
        ccTop.append(j + (ny-1)*nx)
    cc=[*ccBot,*ccL,*ccR,*ccTop]
    space = np.linspace(0,(nx*ny-1),36)
    for k in space:
        if k not in cc:
            inner.append(int(k))
    print(f"ccL = {ccL}\nccR = {ccR}\nccTop = {ccTop}\nccBot = {ccBot}\ncc = {cc}\ninner = {inner}")
    return ccL,ccR,ccTop,ccBot,cc, inner

def Drichlet2D(npoints,K,M,ccL,ccR,ccTop,ccBot,q=1.,rho=1.,cv=1.):
    Q = (q/(rho*cv))*np.ones((npoints),dtype='float')
    A = K.copy()
    b = M@Q
    for i in ccL:
        A[i,:] = 0.0
        A[i,i] = 1.0
        b[i] = 0.
    for i in ccR:
        A[i,:] = 0.0
        A[i,i] = 1.0
        b[i] = 0.
    for i in ccTop:
        A[i,:] = 0.0
        A[i,i] = 1.0
        b[i] = 0.
    for i in ccBot:
        A[i,:] = 0.0
        A[i,i] = 1.0
        b[i] = 0.
    return A,Q,b
    
def solveWithTheta(theta = 1,dt=0.01,lim_e=1e-5):
    """
    Equation
    --------
    (M/dt)*T+theta*alpha*K*T = (M/dt)*T0 + Q - (1-theta)*alpha*K*T0
    

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
    alpha= 1.
    rho=1.
    cv=1.
    q=1000.
    e=1.
    Ti = 27.
    Tc = Ti

    X,Y,IEN = pm.getXY(nx,ny,1,1,"triangle")
    K,M,npoints = montaKM(X,Y,IEN)
    ccL,ccR,ccTop,ccBot,cc,inner = contornoPlaca(nx, ny)
    A,Q,beta = Drichlet2D(npoints,K,M,ccL,ccR,ccTop,ccBot)
    
    T = Ti*np.ones((npoints),dtype='float')
    for c in cc:
        T[c] = Tc
    
    Q = 0.1*np.ones((npoints),dtype='float')
    for i in range(0,npoints):
        if X[i] <=0.4 and X[i] >= 0.2 and Y[i] <= 0.4 and Y[i]>=0.2: 
            Q[i] = q/(rho*cv)
        if X[i] <=0.8 and X[i] >= 0.6 and Y[i] <= 0.4 and Y[i]>=0.2: 
            Q[i] = (q/(rho*cv))
        if X[i] <=0.4 and X[i] >= 0.2 and Y[i] <= 0.8 and Y[i]>=0.6: 
            Q[i] = (q/(rho*cv))
        if X[i]<=0.8 and X[i] >= 0.6 and Y[i] <= 0.8 and Y[i]>=0.6: 
            Q[i] = (q/(rho*cv))
    index = 0
    while(e > lim_e):
        Tpast = T
        A = (M/dt + theta*K)
        b = ((M/dt)@T + M@Q - ((1-theta)*K)@T)
        for i in ccL:
            A[i,:] = 0.0
            A[i,i] = 1.0
            b[i] = Tc
        for i in ccR:
            A[i,:] = 0.0
            A[i,i] = 1.0
            b[i] = Tc
        for i in ccTop:
            A[i,:] = 0.0
            A[i,i] = 1.0
            b[i] = Tc
        for i in ccBot:
            A[i,:] = 0.0
            A[i,i] = 1.0
            b[i] = Tc
        Ainv = np.linalg.inv(A)
        T = Ainv@b
        Z = T.reshape(ny,nx)
        xy = np.stack((X, Y), axis=-1)
        verts = xy[IEN]
        fig, ax = plt.subplots()
        fig.set_size_inches(8, 8)
        ax=plt.gca()
        pc = matplotlib.collections.PolyCollection(verts,edgecolors=('black',),
                                                      facecolor='None',
                                                      linewidths=(0.7,))
        ax.add_collection(pc)
        plt.plot(X,Y,'k.')
        plt.title('Distrubuição de calor no regime transiente 2D')
        plt.xlabel('comprimento da placa no eixo X [m]')
        plt.ylabel('comprimento da placa no eixo Y [m]')
        surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                          cmap=matplotlib.cm.jet, extent=(X.min(),
                          X.max(), Y.min(), Y.max()))
        plt.clim(0, 100)
        cbar = plt.colorbar(surf,shrink=1.0, aspect=20)
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
        #plt.savefig(os.path.join(name, subname, f'{index:04}.png'))
        plt.show()
        index += 1
        ## mean squared error MSE
        e = np.sum((T - Tpast)**2)/T.shape[0]
    # fp_in = f"./{name}/{subname}/*.png"
    # fp_out = f"./{name}/{subname}/simulacao.gif"
    
    # img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    # img.save(fp=fp_out, format='GIF', append_images=imgs,
    #          save_all=True, duration=dt*1000, loop=0)
solveWithTheta(1)
    