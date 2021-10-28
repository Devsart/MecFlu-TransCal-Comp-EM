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
import meshio



name = 'calor_transiente_mef_2d_triangular_final'
Path(os.path.join(name,'resultados_implicito_placa_dell_v6')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_explicito_placa_dell')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_crank_nicholson_placa_dell')).mkdir(parents=True, exist_ok=True)

def montaKM(X,Y,IEN,regions):
    npoints = X.shape[0]
    q = 65
    cooler = 1
    if(cooler == 1):
            q -= 33.4
    cpSi = .867  ## 
    rhoSi = 2.659
    kappa_Si = 1.61
    alpha_SiAl = kappa_Si/(rhoSi*cpSi)#(7.172*1e-5)
    cpPet = .712 ## https://www.plastmetal.com.br/tabelas/820116434bbd/tabela_de_propriedades_polietileno.pdf SINTMID
    rhoPet = 1.790
    kappa_Pet = 0.22
    alpha_ladeVidro = kappa_Pet/(cpPet*rhoPet)#(22.6*1e-7) # https://www.braskem.com.br/Portal/Principal/Arquivos/html/boletm_tecnico/Tabela_de_Propriedades_de_Referencia_dos_Compostos_de_PVC.pdf
    cpCu = .383
    rhoCu = 8.960
    kappa_Cu = 3.86
    alpha_Cu = kappa_Cu/(rhoCu*cpCu)#(11.234*1e-5)
    cpAl = .921
    rhoAl = 2.700
    kappa_Al = 2.04
    alpha_Al = kappa_Al/(rhoAl*cpAl)#(8.418*1e-5)
    alpha = np.ones(len(IEN), dtype='float')
    kappa = np.ones(len(IEN), dtype='float')
    rho = np.ones(len(IEN), dtype='float')
    cv = np.ones(len(IEN), dtype='float')
    Q = 0*np.ones((npoints),dtype='float')
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
        if(regions[elem] == 1): ##Base da Placa
            kappa[elem] = kappa_Pet
            cv[elem] = cpPet
            rho[elem] = rhoPet
            alpha[elem] = alpha_ladeVidro
        if(regions[elem] == 2): ## Alimentador
            kappa[elem] = kappa_Al
            cv[elem] = cpAl
            rho[elem] = rhoAl
            alpha[elem] = alpha_Al
            Q[v1] = q*0.1/(rho[elem]*cv[elem])
            Q[v2] = q*0.1/(rho[elem]*cv[elem])
            Q[v3] = q*0.1/(rho[elem]*cv[elem])
        if(regions[elem] in [3,4]): ## Processador AMD4
            kappa[elem] = kappa_Si
            cv[elem] = cpSi
            rho[elem] = rhoSi
            alpha[elem] = alpha_SiAl
            Q[v1] = q*0.15/(rho[elem]*cv[elem])
            Q[v2] = q*0.15/(rho[elem]*cv[elem])
            Q[v3] = q*0.15/(rho[elem]*cv[elem])
        if(regions[elem] == 5): ## SSD M2
            kappa[elem] = kappa_Cu
            cv[elem] = cpCu
            rho[elem] = rhoCu
            alpha[elem] = alpha_Cu
            Q[v1] = q*0.1/(rho[elem]*cv[elem])
            Q[v2] = q*0.1/(rho[elem]*cv[elem])
            Q[v3] = q*0.1/(rho[elem]*cv[elem])
        if(regions[elem] in [6,7]): ## Memoria RAM DDR4
            kappa[elem] = kappa_Cu
            cv[elem] = cpCu
            rho[elem] = rhoCu
            alpha[elem] = alpha_Cu
            Q[v1] = q*0.1/(rho[elem]*cv[elem])
            Q[v2] = q*0.1/(rho[elem]*cv[elem])
            Q[v3] = q*0.1/(rho[elem]*cv[elem])
        if(regions[elem] == 8): ## Placa de Video SATA 6Gbs
            kappa[elem] = kappa_Si
            cv[elem] = cpSi
            rho[elem] = rhoSi
            alpha[elem] = alpha_SiAl
            Q[v1] = q*0.15/(rho[elem]*cv[elem])
            Q[v2] = q*0.15/(rho[elem]*cv[elem])
            Q[v3] = q*0.15/(rho[elem]*cv[elem])
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
    return K,M,Q,npoints

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
    
def solveWithTheta(theta = 1,dt=0.5,lim_e=1e-5):
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
    # alpha= 1.
    # rho=1.
    # cv=1.
    # q=1000.
    e=1.
    Ti = 27.
    Tc = Ti

    msh = meshio.read('placa-mae-dell-inspiron-5547.msh')
    X = msh.points[:,0]
    Y = msh.points[:,1]
    IEN = msh.cells['triangle']
    npoints = len(X) 
    ne = IEN.shape[0]
    regions = msh.cell_data['triangle']['gmsh:geometrical']
    contorno_mae_v2 = [0,36,37,38,3,39,40,41,42,43,44,45,2,46,47,48,1,49,50,51,52,53,54,55]
    contorno_dell = [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,2,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,3,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,3,213,214,30,219,220,0,244,28,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,1]


    K,M,Q,npoints = montaKM(X,Y,IEN,regions)
    cc = IENbound = msh.cells['line']
    # A,Q,beta = Drichlet2D(npoints,K,M,ccL,ccR,ccTop,ccBot)
    
    T = Ti*np.ones((npoints),dtype='float')
    for c in contorno_dell:
        T[c] = Tc
    
    # for i in range(0,npoints):
    #     if X[i] <=0.4 and X[i] >= 0.2 and Y[i] <= 0.4 and Y[i]>=0.2: 
    #         Q[i] = q/(rho*cv)
    #     if X[i] <=0.8 and X[i] >= 0.6 and Y[i] <= 0.4 and Y[i]>=0.2: 
    #         Q[i] = (q/(rho*cv))
    #     if X[i] <=0.4 and X[i] >= 0.2 and Y[i] <= 0.8 and Y[i]>=0.6: 
    #         Q[i] = (q/(rho*cv))
    #     if X[i]<=0.8 and X[i] >= 0.6 and Y[i] <= 0.8 and Y[i]>=0.6: 
    #         Q[i] = (q/(rho*cv))
    index = 0
    t = 0
    while(index < 60):
        Tpast = T
        # A = (M/dt + theta*K)
        # b = ((M/dt)@T + M@Q - ((1-theta)*K)@T)
        A = (K)
        b = (M@Q)
        for i in contorno_dell:
            A[i,:] = 0.0
            A[i,i] = 1.0
            b[i] = Tc
        # for i in ccR:
        #     A[i,:] = 0.0
        #     A[i,i] = 1.0
        #     b[i] = Tc
        # for i in ccTop:
        #     A[i,:] = 0.0
        #     A[i,i] = 1.0
        #     b[i] = Tc
        # for i in ccBot:
        #     A[i,:] = 0.0
        #     A[i,i] = 1.0
        #     b[i] = Tc
        Ainv = np.linalg.inv(A)
        T = Ainv@b
        Z = T
        xy = np.stack((X, Y), axis=-1)
        verts = xy[IEN]
        verts2 = xy[cc]
        indice, = np.where(T == T.max())
        XTmax = X[indice]
        YTmax = Y[indice]
        Tmax = T[indice]
        t += dt
        fig, ax = plt.subplots()
        surf = ax.tricontourf(X,Y,Z,200,cmap=matplotlib.cm.jet, extent=(X.min(),
                           X.max(), Y.min(), Y.max()),vmin=0,vmax=350)
        fig.set_size_inches(10, 10)
        ax=plt.gca()
        pc = matplotlib.collections.PolyCollection(verts,edgecolors=('lightgray',),
                                                      facecolor='None',
                                                      linewidths=(0.7,))
        ax.add_collection(pc)
        pc2 = matplotlib.collections.PolyCollection(verts2,edgecolors=('black',),
                                                      facecolor='None',
                                                      linewidths=(0.7,))
        ax.add_collection(pc2)
        #plt.plot(X,Y,'w.')
        plt.title(r'Solution of $\frac{M}{dt} T^{n+1} + \theta \alpha K T^{n+1} = \frac{M}{dt} T^{n} + Q^{n} + ( 1 - \theta ) \alpha K T^{n} $'+"\n$T_{max}$"+f" = {Tmax[0]}, "+"$X_{Tmax}$ = "+f"{XTmax[0]}, "+"$Y_{Tmax}$ = "+f"{YTmax[0]}"+f"\n # de nós = {npoints}; # de elementos = {ne}\ndt = {dt};  t = {t}")
        plt.xlabel('comprimento da placa no eixo X [cm]')
        plt.ylabel('comprimento da placa no eixo Y [cm]')
        #surf = plt.imshow(X,Y,Z, interpolation='quadric', origin='lower',
        #                  cmap=matplotlib.cm.jet, extent=(X.min(),
        #                  X.max(), Y.min(), Y.max()))
        
        # plt.clim(0, 300)
        cbar = plt.colorbar(surf,shrink=1.0, aspect=20)
        cbar.set_label('Temperatura [°C]')
        labx = np.linspace(X.min(),X.max(),nx)
        laby = np.linspace(Y.min(),Y.max(),ny)
        plt.xticks(labx)
        plt.yticks(laby)
        subname = '???'
        if(theta == 0):
            subname = 'resultados_explicito_placa_dell'
        if(theta == 1):
            subname = 'resultados_implicito_placa_dell_v6'
        if(theta == 1/2):
            subname = 'resultados_crank_nicholson_placa_dell'
        plt.savefig(os.path.join(name, subname, f'{index:04}.png'))
        plt.show()
        index += 1
        ## mean squared error MSE
        e = np.sum((T - Tpast)**2)/T.shape[0]
    fp_in = f"./{name}/{subname}/*.png"
    fp_out = f"./{name}/{subname}/simulacao.gif"
    
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(fp=fp_out, format='GIF', append_images=imgs,
              save_all=True, duration=dt*1000, loop=0)
solveWithTheta(1,1)
# msh = meshio.read('placa-mae-dell-inspiron-5547.msh')
# X = msh.points[:,0]
# Y = msh.points[:,1]
# IEN = msh.cells['triangle']
# npoints = len(X) 
# ne = IEN.shape[0]
# regions = msh.cell_data['triangle']['gmsh:geometrical']
# contorno = [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,2,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,3,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,3,213,214,30,219,220,0,244,28,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,1]

# print(X)  