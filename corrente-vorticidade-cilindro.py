# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 08:38:53 2021

@author: matheus.sartor
"""

import plotMalhas as pm
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.tri as mtri
import numpy as np
import os
from pathlib import Path
import glob
from PIL import Image
import meshio
import copy



name = 'funcao_corrente_vorticidade_mef_2d_triangular'
Path(os.path.join(name,'resultados_cilindro_implicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_explicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_crank_nicholson')).mkdir(parents=True, exist_ok=True)

def getContorno(IENBound):
    contorno = []
    entrada = []
    top = []
    bot = []
    saida = []
    for i in range(len(IENBound)):
        contorno.append(IENBound[i][1])
    for index,j in enumerate(contorno):
        if (j == 3):
            entrada = contorno[:index]
            i_ent = index
        if (j == 2):
            top = contorno[i_ent:index]
            i_top = index
        if (j == 1):
            saida = contorno[i_top:index]
            i_ext = index
        if (j == 0):
            bot = contorno[i_ext:]
    return contorno,entrada,top,saida,bot


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
        ke_y = (1.0/(4.0*Area))*np.array([[ci*ci,ci*cj,ci*ck],[cj*ci,cj*cj,cj*ck],[ck*ci,ck*cj,ck*ck]])
        # ke = ke_x + ke_y
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
    
def solveWithTheta(theta = 1,dt=0.1,lim_e=1e-5):
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
    
    v = .001
    e=1.
    msh = meshio.read('duto-furo.msh')
    X = msh.points[:,0]
    Y = msh.points[:,1]
    IEN = msh.cells['triangle']
    cc = msh.cells['line']
    npoints = len(X) 
    ne = IEN.shape[0]
    regions = msh.cell_data['triangle']['gmsh:geometrical']
    ##contorno_mae_v2 = [0,36,37,38,3,39,40,41,42,43,44,45,2,46,47,48,1,49,50,51,52,53,54,55]
    ##contorno_dell = [31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,2,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,3,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,236,237,238,239,240,241,242,243,3,213,214,30,219,220,0,244,28,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,260,261,262,263,264,265,266,267,1]
    parede_top = [1,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,2]
    parede_bot = [3,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,127,128,0]
    entrada = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19]
    saida = [67,68,69,70,71,72,73,74,75,76,77,78,79,80,81]
    furo = [4]
    for i in range(129,208):
        furo.append(i)
    contorno = [*entrada,*parede_top,*saida,*parede_bot,*furo]
    Yfuro = Y[furo]
    centro = (Yfuro.max() + Yfuro.min())/2
    # contorno,entrada,parede_top,saida,parede_bot = getContorno(IENBound)
    w = np.zeros((npoints),dtype='float64')
    vx = np.zeros((npoints),dtype='float64')
    vy = np.zeros((npoints),dtype='float64')
    
    for i in parede_top:
        vx[i] = 0.0
        vy[i] = 0.0
    for i in parede_bot:
        vx[i] = 0.0
        vy[i] = 0.0
    for i in entrada:
        vx[i] = 1.0
        vy[i] = 0.0
    psi = np.zeros((npoints),dtype='float64')
    K,M,Gx,Gy = montaKM(X,Y,IEN,regions)
    # A,Q,beta = Drichlet2D(npoints,K,M,ccL,ccR,ccTop,ccBot)
    index = 0
    Minv = np.linalg.inv(M)
    gv = Gx@vy-Gy@vx
    w = Minv@gv
    while(index <= 999):
        #w = omegacc.copy()
        #mapear um n ́o e quais elementos ele pertence e
        #anotar quais n ́os s~ao esses[
        # R = []
        # for NODE in range(0,npoints):
        #     V=[]
        #     for i in range(0,ne):
        #         if IEN[i,0] == NODE :
        #             V.append(IEN[i,1])
        #             V.append(IEN[i,2])
        #         if IEN[i,1] == NODE :
        #             V.append(IEN[i,0])
        #             V.append(IEN[i,2])
        #         if IEN[i,2] == NODE :
        #             V.append(IEN[i,0])
        #             V.append(IEN[i,1])
        #     V = list(set(V))
        #     R.append(V)
        
        # p = np.zeros((npoints,1),dtype = 'float')
        # ui = np.zeros((npoints,1), dtype = 'float')
        # vi = np.zeros((npoints,1), dtype = 'float')
        
        # #calculo da velocidade ponderada pela distancia
        # #de cada n ́o vizinho
        # for i in range(0,len(R)):
        #     VIZINHOS = R[i]
        #     sum_xj = 0.0
        #     sum_yj = 0.0
        #     sumdist = 0.0
        #     for j in VIZINHOS:
        #         dist = np.sqrt(((X[j]-X[i])**2)+((Y[j]-Y[i])**2))
        #         sumdist += dist
        #         sum_xj += X[j]*dist
        #         sum_yj += Y[j]*dist
    
        #     ui[i] = 0.2*((sum_xj/sumdist)-X[i])/dt
        #     #valor da velocidade da malha
        #     vi[i] = 0.2*((sum_yj/sumdist)-Y[i])/dt
        #     #valor da velocidade da malha
        
        # unew = vx - ui
        # vnew = vy - vi
        
        I = np.eye(len(vx),dtype="float64")
        vxI = vx*I
        vyI = vy*I
        # vxnew = unew*I
        # vynew = vnew*I
        
        VGO = vxI@Gx + vyI@Gy
        
        A_w = (M.copy()/dt+v*K.copy()+theta*(VGO))
        b_w = ((M.copy()/dt) - (1-theta)*(VGO))@w
        for i in contorno:
            A_w[i,:] = 0.0
            A_w[i:i] = 1.0
            b_w[i] = w[i]
        ## Solucao vorticidade
        A_winv = np.linalg.pinv(A_w)
        w = A_winv@b_w
        ## Solucao funcao corrente
        A_psi = K.copy()
        b_psi = M@w
        for i in parede_top:
            A_psi[i,:] = 0.0
            A_psi[i,i] = 1.0
            b_psi[i] = 1.0
        for i in parede_bot:
            A_psi[i,:] = 0.0
            A_psi[i,i] = 1.0
            b_psi[i] = 0.0
        for i in entrada:
            A_psi[i,:] = 0.0
            A_psi[i,i] = 1.0
            b_psi[i] = Y[i]
        for i in furo:
            A_psi[i,:] = 0.0
            A_psi[i,i] = 1.0
            b_psi[i] = centro
        # for i in saida:
        #     b_psi[i] = 0.0
        # for l in saida:
        #     S[l,:] = 0.0
        #     S[l,l] = 1.0
        #     Mw[l] = 0.0 #Y[l]
        A_psiinv = np.linalg.inv(A_psi)
        psi = A_psiinv@b_psi
        ## Campo de velocidades
        Gypsi = Gy@psi
        vx = Minv@Gypsi
        Gxpsi = Gx@psi
        vy = -Minv@Gxpsi
        ## resolvendo contornos
        for i in parede_top:
            vx[i] = 0.0
            vy[i] = 0.0
        for i in parede_bot:
            vx[i] = 0.0
            vy[i] = 0.0
        for i in entrada:
            vx[i] = 1.0
            vy[i] = 0.0
        # for k in saida:
        #     vx[k] = 1.0
        #     vy[k] = 0.0
        Minv = np.linalg.inv(M)
        gv = Gx@vy-Gy@vx
        omegacc = Minv@gv
        for i in contorno:
            w[i] = omegacc[i]
        Z1 = w
        Z2 = psi
        Z3 = vx
        Z4 = vy
        xy = np.stack((X, Y), axis=-1)
        verts = xy[IEN]
        verts2 = xy[cc]
        triang = mtri.Triangulation(X,Y,IEN)
        fig, ax = plt.subplots(2,2)
        fig.subplots_adjust(hspace=0.5,top=0.8)
        surf = ax[0,0].tricontourf(triang,Z1,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                           X.max(), Y.min(), Y.max()))
        ax[0,0].set_title(r'$\omega$')
        fig.colorbar(surf, ax=ax[0,0])
        surf2 = ax[0,1].tricontourf(triang,Z2,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                           X.max(), Y.min(), Y.max()))
        ax[0,1].set_title(r'$\psi$')
        fig.colorbar(surf2, ax=ax[0,1])
        surf3 = ax[1,0].tricontourf(triang,Z3,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                           X.max(), Y.min(), Y.max()))
        ax[1,0].set_title(r'$v_{x}$')
        fig.colorbar(surf3, ax=ax[1,0])
        surf4 = ax[1,1].tricontourf(triang,Z4,50,cmap=matplotlib.cm.jet, extent=(X.min(),
                           X.max(), Y.min(), Y.max()))
        ax[1,1].set_title(r'$v_{y}$')
        fig.colorbar(surf4, ax=ax[1,1])
        fig.set_size_inches(9, 4)
        pc = matplotlib.collections.PolyCollection(verts,edgecolors=('lightgray',),
                                                      facecolor='None',
                                                      linewidths=(0.1,))
        # ax.add_collection(pc)
        pc2 = matplotlib.collections.PolyCollection(verts2,edgecolors=('black',),
                                                      facecolor='None',
                                                      linewidths=(.7,))
        pc2_c = copy.copy(pc2)
        pc2_cc = copy.copy(pc2_c)
        pc2_ccc = copy.copy(pc2_cc)
        pc2_cccc = copy.copy(pc2_ccc)
        pc_c = copy.copy(pc)
        pc_cc = copy.copy(pc_c)
        pc_ccc = copy.copy(pc_cc)
        pc_cccc = copy.copy(pc_ccc)
        ax[0,0].add_collection(pc_c)
        ax[0,0].add_collection(pc2_c)
        ax[0,1].add_collection(pc_cc)
        ax[0,1].add_collection(pc2_cc)
        ax[1,0].add_collection(pc_ccc)
        ax[1,0].add_collection(pc2_ccc)
        ax[1,1].add_collection(pc_cccc)
        ax[1,1].add_collection(pc2_cccc)
        #plt.plot(X,Y,'w.')
        fig.suptitle(r'Solução da Equação de Corrente-Vorticidade: $(\frac{M}{dt} + \nu K + v.G)w^{n+1} = (\frac{M}{dt})w^{n} $'+"\n\n")
        #plt.xlabel('comprimento da placa no eixo X [cm]')
        #plt.ylabel('comprimento da placa no eixo Y [cm]')
        #surf = plt.imshow(X,Y,Z, interpolation='quadric', origin='lower',
        #                  cmap=matplotlib.cm.jet, extent=(X.min(),
        #                  X.max(), Y.min(), Y.max()))
        
        # plt.clim(0, 300)
        # cbar = plt.colorbar(surf,shrink=1.0, aspect=20)
        # cbar.set_label('Temperatura [°C]')
        # labx = np.linspace(X.min(),X.max(),nx)
        # laby = np.linspace(Y.min(),Y.max(),ny)
        # plt.xticks(labx)
        # plt.yticks(laby)
        subname = '???'
        if(theta == 0):
            subname = 'resultados_explicito'
        if(theta == 1):
            subname = 'resultados_cilindro_implicito'
        if(theta == 1/2):
            subname = 'resultados_crank_nicholson'
        plt.savefig(os.path.join(name, subname, f'{index:04}.png'))
        plt.show()
        index += 1
        ## mean squared error MSE
        # e = np.sum((T - Tpast)**2)/T.shape[0]
    fp_in = f"./{name}/{subname}/*.png"
    fp_out = f"./{name}/{subname}/simulacao.gif"
    
    img, *imgs = [Image.open(f) for f in sorted(glob.glob(fp_in))]
    img.save(fp=fp_out, format='GIF', append_images=imgs,
              save_all=True, duration=dt*1000, loop=0)
solveWithTheta(1,0.01)
# msh = meshio.read('cavidade.msh')
# X = msh.points[:,0]
# Y = msh.points[:,1]
# IEN = msh.cells['triangle']
# npoints = len(X) 
# ne = IEN.shape[0]
# regions = msh.cell_data['triangle']['gmsh:geometrical']
# IENBound = cc = msh.cells['line']
# contorno,entrada,parede_top,saida,parede_bot = getContorno(IENBound)

# print(parede_top)  