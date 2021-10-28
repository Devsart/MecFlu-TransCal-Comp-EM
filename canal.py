# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 00:11:19 2021

@author: matheus.sartor
"""

import meshio
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import matplotlib
import os
from pathlib import Path
import glob
from PIL import Image
import meshio
import copy

name = 'funcao_corrente_vorticidade_mef_2d_triangular'
Path(os.path.join(name,'resultados_canal_v7_implicito')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_explicito_v3')).mkdir(parents=True, exist_ok=True)
Path(os.path.join(name,'resultados_crank_nicholson_v3')).mkdir(parents=True, exist_ok=True)


def montaKMG(X,Y,IEN):
	npoints = X.shape[0]
	K = np.zeros( (npoints,npoints), dtype='float' )
	M = np.zeros( (npoints,npoints), dtype='float' )
	Gx = np.zeros( (npoints,npoints), dtype='float' )
	Gy = np.zeros( (npoints,npoints), dtype='float' )

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

		area = (1/2.0)*np.linalg.det([[1, xi, yi], 
									 [1, xj, yj], 
									 [1, xk, yk]])

		ke_x = (1/(4.0*area))*np.array([[bi**2, bi*bj, bi*bk], 
										[bj*bi, bj**2, bj*bk], 
										[bk*bi, bk*bj, bk**2]])

		ke_y = (1/(4.0*area))*np.array([[ci**2, ci*cj, ci*ck], 
										[cj*ci, cj**2, cj*ck], 
										[ck*ci, ck*cj, ck**2]])
		ke = ke_x + ke_y  

		me = (area/12.0)*np.array([[2.0, 1.0, 1.0], 
							   [1.0, 2.0, 1.0], 
							   [1.0, 1.0, 2.0]])

		ge_x = (1/6.0)*np.array([[bi, bj, bk], 
								[bi, bj, bk], 
								[bi, bj, bk]])

		ge_y = (1/6.0)*np.array([[ci, cj, ck], 
								[ci, cj, ck], 
								[ci, cj, ck]])


		for i_loc in range(0,3):
			i_glb = IEN[elem,i_loc]
			for j_loc in range(0,3):
				j_glb = IEN[elem,j_loc]
				
				K[i_glb, j_glb] += ke[i_loc, j_loc]
				M[i_glb, j_glb] += me[i_loc, j_loc]
				Gx[i_glb, j_glb] += ge_x[i_loc, j_loc]
				Gy[i_glb, j_glb] += ge_y[i_loc, j_loc]

	return K, M, Gx, Gy

def montaKest(X, Y, IEN, vx, vy, dt):

	npoints = X.shape[0]
	ne = len(IEN)

	Kest = np.zeros( (npoints,npoints), dtype='float' )

	for elem in range(ne):

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

		vxm = (vx[v1] + vx[v2] + vx[v3])/3
		vym = (vy[v1] + vy[v2] + vy[v3])/3

		area = (1/2.0)*np.linalg.det([[1, xi, yi], 
									 [1, xj, yj], 
									 [1, xk, yk]])

		ke_x = (1/(4.0*area))*np.array([[bi**2, bi*bj, bi*bk], 
										[bj*bi, bj**2, bj*bk], 
										[bk*bi, bk*bj, bk**2]])

		ke_y = (1/(4.0*area))*np.array([[ci**2, ci*cj, ci*ck], 
										[cj*ci, cj**2, cj*ck], 
										[ck*ci, ck*cj, ck**2]])
		ke_xy = ke_x + ke_y  

		ke_est = vxm * dt/2 * (vxm*ke_x + vym*ke_xy) + vym * dt/2 * (vxm*ke_xy + vym*ke_y)

		for i_loc in range(0,3):
			i_glb = IEN[elem,i_loc]
			for j_loc in range(0,3):
				j_glb = IEN[elem,j_loc]
				
				Kest[i_glb, j_glb] += ke_est[i_loc, j_loc]
	return Kest

arquivo = 'canal'


msh = meshio.read(f'{arquivo}.msh')
X = msh.points[:, 0]
Y = msh.points[:, 1]
IEN = msh.cells['triangle']  # malha bidimensional (triangulos)
IENbound = msh.cells['line']  # malha de contorno (IEN de seguimento de retas)
# IENboundTypeElem = list(msh.cell_data['line']['gmsh:physical'] - 1)
# boundNames = list(msh.field_data.keys())
# IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]

npoints = len(X)  # obtendo o numero de pontos
ne = IEN.shape[0]  # obtendo o numero de elementos

parede_top = [0,5]
for i in range(13,52):
    parede_top.append(i)
parede_bot = [4,3,2,1]
for i in range(67,116):
    parede_bot.append(i)
entrada = []
for i in range(6,13):
    entrada.append(i)
saida = []
for i in range(52,67):
    saida.append(i)
contorno = [*entrada,*parede_top,*parede_bot,*saida]

K, M, Gx, Gy = montaKMG(X,Y,IEN)

theta = 1
t = 0
index = 0
Re = 30
dt = 0.05
sigma = 1/Re
nIter = 40

w_z = np.zeros((npoints), dtype = 'float64')
psi = np.zeros(npoints, dtype = 'float64')
vx = np.zeros((npoints), dtype = 'float64')
vy= np.zeros((npoints), dtype = 'float64')
T = np.zeros((npoints), dtype = 'float64')

for i in parede_top:
	vx[i] = 0.0
	vy[i] = 0.0
for i in parede_bot:
	vx[i] = 0.0
	vy[i] = 0.0
for i in entrada:
	vx[i] = 1.0
	vy[i] = 0.0
	# if ccName[i] == 'right':
	#   vx[i] = 0.0
	#   vy[i] = 0.0
# for i in saida:
# 	vx[i] = 0.0
# 	vy[i] = 0.0

w_z = np.linalg.solve(M, (Gx@vy - Gy@vx))
I = np.identity(npoints)

for n in range(nIter):
	t += dt
	w_zc = np.linalg.solve(M, (Gx@vy - Gy@vx))

	# Kest = montaKest(X, Y, IEN, vx, vy, dt)

	A_wz = M.copy()/dt + sigma*K.copy() + ((vx*I)@Gx + (vy*I)@Gy)
	b_wz = M.copy()/dt @ w_z

	for i in contorno:
		A_wz[i,:] = 0.0
		A_wz[i,i] = 1.0
		b_wz[i] = w_zc[i]

	w_z = np.linalg.solve(A_wz, b_wz)

	A_psi = K.copy()
	b_psi = M.copy() @ w_z

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
		b_psi[i] = 2*(Y[i]-0.5)
# 	for i in saida:
# 		A_psi[i,:] = 0.0
# 		A_psi[i,i] = 1.0
# 		b_psi[i] = 0.0

	psi = np.linalg.solve(A_psi, b_psi)

	vx =  np.linalg.solve(M, (Gy@psi))
	vy = np.linalg.solve(M, -(Gx@psi))

	for i in parede_top:
		vx[i] = 0.0
		vy[i] = 0.0
	for i in parede_bot:
		vx[i] = 0.0
		vy[i] = 0.0
	for i in entrada:
		vx[i] = 1.0
		vy[i] = 0.0
	# if ccName[i] == 'right':
	#   vx[i] = 0.0
# 	#   vy[i] = 0.0
# 	for i in saida:
# 		vx[i] = 0.0
# 		vy[i] = 0.0
	Z1 = w_zc
	Z2 = psi
	Z3 = vx
	Z4 = vy
	xy = np.stack((X, Y), axis=-1)
	verts = xy[IEN]
	verts2 = xy[IENbound]
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
	fig.suptitle(r'Solução da Equação de Corrente-Vorticidade: $(\frac{M}{dt} + \nu K + v.G)w^{n+1} = (\frac{M}{dt})w^{n} $'+"\n"+f'Re = {Re}; dt = {dt}s; t = {t}s')
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
		subname = 'resultados_explicito_v3'
	if(theta == 1):
		subname = 'resultados_canal_v7_implicito'
	if(theta == 1/2):
		subname = 'resultados_crank_nicholson_v3'
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
	# A_T = M.copy()/dt + sigma*K.copy() + ((vx*I)@Gx + (vy*I)@Gy)
	# b_T = M.copy()/dt @ T

	# for i in cc:
	# 	if ccName[i] == 'hole':
	# 		A_T[i,:] = 0.0
	# 		A_T[i,i] = 1.0
	# 		b_T[i] = 100

	# T = np.linalg.solve(A_T,b_T)

print(n)