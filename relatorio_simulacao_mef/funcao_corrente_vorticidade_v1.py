import meshio
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

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

arquivo = 'degrau'

msh = meshio.read(f'{arquivo}.msh')
X = msh.points[:, 0]
Y = msh.points[:, 1]
IEN = msh.cells['triangle']  # malha bidimensional (triangulos)
IENbound = msh.cells['line']  # malha de contorno (IEN de seguimento de retas)
IENboundTypeElem = list(msh.cell_data['line']['gmsh:physical'] - 1)
boundNames = list(msh.field_data.keys())
IENboundElem = [boundNames[elem] for elem in IENboundTypeElem]

npoints = len(X)  # obtendo o numero de pontos
ne = IEN.shape[0]  # obtendo o numero de elementos

# Cria lista de nos do contorno
cc = np.unique(IENbound.reshape(IENbound.size))  # lista de nos do contorno
ccName = [[] for i in range(len(X))]  # lista com os nomes das cc's
for elem in range(0, len(IENbound)):
	ccName[IENbound[elem][0]] = IENboundElem[elem]
	ccName[IENbound[elem][1]] = IENboundElem[elem]

K, M, Gx, Gy = montaKMG(X,Y,IEN)

dt = 0.05
sigma = 1/30
nIter = 200

w_z = np.zeros((npoints), dtype = 'float64')
psi = np.zeros(npoints, dtype = 'float64')
vx = np.zeros((npoints), dtype = 'float64')
vy= np.zeros((npoints), dtype = 'float64')
T = np.zeros((npoints), dtype = 'float64')

for i in cc:
	if ccName[i] == 'top':
		vx[i] = 0.0
		vy[i] = 0.0
	if ccName[i] == 'bottom':
		vx[i] = 0.0
		vy[i] = 0.0
	if ccName[i] == 'left':
		vx[i] = 1.0
		vy[i] = 0.0
	# if ccName[i] == 'right':
	#   vx[i] = 0.0
	#   vy[i] = 0.0
	# if ccName[i] == 'hole':
	#   vx[i] = 0.0
	#   vy[i] = 0.0

w_z = np.linalg.solve(M, (Gx@vy - Gy@vx))
I = np.identity(npoints)

for n in range(nIter):

	w_zc = np.linalg.solve(M, (Gx@vy - Gy@vx))

	# Kest = montaKest(X, Y, IEN, vx, vy, dt)

	A_wz = M.copy()/dt + sigma*K.copy() + ((vx*I)@Gx + (vy*I)@Gy)
	b_wz = M.copy()/dt @ w_z

	for i in cc:
		A_wz[i,:] = 0.0
		A_wz[i,i] = 1.0
		b_wz[i] = w_zc[i]

	w_z = np.linalg.solve(A_wz, b_wz)

	A_psi = K.copy()
	b_psi = M.copy() @ w_z

	for i in cc:
		if ccName[i] == 'top':
			A_psi[i,:] = 0.0
			A_psi[i,i] = 1.0
			b_psi[i] = 1.0
		if ccName[i] == 'bottom':
			A_psi[i,:] = 0.0
			A_psi[i,i] = 1.0
			b_psi[i] = 0.0
		if ccName[i] == 'left':
			A_psi[i,:] = 0.0
			A_psi[i,i] = 1.0
			b_psi[i] = Y[i]
		# if ccName[i] == 'right':
			# A_psi[i,:] = 0.0
			# A_psi[i,i] = 1.0
		#   b_psi[i] = 0
		# if ccName[i] == 'hole':
		# 	A_psi[i,:] = 0.0
		# 	A_psi[i,i] = 1.0
		# 	b_psi[i] = 0.5

	psi = np.linalg.solve(A_psi, b_psi)

	vx =  np.linalg.solve(M, (Gy@psi))
	vy = np.linalg.solve(M, -(Gx@psi))

	for i in cc:
		if ccName[i] == 'top':
			vx[i] = 0.0
			vy[i] = 0.0
		if ccName[i] == 'bottom':
			vx[i] = 0.0
			vy[i] = 0.0
		if ccName[i] == 'left':
			vx[i] = 1.0
			vy[i] = 0.0
		# if ccName[i] == 'right':
			# vx[i] = 0.0
			# vy[i] = 0.0
		# if ccName[i] == 'hole':       
		#   vx[i] = 0.0
		#   vy[i] = 0.0

	# A_T = M.copy()/dt + sigma*K.copy() + ((vx*I)@Gx + (vy*I)@Gy)
	# b_T = M.copy()/dt @ T

	# for i in cc:
	# 	if ccName[i] == 'hole':
	# 		A_T[i,:] = 0.0
	# 		A_T[i,i] = 1.0
	# 		b_T[i] = 100

	# T = np.linalg.solve(A_T,b_T)

	print(n)

point_data = {'potencial': w_z}
meshio.write_points_cells(f'{arquivo}_wz.vtk', msh.points, msh.cells,point_data=point_data)
point_data = {'potencial': psi}
meshio.write_points_cells(f'{arquivo}_psi.vtk', msh.points, msh.cells,point_data=point_data)
point_data = {'potencial': vx}
meshio.write_points_cells(f'{arquivo}_vx.vtk', msh.points, msh.cells,point_data=point_data)
point_data = {'potencial': vy}
meshio.write_points_cells(f'{arquivo}_vy.vtk', msh.points, msh.cells,point_data=point_data)
# point_data = {'potencial': T}
# meshio.write_points_cells(f'{arquivo}_T.vtk', msh.points, msh.cells,point_data=point_data)

# triang = mtri.Triangulation(X,Y,IEN)
# ax = plt.axes()
# ax.set_aspect('equal')
# f1 = ax.tricontourf(triang, vx, 200, cmap='jet')
# plt.show()