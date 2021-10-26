import meshio
import numpy as np
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

K, M, Gx, Gy = montaKMG(X,Y,IEN)

c1 = 0.0
c2 = 1.0
delta_c = c2 - c1
delta_y = max(Y) - min(Y)

b = np.zeros( (npoints), dtype = 'float')
A = K.copy()

for i in cc:
	A[i,:] = 0.0
	A[i,i] = 1.0
	if ccName[i] == 'top':
		b[i] = c2
	if ccName[i] == 'bottom':
		b[i] = c1
	if ccName[i] == 'left':
		b[i] = Y[i]
	if ccName[i] == 'right':
		b[i] = Y[i]
	if ccName[i] == 'hole':
		b[i] = (max(Y)-min(Y))/2

psi = np.linalg.solve(A,b)
vx =  np.linalg.solve(M,(Gy@psi))
vy = np.linalg.solve(M,(-Gx@psi))

# point_data = {'potencial': psi}
# meshio.write_points_cells('solucao.vtk', msh.points, msh.cells,point_data=point_data)

triang = mtri.Triangulation(X,Y,IEN)
ax = plt.axes()
ax.set_aspect('equal')
ff = ax.tricontourf(triang,psi,200,cmap='jet')
plt.colorbar(ff)
plt.show()