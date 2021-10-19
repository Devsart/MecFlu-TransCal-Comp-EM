## =================================================================== ##
#  this is file calor2d.py, created at 24-Aug-2021                #
#  maintained by Gustavo R. Anjos                                       #
#  e-mail: gustavo.rabello@coppe.ufrj.br                                #
## =================================================================== ##


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

Lx = 1.0
Ly = 1.0
Lz = 1.0
nx = 4
ny = 4
nz = 4
npoints = nx*ny*nz
ne = (nx-1)*(ny-1)*(nz-1)
dx = Lx/(nx-1)
dy = Ly/(ny-1)
dz = Lz/(nz-1)

# geracao de malha
X = np.array( [0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1, 
               0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1, 
               0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1, 
               0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1,0.0,1/3,2/3,1] )

Y = np.array( [0.0,0.0,0.0,0.0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,2/3,1,1,1,1, 
               0.0,0.0,0.0,0.0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,2/3,1,1,1,1,
               0.0,0.0,0.0,0.0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,2/3,1,1,1,1,
               0.0,0.0,0.0,0.0,1/3,1/3,1/3,1/3,2/3,2/3,2/3,2/3,1,1,1,1] )
Z = np.array( [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
               1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,
               2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,2/3,
               1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1] )


cc1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]
cc2 = [7,11]
cc3 = [12,13,14,15]
cc4 = [4,8]
cc5 = [6,7,9,10]
cc6 = [2,14]
cc = cc1 + cc2 + cc3 + cc4 + cc5 + cc6
inner = [21,22,25,26,37,38,41,42]

#--------------------------------------------------
# plt.plot(X,Y,'bo')
# plt.plot(X[cc],Y[cc],'ko')
# plt.show()
#-------------------------------------------------- 

bc = 0*np.ones( (npoints),dtype='float' )
for i in cc1:
 bc[i] = X[i]
for i in cc2:
 bc[i] = Y[i]*Y[i] + 1
for i in cc3:
 bc[i] = X[i]*X[i] + 1
for i in cc4:
 bc[i] = Y[i]
for i in cc5:
 bc[i] = Y[i]
for i in cc6:
 bc[i] = Y[i]

#--------------------------------------------------
# # # plot 2D cor (quadrilatero)
# Z = bc.reshape(ny,nx)
# surf = plt.imshow(Z, interpolation='quadric', origin='lower',
#                   cmap=matplotlib.cm.jet, extent=(X.min(),
#                   X.max(), Y.min(), Y.max()))
# plt.colorbar(surf,shrink=1.0, aspect=20)
# plt.grid(color='black', linestyle='solid', linewidth=0.5)
# labx = np.linspace(X.min(),X.max(),nx)
# laby = np.linspace(Y.min(),Y.max(),ny)
# plt.xticks(labx)
# plt.yticks(laby)
# plt.show()
#-------------------------------------------------- 

A = np.zeros( (npoints,npoints),dtype='float' )
b = np.zeros( (npoints),dtype='float' )

# Eq. da cond. de contorno de Dirichlet
for i in cc:
 A[i,i] = 1.0
 b[i]   = bc[i]

# Eqs. pontos internos (inner):
for i in inner:
 A[i,i-1] = 1/(dx*dx)                       # direita
 A[i,i] = -2/(dx*dx) -2/(dy*dy) -2/(dz*dz)  # ele mesmo
 A[i,i+1] = 1/(dx*dx)                       # esquerda
 A[i,i-nx] = 1/(dy*dy)                      # baixo
 A[i,i+nx] = 1/(dy*dy)                      # cima
 A[i,i-ny] = 1/(dz*dz)                            # frente
 A[i,i-ny] = 1/(dz*dz)                            # tras

# solucao do sistema linear Ax=b
Ainv = np.linalg.inv(A)
T = Ainv@b

#T = np.linalg.solve(A,b)

# plot 3D (plano) cor (cubo)
Z = T.reshape(ny,nx)
surf = plt.imshow(Z, interpolation='quadric', origin='lower',
                  cmap=matplotlib.cm.jet, extent=(X.min(),
                  X.max(), Y.min(), Y.max()))
plt.colorbar(surf,shrink=1.0, aspect=20)
plt.grid(color='black', linestyle='solid', linewidth=0.5)
labx = np.linspace(X.min(),X.max(),nx)
laby = np.linspace(Y.min(),Y.max(),ny)
plt.xticks(labx)
plt.yticks(laby)
plt.show()
