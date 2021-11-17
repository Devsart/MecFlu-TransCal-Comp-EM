## =================================================================== ##
#  this is file calor2d.py, created at 24-Aug-2021                #
#  maintained by Gustavo R. Anjos                                       #
#  e-mail: gustavo.rabello@coppe.ufrj.br                                #
## =================================================================== ##


import numpy as np
import matplotlib.pyplot as plt
import matplotlib

Lx = 1
Ly = 1
nx = 4
ny = 4
npoints = nx*ny
ne = (nx-1)*(ny-1)
dx = 1.0/ne
dy = 1.0/ne

# geracao de pontos da malha
Xv = np.linspace(0,Lx,nx)
Yv = np.linspace(0,Ly,ny)
X,Y = np.meshgrid(Xv,Yv)
X = X.reshape(nx*ny)
Y = Y.reshape(nx*ny)

# captura dos pontos de contorno (automatizar)
cc1 = [0,1,2,3]
cc2 = [7,11]
cc3 = [12,13,14,15]
cc4 = [4,8]
cc = cc1 + cc2 + cc3 + cc4
inner = [5,6,9,10]

#--------------------------------------------------
# plt.plot(X,Y,'bo')
# plt.plot(X[cc],Y[cc],'ko')
# plt.show()
#-------------------------------------------------- 

# condicoes de contorno espaciais
bc = 0*np.ones( (npoints),dtype='float' )
for i in cc1:
 bc[i] = X[i]
for i in cc2:
 bc[i] = Y[i]*Y[i] + 1
for i in cc3:
 bc[i] = X[i]*X[i] + 1
for i in cc4:
 bc[i] = Y[i]

# condicao inicial (temporal)
T = 0*np.ones( (npoints),dtype='float' )
for i in cc:
 T[i] = bc[i]

# distribuicao de fonte de calor
Q = 0*np.ones( (npoints),dtype='float' )
for i in range(0,npoints):
 # cpu 1
 if X[i]<0.4 and X[i] > 0.2 and Y[i] < 0.4 and Y[i]>0.2: 
  Q[i] = 100.0
#--------------------------------------------------
#  # cpu 2
#  # cpu 3
#  # cpu 4
#-------------------------------------------------- 


# # plot 2D cor (quadrilatero) - condicao inicial do problema
Z = Q.reshape(ny,nx)
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
