## =================================================================== ##
#  this is file calor1d-mdf.py, created at 10-Aug-2021                #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##


import numpy as np
import matplotlib.pyplot as plt

# parametros da simulacao
L = 1.0
npoints = 50
ne = npoints-1
dx = L/ne
#Q = 0.0 # fonte de calor
#k = 1.0  # condutividade termica do material

# cond. termica e fonte de calor variaveis
k = np.ones( (npoints),dtype='float' )
Q = np.ones( (npoints),dtype='float' )
X = np.linspace(0,L,npoints)
#--------------------------------------------------
for i in range(0,npoints):
    if X[i] > 0.3:
        k[i] = 0.01
    if X[i] > 0.5:
        Q[i] = 3.0
#-------------------------------------------------- 
  



# condicao de contorno de Dirichlet
Te = 1.0
Td = 0.0

# geracao dos pontos

# geracao da matriz de conectividade
IEN = np.zeros( (ne,2),dtype='int' )
for e in range(0,ne):
 IEN[e] = [e,e+1]

# vetor de indices de contorno
cc = [0,npoints-1]
# vetor dos valores do contorno
bcc = np.zeros( (npoints),dtype='float' )
bcc[0] = Te
bcc[-1] = Td

#--------------------------------------------------
# plt.plot(X,0*X,'ko-')
# plt.plot(X[cc],0*X[cc],'ro')
# plt.show()
#-------------------------------------------------- 

# inicializar a matriz A e o vetor b
A = np.zeros( (npoints,npoints),dtype='float' )
b = np.zeros( (npoints),dtype='float' )
# populando os valores dos pontos internos
for i in range(1,npoints-1):
 A[i,i]   = -2.0/(dx*dx) # diagonal principal
 A[i,i-1] =  1.0/(dx*dx) # diagonal inferior
 A[i,i+1] =  1.0/(dx*dx) # diagonal superior
 b[i]     = -Q[i]/k[i]


# populando os valores dos pontos de contorno
for i in cc:
 A[i,i] = 1.0
 b[i]   = bcc[i]

# solucao do sistema linear Ax=b
#Ainv = np.linalg.inv(A)
#T = Ainv@b

T = np.linalg.solve(A,b)

plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.show()
