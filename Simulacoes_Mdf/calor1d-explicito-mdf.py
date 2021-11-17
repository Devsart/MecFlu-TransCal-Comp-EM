## =================================================================== ##
#  this is file calor1d-explicito-mdf.py, created at 10-Aug-2021        #
#  maintained by Gustavo Rabello dos Anjos                              #
#  e-mail: gustavo.rabello@gmail.com                                    #
## =================================================================== ##

# tarefas:
# 1) usar while no loop do tempo
# 2) verificar qual eh a relacao de dt e dx para que a simulacao 
# funcione corretamente

import numpy as np
import matplotlib.pyplot as plt

# parametros da simulacao
L = 1.0
npoints = 10
ne = npoints-1
dx = L/ne
Q = 0.0 # fonte de calor
k = 1.0  # condutividade termica do material
cv = 1.0 # capacidade termica
rho = 1.0 # densidade
dt = 0.001 # passo de tempo
time = 0.0 # tempo de simulacao

# condicao de contorno de Dirichlet (espacial)
Te = 100.0
Td = 0.0

# condicao de inicial (temporal)
T = np.zeros( (npoints),dtype='float' )
T[0] = Te
T[-1] = Td

# geracao dos pontos
X = np.linspace(0,L,npoints)

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

# plot da condicao inicial (time=0.00)
plt.plot(X,T,'r-')

# Eq. do ponto i (pontos do miolo) para time=0.01
for n in range(0,10000):
 Tp = T[1]
 for i in range(1,npoints-1):
  T[i] = T[i] + \
         ( dt/(dx*dx) ) * (k/(rho*cv)) * ( T[i+1] - 2*T[i] + T[i-1] ) + \
         ( dt/(rho*cv) ) * Q
 Tf = T[1]
 if (Tf-Tp) < 1E-8:
  break
print (n,Tf-Tp)

plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.show()
