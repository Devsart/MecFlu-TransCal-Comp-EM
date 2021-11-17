# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 01:57:14 2021

@author: matheus.sartor
"""
import sympy as sp
import numpy as np
import matplotlib.pyplot as plt

x = sp.symbols('x')
def Galerkin1D(X):
    k = np.zeros((len(X),len(X)))
    for i in range(len(X)-1):
        N0 = -(1/(X[i+1]-X[i]))*x+1
        N1 = (1/(X[i+1]-X[i]))*x
        dN0 = sp.diff(N0,x)
        dN1 = sp.diff(N1,x)
        N00 = sp.integrate(dN0*dN0,(x,X[0],X[1]))
        N01 = sp.integrate(dN0*dN1,(x,X[0],X[1]))
        N10 = sp.integrate(dN1*dN0,(x,X[0],X[1]))
        N11 = sp.integrate(dN1*dN1,(x,X[0],X[1]))
        k[i][i] += N00
        k[i][i+1] += N01
        k[i+1][i] += N10
        k[i+1][i+1] += N11
    return k

def Drichlet1D(k,Ti,Tf):
    k[0][0] = 1
    k[0][1] = 0
    k[-1][-1] = 1
    k[-1][-2] = 0
    b = np.zeros((k.shape[0],1))
    b[0][0] = Ti
    b[-1][-1] = Tf
    return k,b

def SolveLinear1D(X,Ti,Tf):
    k_1 = Galerkin1D(X)
    k_2,b = Drichlet1D(k_1,Ti,Tf)
    print(k_2,b)
    T = np.linalg.solve(k_2,b)
    return T

X = [0.,.167,.334,.5,1.]
T = SolveLinear1D(X,100,0)

plt.plot(X,T,'ko-')
plt.xlabel('comprimento da barra [m]')
plt.ylabel('temperatura [oC]')
plt.show()


