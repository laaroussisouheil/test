import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import random
### la bijection et sa reciproque qui transforme carte <----------> vecteur colonne
def bij2D(i,j,n):
    return (i+(j-1)*n)
def bij2D_1(k,n) :
    j=k//n + 1
    i=k-(j-1)*n
    return(i,j)
### Initialisation (valeurs aleatoires)
n=15 # ligne
m=n # colonne
M=np.eye(n*m)
D= 2.9*10**(-9)
Lx=200*10**(-9) ### profondeur de la couche
Nx=n ## nbr des points
dx=Lx/(Nx-1)
dt= float ((dx**2)/(1*D))
alpha = dt*D/(2*(dx**2))
beta = 2*(1+2*alpha)
gamma = 2*(1-2*alpha)
point_inter=[]
### boucle qui donne la liste K des pts qui ont 4 voisins
for j in range(2,n) :
    for i in range(2,m) :
        point_inter.append(bij2D(i,j,n))
### boucle qui genere la matrice
for k in point_inter :
    M[k-1,k-1]=beta
    M[k-1,k]=-alpha
    M[k-1,k-2]=-alpha
    M[k-1,k+n-1]=-alpha
    M[k-1,k-n-1]=-alpha
for k in range (n) :
    M[k][k+n]=-1
for k in range (1,n) :
    M[k*n][k*n+1]=-1
    M[n-1+k*n][n-2+k*n]=-1
### implementation des CL
### la bijection ksi entre [1,n]*[1,n] --------->[1,n]*[1,n] qui verifie M*X(t+1)=ksi(X0)

def ksi(n,X):
    Ximg=np.zeros((n*n,1))
    for k in point_inter:
        Ximg[k-1]=alpha*(X[k-2]+X[k]+X[k+n-1]+X[k-n-1])+gamma*X[k-1]
    return(Ximg)
X0=np.zeros((n*n,1))
for i in range(n**2-1):
    X0[i]=5 ##random.uniform(0,0.2)