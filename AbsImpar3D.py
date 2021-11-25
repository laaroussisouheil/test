
import numpy as np

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

def bij3D(i,j,k,n):
    return (i+(j-1)*n+(k-1)*(n*n))
def bij3D_1(l,n) :
    k=l//(n*n) + 1
    j=(l-(k-1)*n*n)//n + 1
    i=l-(j-1)*n-(k-1)*n*n
    return(i,j)
N=2000
n=4 # ligne
m=n # colonne
M=np.eye(n**3)
D= 10**(-10)
Lx=200*10**(-9) ### profondeur de la couche
Nx=n ## nbr des points
dx=Lx/(Nx-1)
P=0.5
PRO=dx*P/(1-P)
dt= (dx**2)/(6*D)
alpha = dt*D/((dx**2))
beta =(1+3*alpha)*2
gamma =(1-3*alpha)*2
point_inter=[]
### boucle qui donne la liste K des pts qui ont 6 voisins
for k in range(2,n):
    for j in range(2,n) :
        for i in range(2,n) :
            point_inter.append(bij3D(i,j,k,n))
### boucle qui genere la matrice
for k in point_inter :
    M[k-1,k-1]=beta
    M[k-1,k]=-alpha
    M[k-1,k-2]=-alpha
    M[k-1,k+n-1]=-alpha
    M[k-1,k-n-1]=-alpha
    M[k-1,k+n*n-1]=-alpha
    M[k-1,k-n*n-1]=-alpha

for k in range (n**2):
    M[k][k+n**2]=-1
    M[k+(n-1)*(n**2)][k+(n-1)*(n**2)]=1+dx/PRO
    M[k+(n-1)*(n**2)][k+(n-2)*(n**2)]=-1

for k in range(1,n-1):
    for i in range (n-1):
        M[i+k*(n**2)][i+n+k*(n**2)]=-1
    for i in range(0,n-1):
        M[i+(n-1)*n+k*(n**2)][i+(n-2)*n+k*(n**2)]=-1

for k in range(1,n-1):
    for j in range (1,n-1):
        M[j*n+k*(n**2)][1+j*n+k*(n**2)]=-1
    for j in range( n):
        M[n-1+j*n+k*(n**2)][n-2+j*n+k*(n**2)]=-1
    
def ksi(X):
    Ximg=np.zeros((n*n*n,1))
    for k in point_inter:
        Ximg[k-1]=alpha*(X[k-2]+X[k]+X[k+n-1]+X[k-n-1]+X[k+n**2-1]+X[k-n**2-1])+gamma*X[k-1]
    return(Ximg)

X0=np.zeros((n*n*n,1))
for i in range(n**3):
    X0[i]=10
col=X0
k=np.sum(col)
y=np.zeros(N)
iM=np.linalg.inv(M)
for p in range(N):
     
        col = ksi(col)
        col = np.dot(iM, col)
        z=np.sum(col)
        y[p]=k-z
#z=np.reshape(z,(n,))
x=[]
for j in range(N):
    x.append(j*dt)
np.savetxt('0-5,3D.txt',y,fmt='%.2f')
plt.title("Concentration totale")
plt.xlabel("temps")
plt.ylabel("Concentration")
plt.plot(x, y)
###plt.show()