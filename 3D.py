
import functools
import itertools as IT
import numpy as np
import scipy.ndimage as ndimage
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

n=20 # ligne
m=n # colonne
M=np.eye(n**3)
D= 10**(-10)
Lx=200*10**(-9) ### profondeur de la couche
Nx=n ## nbr des points
dx=Lx/(Nx-1)
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
for i in range(n**2):
    X0[i]=10
col=X0

iM=np.linalg.inv(M)
iM=np.linalg.inv(M)
col =np.dot(iM,ksi(col))
col=np.reshape(col,(n,n,n))
def cartesian_product_broadcasted(*arrays):
    """
    http://stackoverflow.com/a/11146645/190597 (senderle)
    """
    broadcastable = np.ix_(*arrays)
    broadcasted = np.broadcast_arrays(*broadcastable)
    dtype = np.result_type(*arrays)
    rows, cols = functools.reduce(np.multiply, broadcasted[0].shape), len(broadcasted)
    out = np.empty(rows * cols, dtype=dtype)
    start, end = 0, rows
    for a in broadcasted:
        out[start:end] = a.reshape(-1)
        start, end = end, end + rows
    return out.reshape(cols, rows).T

# @profile  # used with `python -m memory_profiler script.py` to measure memory usage
fig , ax = plt.subplots()
ax = fig.add_subplot(1, 1, 1, projection='3d')
x, y, z = cartesian_product_broadcasted(*[np.arange(n, dtype='int16')]*3).T
surf=ax.scatter(z, y, x, c=col, cmap=plt.get_cmap('plasma'),vmax=0.1,vmin=0)
fig.colorbar(surf, ax=ax)

def update(i):
    ax.clear()
    global col
    col = np.reshape(col, (n * n * n, 1))
    col = np.dot(iM, ksi(col))
    col = np.reshape(col, (n, n, n))
    x, y, z = cartesian_product_broadcasted(*[np.arange(n, dtype='int16')] * 3).T
    col = col.ravel()
    ax.scatter(z, y, x, c=col, cmap=plt.get_cmap('plasma'),vmin=0,vmax=0.1)
    plt.title(str(i+1))
ani = animation.FuncAnimation(fig, update,frames=500,interval=20,repeat = False)
plt.show()