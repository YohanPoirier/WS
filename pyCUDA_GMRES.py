import numpy as np
import time
import pycuda.autoinit
import skcuda.cublas as cublas
from pycuda import gpuarray as ga
import random as rd

def product_Ax(A_d, X, N, handler):
    
    X_d = ga.to_gpu(X.reshape((N,1)).astype(RP))
    res_d = ga.to_gpu(np.zeros((N,1)).astype(RP))
    
    alpha = 1
    if precision == 1 :
        cublas.cublasSgemv(handler, "t", N, N, alpha, A_d.gpudata, N, X_d.gpudata, 1, 0, res_d.gpudata, 1)
    else :
        cublas.cublasDgemv(handler, "t", N, N, alpha, A_d.gpudata, N, X_d.gpudata, 1, 0, res_d.gpudata, 1)

    result = res_d.get()
    
    return result[:,0]
    
    
    
def rotation(H,cs,sn,k):
    
    if H[k,k] == 0 :
        cs[k] = 0
        sn[k] = 1
    else:
        cs[k] = abs(H[k,k])/np.sqrt(H[k,k]**2 + H[k+1,k]**2)
        sn[k] = cs[k]*H[k+1,k]/H[k,k]
        
    
def apply_rotation(H,cs,sn,k):
    
    for i in range(k-1):
        temp = cs[i]*H[i,k] + sn[i]*H[i+1,k]
        H[i+1,k] = -sn[i]*H[i,k] + cs[i]*H[i+1,k]
        H[i,k] = temp
        
    rotation(H,cs,sn,k)
    
    H[k,k] = cs[k]*H[k,k] + sn[k]*H[k+1,k]
    H[k+1,k] = 0
  
def Arnoldi_d(Q, A_d, H, k, handler):
    
    N = Q.shape[0]
    ta = time.time()
    v = product_Ax(A_d, Q[:,k], N, handler)

    tb = time.time()
    for i in range(k+1):
        H[i,k] = np.dot(Q[:,i], v)
        v = v - H[i,k]*Q[:,i]
    
    tc = time.time()
    H[k+1,k] = np.linalg.norm(v)
    
    Q[:, k+1] = v/H[k+1,k]
    td = time.time()
    
    print(tb-ta,tc-tb,td-tc)
    
    
def Arnoldi(Q, A, H, k):
    

    ta = time.time()
    v = np.matmul(A, Q[:,k])
    
    tb = time.time()
    for i in range(k+1):
        H[i,k] = np.dot(Q[:,i], v)
        v = v - H[i,k]*Q[:,i]
    
    tc = time.time()
    H[k+1,k] = np.linalg.norm(v)
    
    Q[:, k+1] = v/H[k+1,k]
    td = time.time()
    
    print(tb-ta,tc-tb,td-tc)    

def GMRES(A, B, X0, tol, Nmax) :
    
    
    # Initialisation
    N = A.shape[0]
    Q = np.zeros((N, Nmax+1))
    H = np.zeros((Nmax+1,Nmax+1))
    cs = np.zeros(Nmax)
    sn = np.zeros(Nmax)
    e1 = np.zeros(Nmax+1)
    
    # Residu initial
    Ax = np.dot(A,X0)
    r0 = B - Ax
    Q[:,0] = r0/np.linalg.norm(r0)
    e1[0] = 1
    
    error = np.linalg.norm(r0)/np.linalg.norm(B)
    beta = np.linalg.norm(r0)*e1
    

    k = 0
    
    while (error > tol) and k < min(Nmax,N) - 1 :
        
        Arnoldi(Q, A, H, k)

        apply_rotation(H, cs, sn, k)
        beta[k+1] = -sn[k]*beta[k]
        beta[k] = cs[k]*beta[k]
        
        error = abs(beta[k+1])/np.linalg.norm(B)

        k+= 1
    
    
    print("Nombre d'iterations : {}".format(k))
    
    Hinv = np.linalg.inv(H[:k,:k])
    
    y = np.matmul(Hinv[:k,:k],beta[:k])
    
    X = np.dot(Q[:, :k], y) + X0
    
    td = time.time()
    
    
    return X

def GMRES_d(A_d, B, X0, tol, Nmax) :
    
    handler = cublas.cublasCreate()

    # Initialisation
    N = B.shape[0]
    Q = np.zeros((N, Nmax+1))
    H = np.zeros((Nmax+1,Nmax+1))
    cs = np.zeros(Nmax)
    sn = np.zeros(Nmax)
    e1 = np.zeros(Nmax+1)
    
    # Residu initial
    Ax = product_Ax(A_d, X0, N, handler)
    
    print(Ax[:10])
    r0 = B - Ax
    Q[:,0] = r0/np.linalg.norm(r0)
    e1[0] = 1
    
    error = np.linalg.norm(r0)/np.linalg.norm(B)
    beta = np.linalg.norm(r0)*e1
    

    k = 0
    
    while (error > tol) and k < min(Nmax,N) - 1 :
        
        Arnoldi_d(Q, A_d, H, k, handler)

        apply_rotation(H, cs, sn, k)
        beta[k+1] = -sn[k]*beta[k]
        beta[k] = cs[k]*beta[k]
        
        error = abs(beta[k+1])/np.linalg.norm(B)

        k+= 1
    
    
    print("Nombre d'iterations : {}".format(k))
    
    Hinv = np.linalg.inv(H[:k,:k])
    
    y = np.matmul(Hinv[:k,:k],beta[:k])
    
    X = np.dot(Q[:, :k], y) + X0
    
    td = time.time()
    
    
    return X
    
## MAin (test)
 
precision = 2

if precision == 2 :
    RP = np.float64
else :
    RP = np.float32
      
N = 1000

X = np.zeros(N, dtype = RP)
A = np.zeros((N,N), dtype = RP)

for i in range(N):
    X[i] = rd.random()
    for j in range(N):
        A[i,j] = rd.random()
        
    A[i,i] = 50 + 50*rd.random()
        
B = np.dot(A,X)


A_d = ga.to_gpu(A)
handler = cublas.cublasCreate()
Bf = product_Ax(A_d, X, N, handler)
    

    

X0 = np.zeros(N)
Nmax = N
tol = 1e-7

Xf = GMRES_d(A_d, B, X0, tol, Nmax) 

for i in range(20) :
    print(i, X[i], Xf[i])

