import numpy as np


def ecriture(L_P, L_ds, L_T, L_G, L_N, L_Rmax):

    
    
    N_n = L_P.shape[0]
    N_f = L_ds.shape[0]
    
    file = open("mesh.txt", "w")
    
    file.write("{}\n".format(N_n))
    for i in range(N_n):
        file.write("{} {} {}\n".format(L_P[i,0], L_P[i,1], L_P[i,2]))
    
    
    file.write("{}\n".format(N_f))
    
    L_ds = L_ds.reshape([N_f,9])
    
    for i in range(N_f):
        file.write("{} {} {} {} {} {} {} {} {}\n".format(L_ds[i,0], L_ds[i,1], L_ds[i,2],L_ds[i,3], L_ds[i,4], L_ds[i,5],L_ds[i,6], L_ds[i,7], L_ds[i,8]))

    for i in range(N_f):
        file.write("{} {} {}\n".format(L_T[i,0], L_T[i,1], L_T[i,2])) 
        
    for i in range(N_f):
        file.write("{} {} {}\n".format(L_G[i,0], L_G[i,1], L_G[i,2]))
        
    for i in range(N_f):
        file.write("{} {} {}\n".format(L_N[i,0], L_N[i,1], L_N[i,2]))
        
    for i in range(N_f):
        file.write("{}\n".format(L_Rmax[i]))
        
    file.close()