import numpy as np



## Cas test

def rot(u,v):
    
    return np.cross(u,v)
    
def norm(u):
    return np.linalg.norm(u)

def calcul_S_et_n(L_T, L_P):
    
    
    N_f = L_T.shape[0]
    
    L_N = np.zeros((N_f,3))
    L_S = np.zeros(N_f)
    
    for i in range(N_f):
        P1 = L_P[L_T[i][0]]
        P2 = L_P[L_T[i][1]]
        P3 = L_P[L_T[i][2]]
        
        L_N[i] = rot(P2-P1, P3-P1)
        L_S[i] = norm(L_N[i])/2
        L_N[i] /= norm(norm(L_N[i]))
        
    return L_S, L_N
    

def calcul_G(L_T, L_P):
    
    N_f = L_T.shape[0]
    
    L_G= np.zeros((N_f,3))

    for i in range(N_f):
        P1 = L_P[L_T[i][0]]
        P2 = L_P[L_T[i][1]]
        P3 = L_P[L_T[i][2]] 
    
        L_G[i] = (P1+P2+P3)/3
    
    return L_G


def calcul_R_max(L_T, L_P, L_G):
    
    N_f = L_T.shape[0]
    
    L_R_max = np.zeros(N_f)
    
    for i in range(N_f):
        for j in range(3):
            if L_T[i][j] >= L_P.shape[0] :
                print(i,j, L_T[i][j])
    

    for i in range(N_f):
        G = L_G[i]
        
        for j in range(3):
            
            
            R = norm(G-L_P[L_T[i][j]])
            if L_R_max[i] < R :
                L_R_max[i] = R

    return L_R_max
    
def calcul_gradient_surfacique(L_T, L_P):
    
    N_f = L_T.shape[0]
    
    L_sigma = np.zeros((N_f,3,3))

    for i in range(N_f):
        P1 = L_P[L_T[i][0]]
        P2 = L_P[L_T[i][1]]
        P3 = L_P[L_T[i][2]] 
    
        A = np.dot(P3-P1,P3-P1)*(P2-P1) - np.dot(P2-P1,P3-P1)*(P3-P1)
        B = -np.dot(P2-P1,P3-P1)*(P2-P1) + np.dot(P2-P1,P2-P1)*(P3-P1)
        
        delta = np.dot(P2-P1,P2-P1)*np.dot(P3-P1,P3-P1) - np.dot(P2-P1,P3-P1)**2
    
        L_sigma[i][:,0] = -(A+B)
        L_sigma[i][:,1] = A
        L_sigma[i][:,2] = B
        
        L_sigma[i] /=delta
        
    return L_sigma
    
def calcul_phi_dphidn(L_P, L_T, L_N, R):
    
    PS = np.array([0,0,2*R])
    
    N_n = L_P.shape[0]
    N_f = L_T.shape[0]
    
    # Calcul des normales sur les points
    n = np.zeros((N_n,3))
    
    for i in range(N_f):
        for j in range(3):
             k = L_T[i][j]
             n[k] += L_N[i]
             
    for i in range(N_n):
        n[i] /= norm(n[i])
        
    # Calcul de phi et dphidn
    
    L_phi = np.zeros(N_n)
    L_dphidn = np.zeros(N_n)
    
    for i in range(N_n):
        
        L_phi[i] = 1/norm(L_P[i] - PS)
        L_dphidn[i] = np.dot(PS-L_P[i], n[i])/norm(L_P[i]-PS)**3
        
        
    return L_phi, L_dphidn
            
    
def calcul_facette_voisines(L_T, N_n):
    
    N_max = 200
    
    L_vois = np.zeros((N_n,N_max,2), dtype = np.int32)
    L_vois_N = np.zeros(N_n, dtype = np.int32)
    
    for i in range(L_T.shape[0]):
        for j in range(3):
            
            t = L_T[i,j]
            n = L_vois_N[t]
            L_vois[t,n] = [i,j]
            L_vois_N[t] += 1
            
    return L_vois, L_vois_N
            
            
    
    
    
    
    
       
def cas_test_sphere(R, N_theta, N_phi):
    
    
    d_theta = np.pi/N_theta
    d_phi = 2*np.pi/N_phi
    
    
    # Calcul des points
    
    L_P = []
    
    L_P.append([0,0,R])
    
    for i in range(1,N_theta):
        for j in range(N_phi):
            
            x = R*np.cos(j*d_phi)*np.sin(i*d_theta)
            y = R*np.sin(j*d_phi)*np.sin(i*d_theta)
            z = R*np.cos(i*d_theta)
            
            L_P.append([x,y,z])
        
    L_P.append([0,0,-R])
    
    
    # Calcul des facettes
     
    L_T = []
     
    for i in range(N_phi):
        L_T.append([0, i+1, (i+1)%N_phi + 1])   
          
    
    for i in range(N_theta-2):
        for j in range(N_phi):
            
            L_T.append([1 + j + N_phi*i, 1 + j + N_phi*(i+1), 1 + (j+1)%N_phi + N_phi*i])
            
            L_T.append([1 + (j+2)%N_phi + N_phi*i, 1 + (j+1)%N_phi + N_phi*(i+1), 1 + (j+2)%N_phi + N_phi*(i+1)])
            
    
    for i in range(N_phi):
        L_T.append([1 + (N_theta-1)*N_phi, 1 + (N_theta-1)*N_phi - i - 1, 1 + (N_theta-1)*N_phi - (i+1)%N_phi - 1])
       
       
    L_P, L_T = np.array(L_P), np.array(L_T)
       
       
    L_S, L_N = calcul_S_et_n(L_T, L_P)
    L_G = calcul_G(L_T, L_P)
    L_R_max = calcul_R_max(L_T, L_P, L_G)
    L_sigma = calcul_gradient_surfacique(L_T, L_P)
    L_phi, L_dphidn = calcul_phi_dphidn(L_P, L_T, L_N, R)
    L_vois, L_vois_N = calcul_facette_voisines(L_T, L_P.shape[0])
    
    return L_P, L_T, L_S, L_N, L_G, L_R_max, L_sigma, L_phi, L_dphidn, L_vois, L_vois_N
    
    

