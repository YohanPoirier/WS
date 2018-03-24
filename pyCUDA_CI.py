
import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
import numpy as np
from pycuda.compiler import SourceModule
import time
    

def init_CI(precision):
    
    code = """
    
    #include <math.h> 
    
    __device__ RP norm_d(RP *u)
    {
        RP res = 0;
        
        for (int i = 0; i < 3; ++i){
            res += u[i]*u[i];
        }
        
        return sqrt(res);
    }
    
    
    __device__ RP dot_product(RP *u, RP *v)
    {
        RP res = 0;
        for (int i = 0; i<3; ++i)
        {
            res += u[i]*v[i];
        }
    
        return res;
    }
    
    
    __device__ RP *rot_d(RP *u, RP *v)
    {
        static RP res[3];
    
        res[0] = u[1]*v[2] - u[2]*v[1];
        res[1] = u[2]*v[0] - u[0]*v[2];
        res[2] = u[0]*v[1] - u[1]*v[0];
        
        return &res[0];
    }
    
    
    __device__ void matmul(RP *X, RP *A, RP *res)
    {
        for (int i = 0; i<3; ++i)
        {
            res[i] = 0;
        }
        
        for (int i = 0; i<3; ++i)
        {
            for (int j = 0; j<3; ++j)
            {
                res[i] += X[j]*A[3*j + i];
            }
        } 
    }
    
    
    
    
    __device__ RP *vec(RP *A, RP *B)
    {
        static RP res[3];
    
        res[0] = B[0] - A[0];
        res[1] = B[1] - A[1];
        res[2] = B[2] - A[2];
        
        return &res[0];
    }
    
    
    __device__ RP S_sigma_sing_d(int t, RP *L_P, RP *M)
    {
        RP P1[3], P2[3], P3[3];
        RP i[3], j[3];
        RP *h;
        RP delta, a, b;
    
        
    
        for(int k=0; k<3; ++k)
        {
            P1[k] = L_P[ (t%3)*3 + k];
            P2[k] = L_P[ ((t+1)%3)*3 + k];
            P3[k] = L_P[ ((t+2)%3)*3 + k];
        
            i[k] = 0.5*(P2[k]-P1[k] + P3[k]-P1[k]);
            j[k] = P3[k]-P2[k];
        }
        h = rot_d(i,j);
        
        delta = dot_product(i,i)*dot_product(j,j) - pow(dot_product(i,j),2);
        
        a = 0.5*dot_product(j,j)/sqrt(delta);
        b = dot_product(i,j)/sqrt(delta);
        
        return norm_d(h)*log((a+b+sqrt(1+(a+b)*(a+b)))/(-a+b+sqrt(1+(-a+b)*(-a+b))))/norm_d(j);
    }
    
    
    
    
    __device__ void S_sigma_mu_d(RP *L_P, RP *M, RP *G, RP *n, RP *mu, RP *sigma)
    {
        RP GM[3], PkPk1[3], PkM[3], Pk1M[3];
        RP Rk, Rk1,dk, Z, Nk1, Nkt, Dk1, Dkt;
        int t;
        RP temp[3];
        
        sigma[0] = 0;
        mu[0] = 0;
        
        /* On regarde si le point etudie est un des noeuds de la facette */
        
        t = -1;
        for (int i=0; i<3; ++i)
        {
            for (int j=0; j<3; ++j)
            {
                temp[j] = L_P[3*i +j] - M[j];
            }
            if (norm_d(temp) == 0)
            {
                t = i;
            }
        }
    
        /* Cas regulier */
        if (t == -1 )
        {
            for (int i=0; i<3; ++i)
            {
                GM[i] = M[i] - G[i];
            }
    
            Z = dot_product(GM,n);
        
            for (int i=0; i<3; ++i)
            {
                for (int j=0; j<3; ++j)
                {
                    PkM[j] = M[j] - L_P[3*i + j];
                    Pk1M[j] = M[j] - L_P[3*((i+1)%3) + j];
                    PkPk1[j] = L_P[3*((i+1)%3) + j]- L_P[3*i + j];
                }
    
                Rk = norm_d(PkM);
                Rk1 = norm_d(Pk1M);
                dk = norm_d(PkPk1);
                Nk1 = Rk1 + Rk + dk;
                Dk1 = Rk1 + Rk - dk;
                Nkt = 2*dot_product(PkM, rot_d(n,PkPk1));
                Dkt = (Rk1+Rk)*(Rk1+Rk) - dk*dk + 2.*abs(Z)*(Rk1+Rk);
            
                if (abs(Z) > 1e-5)
                /* Condition a affiner */
                { 
                    mu[0] += 2.*(Z/abs(Z))*atan(Nkt/Dkt);
                }
    
                sigma[0] += Nkt*log(Nk1/Dk1)/(2.*dk) - 2.*abs(Z)*atan(Nkt/Dkt);
            
            }
        }
        /* Cas singulier */
        else
        {
            sigma[0] = S_sigma_sing_d(t, L_P, M);
        }
    }
    
    
    
            
    __device__ void S_sigma_mu_asympt_d(RP *L_P, RP *M, RP* G, RP *n, RP *mu, RP *sigma)
    {
        RP GM[3], M1M2[3], M1M3[3];
        RP *temp;
        RP A;
    
        for(int j=0; j<3; ++j)
        {
            GM[j] = M[j]-G[j];
            M1M2[j] = L_P[3 + j] - L_P[j];
            M1M3[j] = L_P[6 + j] - L_P[j];
        }
    
        temp = rot_d(M1M2,M1M3);
        A = norm_d(temp)/2.;
        
        sigma[0] = A/norm_d(GM);
        mu[0] = dot_product(GM,n)*A/pow(norm_d(GM),3);
        /* Erreur dans la these de Lucas */
    }
        
    
    /* Teste si le point correspond a une singularite */
    
    __device__ int test_sing_d (RP *L_P, RP *M)
    {
        RP PM[3];
        int sing = 0;
    
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                PM[j] = M[j]-L_P[i*3 + j];
            }
            if (norm_d(PM) == 0)
            {
                sing = 1;
            }
        }    
        return sing;    
    }
    
    
    /*Teste si le point est loin(approximation asymptotique)*/
    
    __device__ int test_asymp_d(RP *G, RP *M, RP R_max, int N_seuil)
    {
        RP GM[3];
    
        for(int j = 0; j<3 ; ++j)
        {
            GM[j] = M[j] - G[j];
        }
        return (norm_d(GM) > N_seuil*R_max);
    }
    
    
    
    __device__ void I_sigma_mu_d(RP *L_P, RP *M, RP *n, RP *mu, RP *sigma, RP *temp, int i_n, int i_f, int N_n_max)
    {
    
        RP A[3], B[3], AB[3], AM[3], BM[3];
        RP *res;
        RP K2;
        RP a_m, b_m;
        RP AB_n, dp_AB_AM;
        
 
        
        
        for(int j = 0; j<3 ; ++j)
        {
            mu[j] = 0;
            sigma[j] = 0;
        }
    
        for(int i = 0; i<3 ; ++i)
        {
            for(int j = 0; j<3 ; ++j)
            {
                A[j] = L_P[3*i + j];
                B[j] = L_P[3*((i+1)%3) + j];
                
                AB[j] = B[j] - A[j];
                AM[j] = M[j] - A[j];
                BM[j] = M[j] - B[j];
            }
            
        
            dp_AB_AM = dot_product(AB,AM);
            
            AB_n = norm_d(AB);
            

            
            K2 = dot_product(AM,AM) - (dp_AB_AM/AB_n)*(dp_AB_AM/AB_n);
            
   
    
            if (K2 > 1e-5)
            {
            
                
                a_m = asinh(-dp_AB_AM /(AB_n*sqrt(K2)));
                b_m = asinh( (AB_n - dp_AB_AM/AB_n )/sqrt(K2) );
                
                
                
                /*b_m = 1.;*/
                /*a_m = 1.;*/
                
                res = rot_d(AM,AB);
                
                
                for(int j = 0; j<3 ; ++j)
                {
                    mu[j] += -res[j]*(b_m-a_m)/AB_n ;
                    /* Erreur de signe dans la these de Lucas */
                }
                
                
                
                res = rot_d(n,AB);
                
                
                
                
                
                for(int j = 0; j<3 ; ++j)
                {

                    sigma[j] += K2/(2.*AB_n)*res[j]*(b_m-a_m+(sinh(2.*b_m)-sinh(2.*a_m))/2.);
  
                }

            
              
                

                
                
         
                
                
            }
            else
            {
                if (norm_d(AM) < 1e-5 || norm_d(BM) < 1e-5 )
                /* Condition a affiner */
                { 
                    res = rot_d(n,AB);
                    
                    for(int j = 0; j<3 ; ++j)
                    {
                        sigma[j] += 0.5*AB_n*res[j];
                    }
                }
                else
                {
                    res = rot_d(n,AB);
                    
                    for(int j = 0; j<3 ; ++j)
                    {
                        sigma[j] += abs(norm_d(AM) - 0.5 * dp_AB_AM/norm_d(AM))*res[j];
                        
                    }

                    /* Erreur dans la these de Lucas (signe "-")*/
        
                     
                }
            }
        }
    }
    
    
    __device__ void coeff_inf_d(RP *L_P, RP *n, RP *G, RP *GS, RP R_max, RP *M, RP *CS, RP *CD, RP *temp, int i_n, int i_f, int N_n_max)
    {

        RP GM[3];
        RP Ssigma, Smu;
        RP Isigma[3], Imu[3];
        RP temp1[3] , temp2[3], temp3[3];
        int N_seuil = 8;
        
        I_sigma_mu_d(L_P, M, n, Imu, Isigma, temp, i_n, i_f, N_n_max);
        
        if (test_asymp_d (G,M,R_max,N_seuil))
        {
            S_sigma_mu_asympt_d(L_P, M, G, n, &Smu, &Ssigma);
        }
        else
        {
            S_sigma_mu_d(L_P, M, G, n, &Smu, &Ssigma);
        }
        
        
        
        
    
        for(int j = 0; j<3 ; ++j)
        {
            GM[j] = M[j] - G[j];
        } 
        
        
        matmul(GM, GS, temp1);
        matmul(Imu, GS, temp2);
        matmul(Isigma, GS, temp3);
        
        for(int j = 0; j<3 ; ++j)
        {
            CD[j] = (1./3. + temp1[j])*Smu - temp2[j];
            CS[j] = (1./3. + temp1[j])*Ssigma - temp3[j];
        }
        

        
    }
    
    
    __global__ void mat_CI_kernel(RP *L_X, int *L_T, RP *L_n, RP *L_GS, RP *L_G, RP *L_R_max, RP *A_CS, RP *A_CD, int i_n_i, int i_n_f, int i_f_i, int i_f_f, int N_n_max, int N_sym, double prof, RP *temp)
    {
        RP CS[3]={0,0,0}, CD[3]={0,0,0};
        RP L_P[9], M[3], GS[9], n[3], G[3];
        int i_n;
        int i;

        i_n = blockIdx.x* blockDim.x + threadIdx.x + i_n_i;
        
        if (i_n <= i_n_f)
        {
         
            /*
            for (int k=0; k<N_n; ++k)
            {
                A_CS[i_n*N_n + k] = 0;
                A_CD[i_n*N_n + k] = 0;
            }
            */
            
            for (int k=0; k<2; ++k)
            {
                M[k] = L_X[i_n*3 + k];
            }
            
            
            for (int i_sym=0; i_sym<=N_sym; ++i_sym)
            {
                M[2] = (1-2*i_sym)*L_X[i_n*3 + 2] - 2.*i_sym*prof;
                
                
    
                for (int i_f=i_f_i; i_f<=i_f_f; ++i_f)
                {
                
                    for (int j=0; j<3; ++j)
                    {
                        i = L_T[3*i_f + j];
                        
                        for (int k=0; k<3; ++k)
                        {
                            L_P[3*j + k] = L_X[3*i + k];
                        }
                    }
                    
                    for (int k=0; k<9; ++k)
                    {
                        GS[k] = L_GS[9*i_f + k];
                    }
                    
                    for (int k=0; k<3; ++k)
                    {
                        n[k] = L_n[3*i_f + k];
                        G[k] = L_G[3*i_f + k];
                    }
                    

                    coeff_inf_d(L_P, n, G, GS, L_R_max[i_f], M, CS, CD, temp, i_n, i_f, N_n_max);
                

     
                    
                    temp[i_n] += CD[0];
                    
                  
                    
                    for (int k=0; k<3; ++k)
                    {
                        i = L_T[3*i_f + k];
                        
                        A_CS[i_n*N_n_max] += CS[k];
                        A_CD[i_n*N_n_max] += CD[k];
                        
                    }
                    
                }
            }
            
            

            
        }
    }
    

    __global__ void angle_solide_kernel(RP *A_CD, int  N_n, int N_n_max)
    {
        int i_n;
        
        
        i_n = blockIdx.x* blockDim.x + threadIdx.x;
        
        A_CD[i_n*N_n_max + i_n] = 0;
        
        for (int k=0; k<i_n; ++k)
        {
            A_CD[i_n*N_n_max + i_n] += -A_CD[i_n*N_n_max + k];
        }
        
        for (int k=i_n+1; k<N_n; ++k)
        {
            A_CD[i_n*N_n_max + i_n] += -A_CD[i_n*N_n_max + k];
        }

    }
        
        
    
    
    __global__ void calcul_GPU(RP *dest, RP *L_P)
    {
    const int i = threadIdx.x;
    dest[i] = norm_d(&L_P[3*i]);
    }
    
    """
    
    
    
    # Reglage de la precision
    if precision == 2 :
        RP = "\n    typedef double RP;"
    else :
        RP = "\n    typedef float RP;"
        
    code = RP + code
    
    mod = SourceModule(code)
    
    return mod.get_function("mat_CI_kernel"), mod.get_function("angle_solide_kernel")
    
    
def creation_matrices_GPU(N_n_max, precision):
    
    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32

    
    A_CS_d = gpuarray.zeros((N_n_max,N_n_max), RP)
    A_CD_d = gpuarray.zeros((N_n_max,N_n_max), RP)
    
    return A_CS_d, A_CD_d
    
    
    
    
def export_mesh_vers_GPU(L_P, L_T, L_N, L_G, L_R_max, L_ds, precision):
    
    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32
        
    N_n = L_P.shape[0]
    N_f = L_T.shape[0]


    L_X_d = gpuarray.to_gpu(np.array(L_P, dtype = RP).reshape(3*N_n))
    L_T_d = gpuarray.to_gpu(np.array(L_T, dtype = np.int32).reshape(3*N_f))
    L_N_d = gpuarray.to_gpu(np.array(L_N, dtype = RP).reshape(3*N_f))
    L_ds_d = gpuarray.to_gpu(np.array(L_ds, dtype = RP).reshape(3*3*N_f))
    L_G_d = gpuarray.to_gpu(np.array(L_G, dtype = RP).reshape(3*N_f))
    L_R_max_d = gpuarray.to_gpu(np.array(L_R_max, dtype = RP))
    # L_double_d = gpuarray.to_gpu(np.array(L_double, dtype = np.int32).reshape(3*N_n))
    # L_double_N_d = gpuarray.to_gpu(np.array(L_double_N, dtype = np.int32))
        

    return L_X_d, L_T_d, L_N_d, L_ds_d, L_G_d, L_R_max_d
    
    
    
def calcul_angle_solide(A_CD_d, N_n, N_n_max, angle_solide_kernel):
    
    BLOCK_SIZE = 64
    
    angle_solide_kernel(
            A_CD_d,
            np.int32(N_n),
            np.int32(N_n_max),
            block=(BLOCK_SIZE,1,1), 
            grid=((N_n-1)//BLOCK_SIZE+1,1)
            )
            
    return A_CD_d
    
    
def calcul_matrice_CI(A_CS_d, A_CD_d, L_X_d, L_T_d, L_N_d, L_G_d, L_R_max_d, L_ds_d, N_sym, prof, mat_CI_kernel, bornes, N_n_max):



    BLOCK_SIZE = 64

    i_n_i = int(bornes[0])
    i_n_f = int(bornes[1])
    i_f_i = int(bornes[2])
    i_f_f = int(bornes[3])
    
    
    print(i_n_i, i_n_f, N_n_max)
    
    
    precision = 2
    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32


    N_sym = 0

    mat_CI_kernel(
            L_X_d, L_T_d, L_N_d, L_ds_d, L_G_d, L_R_max_d,
            A_CS_d, A_CD_d,
            np.int32(i_n_i), np.int32(i_n_f),
            np.int32(i_f_i), np.int32(i_f_f),
            np.int32(N_n_max),
            np.int32(N_sym), np.float64(prof),
            temp_d,
            block=(BLOCK_SIZE,1,1), 
            grid=((i_n_f-i_n_i)//BLOCK_SIZE+1,1)
            )
            
            
    
            
    drv.Context.synchronize()
    
    

    # print("chien")
    # 
    # temp = temp_d.get()
    # 
    # print(temp[0])
    # for i in range(i_n_f):
    #     if temp[i] == 12 :
    #         print(i,temp[i])
    #     
    # input()
    
    print("chien")
    


    return A_CS_d, A_CD_d



