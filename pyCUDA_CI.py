
import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
import numpy as np
from pycuda.compiler import SourceModule
import time
    

def init_CI(precision):
    
    code = """
    
    #include <math.h> 
    
    __constant__ double inv_3 = (double)1./(double)3. ;

    __device__ RP norm_d(RP *u)
    {
        RP res = 0;
        
        #pragma unroll
        for (int i = 0; i < 3; ++i){
            res += u[i]*u[i];
        }
        
        return sqrt(res);
    }
    
    
    __device__ RP dot_product(RP *u, RP *v)
    {
        RP res = 0;
        
        #pragma unroll
        for (int i = 0; i<3; ++i)
        {
            res += u[i]*v[i];
        }
    
        return res;
    }
    
    
    __device__ void rot_d(RP *res, RP *u, RP *v)
    {
        res[0] = u[1]*v[2] - u[2]*v[1];
        res[1] = u[2]*v[0] - u[0]*v[2];
        res[2] = u[0]*v[1] - u[1]*v[0];

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
    
    
    
    
    __device__ RP S_sigma_sing_d(int t, RP *L_P, RP *M)
    {
        RP P1[3], P2[3], P3[3];
        RP i[3], j[3];
        RP h[3];
        RP delta, a, b;
    
        for(int k=0; k<3; ++k)
        {
            P1[k] = L_P[ (t%3)*3 + k];
            P2[k] = L_P[ ((t+1)%3)*3 + k];
            P3[k] = L_P[ ((t+2)%3)*3 + k];
        
            i[k] = 0.5*(P2[k]-P1[k] + P3[k]-P1[k]);
            j[k] = P3[k]-P2[k];
        }
        rot_d(h, i,j);
        
        delta = dot_product(i,i)*dot_product(j,j) - pow(dot_product(i,j),2);
        
        a = 0.5*dot_product(j,j)/sqrt(delta);
        b = dot_product(i,j)/sqrt(delta);
        
        return norm_d(h)*log((a+b+sqrt(1+(a+b)*(a+b)))/(-a+b+sqrt(1+(-a+b)*(-a+b))))/norm_d(j);
    }
    
    __device__ void S_sigma_mu2_d(RP *L_P, RP *M, RP *G, RP *n, RP *mu, RP *sigma)
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
            if (norm_d(temp) < 1e-7)
            /*Condition a affiner*/
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
                rot_d(temp,n,PkPk1);
                Nkt = 2*dot_product(PkM,temp);
                Dkt = (Rk1+Rk)*(Rk1+Rk) - dk*dk + 2.*abs(Z)*(Rk1+Rk);
            
                if (abs(Z) > 1e-7)
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
            if (norm_d(temp) < 1e-7)
            /*Condition a affiner*/
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
                rot_d(temp,n,PkPk1);
                Nkt = 2*dot_product(PkM,temp);
                Dkt = (Rk1+Rk)*(Rk1+Rk) - dk*dk + 2.*abs(Z)*(Rk1+Rk);
            
                if (abs(Z) > 1e-7)
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
    

            
    __device__ void coeff_asympt_d(RP *L_P, RP* GM, RP *n, RP *CD, RP *CS)
    {
        RP M1M2[3], M1M3[3];
        RP temp[3];
        RP A;
        RP Css, Cdd;
        RP inv3;
        
        inv3 = (double)1./(double)3. ;

        #pragma unroll
        for(int j=0; j<3; ++j)
        {
            M1M2[j] = L_P[3 + j] - L_P[j];
            M1M3[j] = L_P[6 + j] - L_P[j];
        }
    
        rot_d(temp,M1M2,M1M3);
        
        A = norm_d(temp)/2.;
        
        Css = inv3*A/norm_d(GM);
        Cdd = inv3*dot_product(GM,n)*A/pow(norm_d(GM),3); /* Erreur dans la these de Lucas */
        

        for(int j=0; j<3; ++j)
        {
            CD[j] += Cdd;
            CS[j] += Css;
        }    
    }
        
    
    /* Teste si le point correspond a une singularite */
    
    __device__ int test_sing_d (RP *L_P, RP *M)
    {
        RP PM[3];

        #pragma unroll
        for(int i=0; i<3; ++i)
        {
            for(int j=0; j<3; ++j)
            {
                PM[j] = M[j]-L_P[i*3 + j];
            }
            if (norm_d(PM) == 0)
            {
                return 1;
            }
        }
        return 0;     
    }
    
    
    /*Teste si le point est loin(approximation asymptotique)*/
    
    __device__ int test_asymp_d(RP *G, RP *M, RP Cr_max)
    {
        RP Cr = 0.;
        
        #pragma unroll
        for(int j = 0; j<3 ; ++j)
        {
            Cr += pow(M[j] - G[j],2);
        }
        return (Cr > Cr_max);
    }
    
    
    __device__ void I_sigma_mu_d(RP *L_P, RP *M, RP *n, RP *mu, RP *sigma)
    {
    
        RP A[3], B[3], AB[3], AM[3], BM[3];
        RP res[3];
        RP K2;
        RP a_m, b_m;
        RP AB_n, AM_n, dp_AB_AM;
        RP coeff;
        
        #pragma unroll
        for(int j = 0; j<3 ; ++j)
        {
            mu[j] = 0;
            sigma[j] = 0;
        }
    
        #pragma unroll
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

            AB_n = norm_d(AB);
            AM_n = norm_d(AM);
            

            if (AM_n < 1e-7 || norm_d(BM) < 1e-7 )
                /* Condition a affiner */
                { 
                    rot_d(res,n,AB);
                    
                    #pragma unroll
                    for(int j = 0; j<3 ; ++j)
                    {
                        sigma[j] += 0.5*AB_n*res[j];
                    }
                }
            else
            {
                dp_AB_AM = dot_product(AB,AM);
                
            
                if (abs(abs(dp_AB_AM) - AB_n*AM_n) < 1e-7)
                {
                    rot_d(res,n,AB);
                    
                    coeff = abs(AM_n - 0.5 * dp_AB_AM/AM_n);
                    
                    #pragma unroll
                    for(int j = 0; j<3 ; ++j)
                    {
                        sigma[j] += coeff*res[j];
                    }
    
                    /* Erreur dans la these de Lucas (signe "-")*/
                }
                else
                {
                    K2 = dot_product(AM,AM) - (dp_AB_AM/AB_n)*(dp_AB_AM/AB_n);
    
                    a_m = asinh(-dp_AB_AM /(AB_n*sqrt(K2)));
                    b_m = asinh( (AB_n - dp_AB_AM/AB_n )/sqrt(K2) );
    
                    rot_d(res,AM,AB);
                    
                    coeff = -(b_m-a_m)/AB_n ;
    
                    #pragma unroll
                    for(int j = 0; j<3 ; ++j)
                    {
                        mu[j] += res[j]*coeff ;
                        /* Erreur de signe dans la these de Lucas */
                        
                    }
                        
                    rot_d(res,n,AB);
                    
                    coeff = K2/(2.*AB_n)*(b_m-a_m+(sinh(2.*b_m)-sinh(2.*a_m))/2.);
                    
                    #pragma unroll
                    for(int j = 0; j<3 ; ++j)
                    {
                        sigma[j] += coeff*res[j];
                    }
                }
            }
        }
    }
    

    __device__ void coeff_inf_d(RP *L_P, RP *n, RP *G, RP *GS, RP Cr_max, RP *M, RP *CS, RP *CD)
    {

        RP GM[3];
        RP Ssigma, Smu;
        RP Isigma[3], Imu[3];
        RP temp1[3] , temp2[3], temp3[3];
        RP inv3 = (double)1./(double)3. ;
        
        #pragma unroll
        for(int j = 0; j<3 ; ++j)
        {
            GM[j] = M[j] - G[j];
        }
            
        if (dot_product(GM,GM) > Cr_max)
        {
            coeff_asympt_d(L_P, GM, n, CD, CS);
            
        }
        else
        {
            I_sigma_mu_d(L_P, M, n, Imu, Isigma);
            S_sigma_mu_d(L_P, M, G, n, &Smu, &Ssigma);
            
            matmul(GM, GS, temp1);
            matmul(Imu, GS, temp2);
            matmul(Isigma, GS, temp3);
            
            #pragma unroll
            for(int j = 0; j<3 ; ++j)
            {
                CD[j] += (inv3 + temp1[j])*Smu - temp2[j];
                CS[j] += (inv3 + temp1[j])*Ssigma - temp3[j];

            }
        }
    }
    
    
    __global__ void mat_CI_kernel(RP *L_X, int *L_T, RP *L_n, RP *L_GS, RP *L_G, RP *L_Cr_max, RP *A_CD, RP *A_CS, int N_n, int N_n_max, int N_f, int *liste_n, int *liste_f, int N_sym, RP prof, RP * temp)
    {
        
        RP L_P[9], M[3];
        int i_n, i_f;
        int i_dx, i;

        i_dx = blockIdx.x* blockDim.x + threadIdx.x;
        
        if (i_dx < N_n)
        {
            i_n = liste_n[i_dx];
            for (int k=0; k<2; ++k)
            {
                M[k] = L_X[i_n*3 + k];
            }

            for (int j_dx = 0; j_dx < N_f; ++j_dx)
            {
                i_f = liste_f[j_dx];
                
                for (int j=0; j<3; ++j)
                {
                    i = L_T[3*i_f + j];
                    
                    for (int k=0; k<3; ++k)
                    {
                        L_P[3*j + k] = L_X[3*i + k];
                    }
                }
                
                RP CS[3] = {0};
                RP CD[3] = {0};
                    
                for (int i_sym=0; i_sym<=N_sym; ++i_sym)

                {
                    M[2] = (1-2*i_sym)*L_X[i_n*3 + 2] - 2.*i_sym*prof;
                    coeff_inf_d(L_P, &L_n[3*i_f], &L_G[3*i_f], &L_GS[9*i_f], L_Cr_max[i_f], M, CS, CD);
                }
                
                for (int k=0; k<3; ++k)
                {
                    i = L_T[3*i_f + k];
              
                
                    A_CS[i_n*N_n_max+i] += CS[k];
                    A_CD[i_n*N_n_max+i] += CD[k];
                } 
            } 
        }
    }
    


    __global__ void angle_solide_kernel(RP *A_CD, int N_n, int N_n_max)
    {
        int i_n;
        
        i_n = blockIdx.x* blockDim.x + threadIdx.x;
        
        if (i_n < N_n)
        
        {
            A_CD[i_n*N_n + i_n] = 0;
            
            for (int k=0; k<i_n; ++k)
            {
                A_CD[i_n*N_n_max + i_n] += -A_CD[i_n*N_n_max + k];
            }
            
            for (int k=i_n+1; k<N_n; ++k)
            {
                A_CD[i_n*N_n_max + i_n] += -A_CD[i_n*N_n_max + k];
            }
        }

    }
    
    __global__ void init_CI_kernel(RP *A_CD, RP *A_CS, int N_n, int N_n_max)
    {

        int i_n;

        i_n = blockIdx.x* blockDim.x + threadIdx.x;
        
        if (i_n < N_n)
        {
            for (int k=0; k<N_n; ++k)
            {
                A_CS[i_n*N_n_max + k] = 0;
                A_CD[i_n*N_n_max + k] = 0;
            }
        }
    }
        
    """
    
    
    
    # Reglage de la precision
    if precision == 2 :
        RP = "\n    typedef double RP;"
    else :
        RP = "\n    typedef float RP;"
        
    code = RP + code
    
    mod = SourceModule(code)
    
    return mod.get_function("init_CI_kernel"), mod.get_function("mat_CI_kernel"), mod.get_function("angle_solide_kernel")
    
    
def creation_matrices_GPU(N_n_max, precision):
    
    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32

    A_CD_d = gpuarray.zeros((N_n_max,N_n_max), RP)
    A_CS_d = gpuarray.zeros((N_n_max,N_n_max), RP)
    A_d = gpuarray.zeros((N_n_max,N_n_max), RP)
    
    return A_CD_d, A_CS_d, A_d
    
    
def calcul_angle_solide(A_CD_d, N_n, N_n_max, angle_solide_kernel):
    
    BLOCK_SIZE = 64

    angle_solide_kernel(
            A_CD_d,
            np.int32(N_n), np.int32(N_n_max),
            block=(BLOCK_SIZE,1,1), 
            grid=((N_n-1)//BLOCK_SIZE+1,1)
            )
            
    drv.Context.synchronize()

    
def calcul_matrice_CI(A_CD_d, A_CS_d, L_X, L_T, L_N, L_G, L_Cr_max, L_ds, N_sym, prof, mat_CI_kernel, precision, N_n_max):


    BLOCK_SIZE = 256
    

    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32

    N_n = L_X.shape[0]
    N_f = L_T.shape[0]
    
    
    L_X_d = gpuarray.to_gpu(np.array(L_X, dtype = RP).reshape(3*N_n))
    L_T_d = gpuarray.to_gpu(np.array(L_T, dtype = np.int32).reshape(3*N_f))
    L_N_d = gpuarray.to_gpu(np.array(L_N, dtype = RP).reshape(3*N_f))
    L_ds_d = gpuarray.to_gpu(np.array(L_ds, dtype = RP).reshape(3*3*N_f))
    L_G_d = gpuarray.to_gpu(np.array(L_G, dtype = RP).reshape(3*N_f))
    L_Cr_max_d = gpuarray.to_gpu(np.array(L_Cr_max, dtype = RP))
    
    liste_n_d = gpuarray.to_gpu(np.arange(N_n, dtype = np.int32))
    liste_f_d = gpuarray.to_gpu(np.arange(N_f, dtype = np.int32))
    # N_block = 10000

   ##   N_iter = (N_f-2)//N_block + 1
    # 

   ##   for k in range(N_iter):
    #     i_f_i = k*N_block
    #     i_f_f = min( (k+1)*N_block, N_f-1)
    #     


    
    mat_CI_kernel(
            L_X_d, L_T_d, L_N_d, L_ds_d, L_G_d, L_Cr_max_d,
            A_CD_d, A_CS_d, 
            np.int32(N_n), np.int32(N_n_max), np.int32(N_f),
            liste_n_d, liste_f_d,
            np.int32(N_sym), RP(prof),

            block=(BLOCK_SIZE,1,1), 
            grid=((N_n-1)//BLOCK_SIZE+1,1)
            )

    drv.Context.synchronize()
    
    

        

def init_matrice_CI(init_CI_kernel, A_CD_d, A_CS_d, N_n, N_n_max):

    BLOCK_SIZE = 64

    init_CI_kernel(
            A_CD_d, A_CS_d,
            np.int32(N_n), np.int32(N_n_max),
            block=(BLOCK_SIZE,1,1), 
            grid=((N_n-1)//BLOCK_SIZE+1,1)
            )

    drv.Context.synchronize()




