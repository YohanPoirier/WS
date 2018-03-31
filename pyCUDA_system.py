import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
from pycuda.compiler import SourceModule

import numpy as np



def init_systeme(precision):
    
    code = """
 
    __global__ void systeme_kernel(float *A, float *B, float *X_connu, float *A_CS, float *A_CD, int *ind, int N_n, int N_n_max)
    {
        int i_n ;
        RP diag;
        i_n = blockIdx.x*blockDim.x + threadIdx.x;
        
        if (i_n < N_n)
        {
            /*Surface libre*/
            for (int j=ind[2]; j<ind[3]; ++j)
            {
                A[i_n*N_n + j] = A_CS[i_n*N_n_max + j];
                B[i_n] += A_CD[i_n*N_n_max + j]*X_connu[j];
            }
            
            /*Corps*/
            for (int j=ind[4]; j<ind[5]; ++j)
            {
                A[i_n*N_n + j] = -A_CD[i_n*N_n_max + j];
                B[i_n] += -A_CS[i_n*N_n_max + j]*X_connu[j];
            }  
        }
        
        
        diag = (RP)1./A[i_n*N_n + i_n];
        B[i_n] *= diag;
        
        for (int j = 0; j<N_n; ++j)
        {
            A[i_n*N_n + j] *= diag;
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
    
    return mod.get_function("systeme_kernel")
    
def construction_systeme(A_CD_d, A_CS_d, X_connu, ind, N_n_max, cuve_ferme, systeme_kernel, precision):
    
    N = Xconnu.shape[0]
    
    BLOCK_SIZE = 64

    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32
        
    A_d = gpuarray.zeros((N,N), RP)
    B_d = gpuarray.zeros(N, RP)

    ind_d = gpuarray.to_gpu(ind)

    Xconnu_d = gpuarray.to_gpu(np.array(X_connu, RP))

    systeme_kernel(
            A_d, B_d,
            Xconnu_d, 
            A_CS_d, A_CD_d,
            ind_d,
            np.int32(N), np.int32(N_n_max),
            block=(BLOCK_SIZE,1,1), 
            grid=((N-1)//BLOCK_SIZE+1,1)
            )
            
    B = B_d.get()
    
    return A_d, B