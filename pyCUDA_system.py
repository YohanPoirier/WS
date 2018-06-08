import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
from pycuda.compiler import SourceModule

import numpy as np



def init_systeme(precision):
    
    code = """
 
    __global__ void systeme_kernel(RP *A, RP *B, RP *X_connu, RP *A_CS, RP *A_CD, int *ind, int N_n, int N_n_max)
    {
        int i_n ;
        RP diag;
        i_n = blockIdx.x*blockDim.x + threadIdx.x;
        
        if (i_n < N_n)
        {
            /*Surface libre*/
            for (int j=ind[2]; j<=ind[3]; ++j)
            {
                A[i_n*N_n_max + j] = A_CS[i_n*N_n_max + j];
                B[i_n] += A_CD[i_n*N_n_max + j]*X_connu[j];
            }
            
            /*Corps*/
            for (int j=ind[4]; j<=ind[5]; ++j)
            {
                A[i_n*N_n_max + j] = -A_CD[i_n*N_n_max + j];
                B[i_n] += -A_CS[i_n*N_n_max + j]*X_connu[j];
            }  
            
            diag = (RP)1./A[i_n*N_n_max + i_n];
            B[i_n] *= diag;
           

            for (int j = 0; j<N_n; ++j)
            {
                A[i_n*N_n_max + j] *= diag;
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
    
    return mod.get_function("systeme_kernel")
    
def construction_systeme(A_d, A_CD_d, A_CS_d, X_connu, ind, N_n_max, cuve_ferme, systeme_kernel, precision):
    
    N_n = X_connu.shape[0]
    
    BLOCK_SIZE = 64

    if precision == 2 :
        RP = np.float64
    else :
        RP = np.float32

    B_d = gpuarray.zeros(N_n, RP)

    ind_d = gpuarray.to_gpu(ind)



    Xconnu_d = gpuarray.to_gpu(np.array(X_connu, RP))
    

    systeme_kernel(
            A_d, B_d,
            Xconnu_d, 
            A_CS_d, A_CD_d,
            ind_d,
            np.int32(N_n), np.int32(N_n_max),
            block=(BLOCK_SIZE,1,1), 
            grid=((N_n-1)//BLOCK_SIZE+1,1)
            )
            
    B = B_d.get()

    return B