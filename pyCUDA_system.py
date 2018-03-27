import pycuda.autoinit
import pycuda.driver as drv
from pycuda import gpuarray
from pycuda.compiler import SourceModule

import numpy as np



def init_systeme():
    
    mod = SourceModule("""
 
    __global__ void systeme_kernel(float *A, float *B, float* Xconnu, int* type, float* A_CS, float *A_CD, int N)
    {
        int idx ;
        
        idx = blockIdx.x*blockDim.x + threadIdx.x;
        
        if (idx < N)
        {
            for (int j=0; j<N; j++)
            {
                if (type[j] == 1)
                {
                    A[idx*N + j] = -A_CD[idx*N + j];
                    B[idx] += -A_CS[idx*N + j]*Xconnu[j];
                }
                else
                {
                    A[idx*N + j] = A_CS[idx*N + j];
                    B[idx] += A_CD[idx*N + j]*Xconnu[j];
                }
            }
        }
    }
    
    """)
    

    
    return mod.get_function("systeme_kernel")
    
def construction_systeme(A_CD, A_CS, Xconnu, L_type, systeme_kernel):
    
    N = Xconnu.shape[0]
    
    A_d = gpuarray.zeros((N,N), np.float32)
    B_d = gpuarray.zeros(N, np.float32)
    diago_A_d = gpuarray.zeros(N, np.float32)
    
    Xconnu_d = gpuarray.to_gpu(np.array(Xconnu, dtype = np.float32))
    L_type_d = gpuarray.to_gpu(np.array(L_type, dtype = np.int32))
    
    
    BLOCK_SIZE = 32

    systeme_kernel(
            A_d, B_d, Xconnu_d, L_type_d, A_CS, A_CD, np.int32(N),
            block=(BLOCK_SIZE,1,1), 
            grid=((N-1)//BLOCK_SIZE+1,1)
            )
            
    return A_d, B_d