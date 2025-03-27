import numpy as np


def C3_T3_mat(C3_stack):
    # Ensure C3_stack has the correct shape
    assert C3_stack.ndim == 4 and C3_stack.shape[0] == 3 and C3_stack.shape[1] == 3, \
        "C3_stack must have shape (3, 3, n, m)"

    nrows = np.size(C3_stack,2)
    ncols = np.size(C3_stack,3)
    C3_stack = C3_stack.transpose(2, 3, 0, 1).reshape(nrows, ncols, 9)
    

    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    TT1 = np.dstack((C3_stack[:,:,0].flatten(), C3_stack[:,:,1].flatten(), C3_stack[:,:,2].flatten(),
                     C3_stack[:,:,3].flatten(), C3_stack[:,:,4].flatten(), C3_stack[:,:,5].flatten(),
                     C3_stack[:,:,6].flatten(), C3_stack[:,:,7].flatten(), C3_stack[:,:,8].flatten()))
    TT1 = TT1[0,:,:]
    TT1 = np.reshape(TT1,(TT1.shape[0],3,3))
    TT1 = np.matmul(np.matmul(D,TT1),D.T)
    TT1 = np.reshape(TT1,(TT1.shape[0],9))
    TT1 = np.reshape(TT1,(nrows,ncols,9)) 
    
    return TT1.reshape(nrows,ncols, 3, 3).transpose(2, 3, 0, 1)

def T3_C3_mat(T3_stack):
    # Ensure T3_stack has the correct shape
    assert T3_stack.ndim == 4 and T3_stack.shape[0] == 3 and T3_stack.shape[1] == 3, \
        "T3_stack must have shape (3, 3, n, m)"

    nrows = np.size(T3_stack,2)
    ncols = np.size(T3_stack,3)
    T3_stack = T3_stack.transpose(2, 3, 0, 1).reshape(nrows, ncols, 9)

    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    TT1 = np.dstack((T3_stack[:,:,0].flatten(), T3_stack[:,:,1].flatten(), T3_stack[:,:,2].flatten(),
                     T3_stack[:,:,3].flatten(), T3_stack[:,:,4].flatten(), T3_stack[:,:,5].flatten(),
                     T3_stack[:,:,6].flatten(), T3_stack[:,:,7].flatten(), T3_stack[:,:,8].flatten()))
    TT1 = TT1[0,:,:]
    TT1 = np.reshape(TT1,(TT1.shape[0],3,3))
    TT1 = np.matmul(np.matmul((D.T),TT1),D)
    TT1 = np.reshape(TT1,(TT1.shape[0],9))
    TT1 = np.reshape(TT1,(nrows,ncols,9)) 
    
    return TT1.reshape(nrows,ncols, 3, 3).transpose(2, 3, 0, 1)
