import numpy as np


def C3_T3_mat(C3_stack):
    D = (1/np.sqrt(2)) * np.array([[1, 0, 1], [1, 0, -1], [0, np.sqrt(2), 0]])
    
    # First expand the dimensions of D and D.T to be broadcastable with C3_stack
    D_expanded = D[:, :, np.newaxis, np.newaxis]  # Shape: (3, 3, 1, 1)
    D_T_expanded = D.T[:, :, np.newaxis, np.newaxis]  # Shape: (3, 3, 1, 1)
    
    result = np.einsum('ij,jklm,ln->iklm', D_expanded, C3_stack, D_T_expanded)

    return result


def C3_T3(C3_Folder):
    # nrows = np.size(C3_stack,0)
    # ncols = np.size(C3_stack,1)
    
    D = (1/np.sqrt(2))*np.array([[1,0,1], [1,0,-1],[0,np.sqrt(2),0]])
    # CC1 = np.dstack((C3_stack[:,:,0].flatten(), C3_stack[:,:,1].flatten(), C3_stack[:,:,2].flatten(),
    #                  C3_stack[:,:,3].flatten(), C3_stack[:,:,4].flatten(), C3_stack[:,:,5].flatten(),
    #                  C3_stack[:,:,6].flatten(), C3_stack[:,:,7].flatten(), C3_stack[:,:,8].flatten()))
    # CC1 = CC1[0,:,:]
    # CC1 = np.reshape(CC1,(CC1.shape[0],3,3))
    # CC1 = np.matmul(np.matmul((D),CC1),D.T)
    # CC1 = np.reshape(CC1,(CC1.shape[0],9))

    return np.einsum('ij,jklm,ln->iklm', D, C3_stack, D.T)