import numpy as np


def C3_T3_mat(C3_stack):
    # Ensure C3_stack has the correct shape
    assert C3_stack.ndim == 4 and C3_stack.shape[0] == 3 and C3_stack.shape[1] == 3, \
        "C3_stack must have shape (3, 3, n, m)"

    # Define the transformation matrix D
    D = (1 / np.sqrt(2)) * np.array([[1, 0, 1], 
                                      [1, 0, -1], 
                                      [0, np.sqrt(2), 0]])

    # Perform the first tensor contraction with D
    # This will contract over the first axis of C3_stack and the second axis of D
    temp_result = np.tensordot(D, C3_stack, axes=([1], [0]))  # Shape: (3, n, m)

    # Now we need to perform the second contraction with D.T
    # Contract temp_result with D.T on the second axis
    result = np.tensordot(temp_result, D.T, axes=([1], [0]))  # Shape: (3, 3, n, m)

    return np.transpose(result, (0, 3, 1, 2))


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

    return 0