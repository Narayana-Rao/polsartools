import numpy as np


def C3_T3_mat(C3_stack):
    # Ensure C3_stack has the correct shape
    assert C3_stack.ndim == 4 and C3_stack.shape[0] == 3 and C3_stack.shape[1] == 3, \
        "C3_stack must have shape (3, 3, n, m)"

    # Define the transformation matrix D
    D = (1 / np.sqrt(2)) * np.array([[1, 0, 1], 
                                      [1, 0, -1], 
                                      [0, np.sqrt(2), 0]])

    temp_result = np.tensordot(D, C3_stack, axes=([1], [0]))  # Shape: (3, n, m)
    result = np.tensordot(temp_result, D.T, axes=([1], [0]))  # Shape: (3, 3, n, m)

    return np.transpose(result, (0, 3, 1, 2))

def T3_C3_mat(T3_stack):
    # Ensure T3_stack has the correct shape
    assert T3_stack.ndim == 4 and T3_stack.shape[0] == 3 and T3_stack.shape[1] == 3, \
        "T3_stack must have shape (3, 3, n, m)"

    # Define the transformation matrix D
    D = (1 / np.sqrt(2)) * np.array([[1, 0, 1], 
                                      [1, 0, -1], 
                                      [0, np.sqrt(2), 0]])
    
    temp_result = np.tensordot(D.T, T3_stack, axes=([1], [0]))  # Shape: (3, n, m)
    result = np.tensordot(temp_result, D, axes=([1], [0]))  # Shape: (3, 3, n, m)

    return np.transpose(result, (0, 3, 1, 2))
