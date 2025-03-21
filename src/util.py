import numpy as np
from scipy import linalg 

def decompose_sysmtic_matrix(N, matrix):
    e_vals,e_vecs = linalg.eig(matrix)
    e_vals = e_vals.real
    e_vecs = e_vecs.real
    AAT = np.zeros((N,N))
    BBT = np.zeros((N,N))
    for i in range(N):
        if e_vals[i] >= 0:
            AAT += np.dot(  np.dot(e_vecs[:,i:i+1], np.diag([e_vals[i]])), np.transpose(e_vecs[:,i:i+1])  )
        else:
            BBT += np.dot(  np.dot(e_vecs[:,i:i+1], np.diag([-e_vals[i]])), np.transpose(e_vecs[:,i:i+1])  )
    
    a_vals_r,a_vecs = np.linalg.eig(AAT)
    b_vals_r,b_vecs = np.linalg.eig(BBT)
    a_vals_r = a_vals_r.real
    b_vals_r = b_vals_r.real
    a_vecs = a_vecs.real
    b_vecs = b_vecs.real
    a_vals = []
    b_vals = []
    for i in a_vals_r:
        if i < 0.00000001:
            a_vals.append(0)
        else:
            a_vals.append(i**0.5)

    for i in b_vals_r:
        if i < 0.00000001:
            b_vals.append(0)
        else:
            b_vals.append(i**0.5)

    A = np.dot( a_vecs, np.diag(a_vals) )
    B = np.dot( b_vecs, np.diag(b_vals) )
    A = A[:,~(A==0).all(0)]
    B = B[:,~(B==0).all(0)]
   
    return A.tolist(),B.tolist()

if __name__ == "__main__":
    N = 3 
    matrix1 = [[1,2,3], 
                [2,3,1], 
                [3,1, 6 ]]

    A, B = decompose_sysmtic_matrix(N, matrix1) 

    print(A@np.transpose(A) - B@np.transpose(B)) 