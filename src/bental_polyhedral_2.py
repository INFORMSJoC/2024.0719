import math
import numpy as np


# how to find the index of the y, \xi, \eta 
# index of y = index - 1
# number of xi: 2*(2**theta) - 1 - 2**theta
# number of eta: (2*(2**theta) - 1 - 2**theta) * (v[0]+1)
def index_xi(theta, l, i):
    begin_index = 2**theta + 1
    for l_index in range(1, l):
        for i_index in range(1, 2**(theta-l_index)+1):
            begin_index += 1
    return begin_index + i - 1



def index_eta(theta, v, l, i, j):
    begin_index = 2*(2**theta)   
    for l_index in range(1, l):
        for i_index in range(1, 2**(theta-l_index)+1):
            begin_index += v[l_index-1] + 1
    for i_index in range(1, i):
        begin_index += v[l-1] + 1
    return begin_index + j



def coefficient_matrix_bental_polyhedral_2(K, epsilon):
    theta = math.ceil(math.log2(K))  
    v = []
    xxx = list(range(1, 21))
    yyy = [ 1/math.cos( math.pi/(2**(item+1)) ) for item in xxx ] 
    total_accuracy = 1
    for layer in range(theta):
        for i in range(len(yyy)):
            if total_accuracy * yyy[i] < 1+epsilon:
                v.append(i+1)
                total_accuracy = total_accuracy * yyy[i]
                break
    # for v_value in range(1, 15):
    #     layer_accuracy = [1/math.cos( math.pi/(2**(v_value+1)) ) for i in range(1, theta+1)]
    #     total_accuracy = 1
    #     for i in layer_accuracy:
    #         total_accuracy = total_accuracy * i 
    #     if total_accuracy <= 1+epsilon:
    #         break
    # v = [v_value  for i in range(1, theta+1)]

    number_variables = 2**theta # layer 0: y 
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_variables += 2
            for j in range(1, v[l-1]+1):
                number_variables += 1
    number_variables += 1 # t

    number_contraints = 0
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            if l == 1:
                number_contraints += 4
            else:
                number_contraints += 2
    
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_contraints += 1
            for j in range(1, v[l-1]+1):
                number_contraints += 2
    number_contraints += 1


    C = np.zeros( (number_contraints, number_variables) )

    begin_index = 0 
    for i in range(1, 2**(theta-1)+1):
        C[begin_index, index_xi(theta, 1, i)] = 1
        C[begin_index, 2*i-2] = -1 

        C[begin_index+1, index_xi(theta, 1, i)] = 1
        C[begin_index+1, 2*i-2] = 1 

        C[begin_index+2, index_eta(theta, v, 1, i, 0)] = 1
        C[begin_index+2, 2*i-1] = -1 

        C[begin_index+3, index_eta(theta, v, 1, i, 0)] = 1
        C[begin_index+3, 2*i-1] = 1 
        begin_index += 4
    
    for l in range(2, theta+1):
        for i in range(1, 2**(theta-l)+1):
            C[begin_index, index_xi(theta, l, i)] = math.sin( math.pi/ (2**(v[l-2]+1)) )
            C[begin_index, index_xi(theta, l-1, 2*i-1)] = -1/ (2**v[l-2])
            for k in range(v[l-2]):
                C[begin_index, index_eta(theta, v, l-1, 2*i-1, k)] = -math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(v[l-2]-1-k) )

            C[begin_index+1, index_eta(theta, v, l, i, 0)] = math.sin( math.pi/ (2**(v[l-2]+1)) )
            C[begin_index+1, index_xi(theta, l-1, 2*i)] = -1/ (2**v[l-2])
            for k in range(v[l-2]):
                C[begin_index+1, index_eta(theta, v, l-1, 2*i, k)] = -math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(v[l-2]-1-k) )

            begin_index += 2

    C[begin_index, 2**theta] = math.sin( math.pi/ (2**(v[theta-1]+1)) )
    C[begin_index, index_xi(theta, theta, 1)] = -1/ (2**v[theta-1])
    for k in range(v[theta-1]):
        C[begin_index, index_eta(theta, v, theta, 1, k)] = -math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(v[theta-1]-1-k) )
    begin_index += 1
    
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            for j in range(1, v[l-1]+1):
                C[begin_index, index_eta(theta, v, l, i, j)] = 2*math.cos( math.pi/ (2**(j+1)) )
                C[begin_index, index_eta(theta, v, l, i, j-1)] = -2*math.cos( math.pi/ (2**(j+1)) )*math.cos( math.pi/ (2**(j+1)) )
                C[begin_index, index_xi(theta, l, i)] = 1/ (2**(j-1))  
                for k in range(j-1):
                    C[begin_index, index_eta(theta, v, l, i, k)] = math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(j-2-k) )

                C[begin_index+1, index_eta(theta, v, l, i, j)] = 2*math.cos( math.pi/ (2**(j+1)) )
                C[begin_index+1, index_eta(theta, v, l, i, j-1)] = 2*math.cos( math.pi/ (2**(j+1)) )*math.cos( math.pi/ (2**(j+1)) )
                C[begin_index+1, index_xi(theta, l, i)] = -1/ (2**(j-1))  
                for k in range(j-1):
                    C[begin_index+1, index_eta(theta, v, l, i, k)] = -math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(j-2-k) )

                begin_index += 2

    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            C[begin_index, index_eta(theta, v, l, i, v[l-1])] = -math.cos( math.pi/ (2**(v[l-1]+1)) )
            C[begin_index, index_xi(theta, l, i)] = 1/ (2**v[l-1])
            for k in range(v[l-1]):
                C[begin_index, index_eta(theta, v, l, i, k)] = math.sin( math.pi/(2**(k+2)) ) * math.sin( math.pi/(2**(k+2)) ) / ( 2**(v[l-1]-1-k) )
            begin_index += 1

    








    
    
    # print(C.shape)
    # print(begin_index) 

    P = C[:,:K]
    Q = np.hstack(( C[:,K:2**theta], C[:, 2**theta+1:] ))
    p = C[:,2**theta]
    # print( P.shape )
    # print( Q.shape )
    # print( p.shape )
    return [P, Q, p]


# coefficient_matrix_bental_polyhedral_2(4, 0.1)