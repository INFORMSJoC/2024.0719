import math
import numpy as np


# how to find the index of the y, \xi, \eta 
def index_y(theta, l, i):
    begin_index = 0
    for index in range(l):
        begin_index += 2**(theta-index)
    return begin_index + i - 1 

def index_xi(theta, v, l, i, j):
    begin_index = 2*(2**theta) - 1
    for l_index in range(1, l):
        for i_index in range(1, 2**(theta-l_index)+1):
            begin_index += 2*v[l_index-1] + 2
    for i_index in range(1, i):
        begin_index += 2*v[l-1] + 2
    return begin_index + j

def index_eta(theta, v, l, i, j):
    return index_xi(theta, v, l, i, j) + v[l-1] + 1



def coefficient_matrix_bental_polyhedral_1(K, epsilon):
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
    
    number_contraints = 0
    number_variables = 2*(2**theta) - 1
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            number_contraints += 4*v[l-1] + 6
            number_variables += 2*v[l-1] + 2

    C = np.zeros( (number_contraints, number_variables) )
    
    begin_index = 0
    for l in range(1, theta+1):
        for i in range(1, 2**(theta-l)+1):
            C[begin_index, index_y(theta, l-1, 2*i-1)] = -1
            C[begin_index, index_xi(theta, v, l, i, 0)] = 1

            C[begin_index+1, index_y(theta, l-1, 2*i-1)] = 1
            C[begin_index+1, index_xi(theta, v, l, i, 0)] = 1

            C[begin_index+2, index_y(theta, l-1, 2*i)] = -1
            C[begin_index+2, index_eta(theta, v, l, i, 0)] = 1

            C[begin_index+3, index_y(theta, l-1, 2*i)] = 1
            C[begin_index+3, index_eta(theta, v, l, i, 0)] = 1

            for j in range(1, v[l-1]+1):
                C[begin_index+3+4*(j-1)+1, index_xi(theta, v, l, i, j)] = 1
                C[begin_index+3+4*(j-1)+1, index_xi(theta, v, l, i, j-1)] = -math.cos( math.pi/(2**(j+1)) ) 
                C[begin_index+3+4*(j-1)+1, index_eta(theta, v, l, i, j-1)] = -math.sin( math.pi/(2**(j+1)) ) 

                C[begin_index+3+4*(j-1)+2, index_xi(theta, v, l, i, j)] = -1
                C[begin_index+3+4*(j-1)+2, index_xi(theta, v, l, i, j-1)] = math.cos( math.pi/(2**(j+1)) ) 
                C[begin_index+3+4*(j-1)+2, index_eta(theta, v, l, i, j-1)] = math.sin( math.pi/(2**(j+1)) ) 

                C[begin_index+3+4*(j-1)+3, index_eta(theta, v, l, i, j)] = 1
                C[begin_index+3+4*(j-1)+3, index_xi(theta, v, l, i, j-1)] = math.sin( math.pi/(2**(j+1)) ) 
                C[begin_index+3+4*(j-1)+3, index_eta(theta, v, l, i, j-1)] = -math.cos( math.pi/(2**(j+1)) ) 

                C[begin_index+3+4*(j-1)+4, index_eta(theta, v, l, i, j)] = 1
                C[begin_index+3+4*(j-1)+4, index_xi(theta, v, l, i, j-1)] = -math.sin( math.pi/(2**(j+1)) ) 
                C[begin_index+3+4*(j-1)+4, index_eta(theta, v, l, i, j-1)] = math.cos( math.pi/(2**(j+1)) ) 
            
            C[begin_index+3+4*v[l-1]+1, index_y(theta, l, i)] = 1
            C[begin_index+3+4*v[l-1]+1, index_xi(theta, v, l, i, v[l-1])] = -1

            C[begin_index+3+4*v[l-1]+2, index_xi(theta, v, l, i, v[l-1])] = math.tan( math.pi/(2**(v[l-1]+1)) )
            C[begin_index+3+4*v[l-1]+2, index_eta(theta, v, l, i, v[l-1])] = -1

            begin_index = begin_index+3+4*v[l-1]+3
    
    # print(C.shape)
    # print(begin_index)
    # print(C)
    P = C[:,:K]
    Q = np.hstack(( C[:,K:2*(2**theta)-2], C[:,2*(2**theta)-1:] ))
    p = C[:,2*(2**theta)-2]
    # print( P.shape )
    # print( Q.shape )
    # print( p.shape )
    return [P, Q, p]


# coefficient_matrix_bental_polyhedral_1(200, 0.1)