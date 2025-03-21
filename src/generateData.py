from scipy.stats import ortho_group
import numpy as np
import random
import json



# number of variables
N = 40
# number of positive eigenvalues
positive_number = 37
# number of constraints
M = 5
# number of convex constraints 
M_1 = 3
# number of nonconvex constraints 
M_2 = 2


# boundary constraint
data_bound = 50

# generate 10 samples
for data_index in range(3):
    A = []
    C = []
    all_Y = []
    all_Z = []
    for m in range(M_1):
        Y_vecs = ortho_group.rvs(dim=N)
        Z_vecs = ortho_group.rvs(dim=N)

        e_vals = [0 for i in range(N)]
        for i in range(N):
            e_vals[i] = random.randrange(1,data_bound)

        Y_vals = []
        Z_vals = []
        for i in range(N):
            if e_vals[i] >= 0:
                Y_vals.append(e_vals[i])
                Z_vals.append(0)
            else:
                Y_vals.append(0)
                Z_vals.append(-e_vals[i])

        Y = np.dot(Y_vecs, np.diag(Y_vals))
        Z = np.dot(Z_vecs, np.diag(Z_vals))
        Y = Y[:,~(Y==0).all(0)]
        Z = Z[:,~(Z==0).all(0)]
        temp = Y.dot(Y.T)-Z.dot(Z.T)
        
        A.append(temp.tolist())
        all_Y.append(Y.tolist())
        all_Z.append(Z.tolist())

        ctemp = [0 for i in range(N)]
        for i in range(N):
            ctemp[i] = random.randrange(-data_bound,data_bound)
        C.append(ctemp)

    for m in range(M_2):
        Y_vecs = ortho_group.rvs(dim=N)
        Z_vecs = ortho_group.rvs(dim=N)
        e_vals = [0 for i in range(N)]

        for i in range(positive_number):
            e_vals[i] = random.randrange(1,data_bound)
        for i in range(positive_number, N):
            e_vals[i] = random.randrange(-data_bound,-1)

        Y_vals = []
        Z_vals = []
        for i in range(N):
            if e_vals[i] >= 0:
                Y_vals.append(e_vals[i])
                Z_vals.append(0)
            else:
                Y_vals.append(0)
                Z_vals.append(-e_vals[i])

        Y = np.dot(Y_vecs, np.diag(Y_vals))
        Z = np.dot(Z_vecs, np.diag(Z_vals))
        Y = Y[:,~(Y==0).all(0)]
        Z = Z[:,~(Z==0).all(0)]
        temp = Y.dot(Y.T)-Z.dot(Z.T)
        
        A.append(temp.tolist())
        all_Y.append(Y.tolist())
        all_Z.append(Z.tolist())

        ctemp = [0 for i in range(N)]
        for i in range(N):
            ctemp[i] = random.randrange(-data_bound,data_bound)
        C.append(ctemp)

    b = [0 for i in range(M)]
    for i in range(M):
        b[i] = random.randrange(-data_bound,data_bound) + 1000

    data_dict = {
        "N" : N,
        "M" : M,
        "A" : A,
        "Y" : all_Y,
        "Z" : all_Z,
        "C" : C,
        "b" : b,
        "u" : 999999
    }
    
    with open("./QCQPData/" + str(N) + "_" + str(M_1) + "_" + str(M_2) + "_" + str(positive_number) + "_" + str(data_index) + ".json","w") as f:
        json.dump(data_dict,f)

