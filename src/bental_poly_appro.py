# This file uses gurobi to solve the Y,Z decomposition qcqp !(2)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np
from scipy import sparse

from bental_polyhedral_2 import coefficient_matrix_bental_polyhedral_2


# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, all_Y, all_Z, C, b, u, epslong):
    M = M + 1
    A.append(np.identity(N).tolist())
    all_Y.append(np.identity(N).tolist())
    all_Z.append(np.zeros((N,0)).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    model = gp.Model() 
    vars = model.addMVar(shape=N, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)   
    model.setObjective(np.ones(N) @ vars, GRB.MINIMIZE)


    epslong = float(epslong)

    for m in range(1,M):  
        Y = all_Y[m]
        Z = all_Z[m]
        YT = np.transpose(Y)
        ZT = np.transpose(Z)
        k_upper = len(Y[0]) +1 
        k_low = len(Z[0]) +1 

        y_i = []
        z_i = []
        
        for i in range(k_upper):
            temp = gp.LinExpr()
            if i == k_upper-1:
                number = (b[m]-1) / 2
                for j in range(N):
                    temp += C[m][j]*vars[j]
                y_i.append(number - temp)
            else:
                for j in range(N):
                    temp += YT[i][j]*vars[j]
                y_i.append(temp)

        for i in range(k_low):
            temp = gp.LinExpr()
            if i == k_low-1:
                number = (b[m]+1) / 2
                for j in range(N):
                    temp += C[m][j]*vars[j]
                z_i.append(number - temp)
            else:
                for j in range(N):
                    temp += ZT[i][j]*vars[j]
                z_i.append(temp)
        
        yiMlin = gp.MLinExpr.zeros(len(y_i))
        for i in range(len(y_i)):
            yiMlin[i] = y_i[i] 

        ziMlin = gp.MLinExpr.zeros(len(z_i))
        for i in range(len(z_i)):
            ziMlin[i] = z_i[i]

        t = model.addMVar(shape=1, lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) 

        # temp = gp.QuadExpr()
        # for i in range(k_upper):
        #     temp += y_i[i] * y_i[i]
        # model.addConstr(temp <= t*t)

        ycoe_matrix = coefficient_matrix_bental_polyhedral_2(k_upper, epslong) 
        yP = ycoe_matrix[0]
        yQ = ycoe_matrix[1]
        yp = ycoe_matrix[2]
        # yP = sparse.csr_matrix(yP)
        # yQ = sparse.csr_matrix(yQ)
        # yp = sparse.csr_matrix(yp) 
        y_extra_vars = model.addMVar(shape=yQ.shape[1], vtype=GRB.CONTINUOUS)
        model.addConstr(yP @ yiMlin + yQ @ y_extra_vars + yp.reshape(yP.shape[0],1) @ t >= 0)


        # temp = gp.QuadExpr()
        # for i in range(k_low):
        #     temp += z_i[i] * z_i[i]
        # model.addConstr(temp >= t*t) 

        if k_low == 1:
            model.addConstr(z_i[0]-t >= 0)
            continue

        zcoe_matrix = coefficient_matrix_bental_polyhedral_2(k_low, epslong)
        zP = zcoe_matrix[0]
        zQ = zcoe_matrix[1]
        zp = zcoe_matrix[2]
        # zP = sparse.csr_matrix(zP)
        # zQ = sparse.csr_matrix(zQ)
        # zp = sparse.csr_matrix(zp)
        extra_w = model.addMVar(shape=k_low, vtype=GRB.CONTINUOUS)
        extra_vars = model.addMVar(shape=zQ.shape[1], vtype=GRB.CONTINUOUS)
        extra_v = model.addMVar(shape=zP.shape[0], lb=0, vtype=GRB.CONTINUOUS) 
        extra_bin = model.addMVar(shape=zP.shape[0], vtype=GRB.BINARY)


        model.addConstr(zp.T @ extra_v >= t)
        model.addConstr(zP @ extra_w + zQ @ extra_vars + zp >= 0)
        model.addConstr(zP.T @ extra_v + ziMlin == 0) 
        model.addConstr(zQ.T @ extra_v == 0)
        # model.addConstr(extra_v.T @ (zP @ extra_w + zQ @ extra_vars + zp) == 0)
        for i in range(zP.shape[0]):
            model.addConstr(extra_v[i] - extra_bin[i] * 99999 <= 0)
            model.addConstr( (zP[i,:] @ extra_w + zQ[i,:] @ extra_vars + zp[i]) - (1-extra_bin[i]) * 99999 <= 0)
  

        


    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()

    solution = []
    if model.status == GRB.OPTIMAL:
        for i in range(N):
            solution.append(float(vars[i].X))
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        return (model.objVal, solution)
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return
    



if __name__ == "__main__":
    with open("../data/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"], sys.argv[2])

    