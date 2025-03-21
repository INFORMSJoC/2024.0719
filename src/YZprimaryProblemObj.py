# This file uses gurobi to solve the Y,Z decomposition qcqp !(2)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np
import math

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, all_Y, all_Z, C, b, u):
    M = M + 1
    A.append(np.identity(N).tolist())
    all_Y.append(np.identity(N).tolist())
    all_Z.append(np.zeros((N,0)).tolist())
    C.append([0 for i in range(N)])
    b.append(u)

    model = gp.Model()
    # Variables
    vars = [model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(N)]
       
    # Objective constraint
    # obj = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    

    temp = gp.QuadExpr()
    for i in range(N):
        for j in range(N):
            temp += A[0][i][j]*vars[i]*vars[j]
    for i in range(N):
        temp += 2*C[0][i]*vars[i]

    # model.addConstr(temp, GRB.LESS_EQUAL, obj)

    model.setObjective(temp, GRB.MINIMIZE)

    # epslong = float(epslong)
    # multip_number = math.ceil(math.log(2/epslong))
    t = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(M)]

    for m in range(1,M):
        Y = all_Y[m]
        Z = all_Z[m]
        YT = np.transpose(Y)
        ZT = np.transpose(Z)
        k_upper = len(Y[0]) +1 
        k_low = len(Z[0]) +1 

        y_i = []
        z_i = []
        if m != 0:
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
        else:
            for i in range(k_upper):
                temp = gp.LinExpr()
                if i == k_upper-1:
                    for j in range(N):
                        temp += C[m][j]*vars[j]
                    y_i.append(((obj-1) / 2) - temp)
                else:
                    for j in range(N):
                        temp += YT[i][j]*vars[j]
                    y_i.append(temp)

            for i in range(k_low):
                temp = gp.LinExpr()
                if i == k_low-1:
                    for j in range(N):
                        temp += C[m][j]*vars[j]
                    z_i.append(((obj+1) / 2) - temp)
                else:
                    for j in range(N):
                        temp += ZT[i][j]*vars[j]
                    z_i.append(temp)

        temp = gp.QuadExpr()
        for i in range(k_upper):
            temp += y_i[i] * y_i[i]
        model.addConstr(temp, GRB.LESS_EQUAL, t[m]*t[m])
        
        temp = gp.QuadExpr()
        for i in range(k_low):
            temp += z_i[i] * z_i[i]
        model.addConstr(temp, GRB.GREATER_EQUAL, t[m]*t[m])



    model.setParam('TimeLimit', 60*60*3) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()

    solution = []
    if model.status == GRB.OPTIMAL:
        for i in range(N):
            solution.append(vars[i].X)
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        return
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return
    



if __name__ == "__main__":
    with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"])

    