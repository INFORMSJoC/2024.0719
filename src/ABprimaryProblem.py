# This file uses gurobi to solve the A,B decomposition qcqp !(1)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, all_Y, all_Z, C, b, u):
    M = M + 1
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
    AAtranspose = np.dot(all_Y[0], np.transpose(all_Y[0]))
    BBtranspose = np.dot(all_Z[0], np.transpose(all_Z[0]))
    for i in range(N):
        for j in range(N):
            temp += AAtranspose[i][j]*vars[i]*vars[j]
    for i in range(N):
        for j in range(N):
            temp += -BBtranspose[i][j]*vars[i]*vars[j]
    for i in range(N):
        temp += 2*C[0][i]*vars[i]

    # model.addConstr(temp, GRB.LESS_EQUAL, obj)
    model.setObjective(temp, GRB.MINIMIZE)
    # Constraints
    for m in range(1, M):
        constrain = gp.QuadExpr()
        AAtranspose = np.dot(all_Y[m], np.transpose(all_Y[m]))
        BBtranspose = np.dot(all_Z[m], np.transpose(all_Z[m]))
        for i in range(N):
            for j in range(N):
                constrain += AAtranspose[i][j]*vars[i]*vars[j]
        for i in range(N):
            for j in range(N):
                constrain += -BBtranspose[i][j]*vars[i]*vars[j]
        for i in range(N):
            constrain += 2*C[m][i] * vars[i]
        model.addConstr(constrain, GRB.LESS_EQUAL, b[m])


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
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"])

    