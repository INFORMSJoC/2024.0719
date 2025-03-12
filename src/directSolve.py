#!/usr/bin/env/python3
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, Q, c, A, C, b, u):
    M = M + 1
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    model = gp.Model()
    # Variables
    vars = []
    for i in range(N): 
        vars.append(model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))

    # Objective 
    obj = gp.QuadExpr()
    for i in range(N):
        for j in range(N):
            if Q[i][j] != 0:
                obj += Q[i][j]*vars[i]*vars[j]
    for i in range(N):
        if c[i] != 0:
            obj += 2*c[i]*vars[i]

    model.setObjective(obj, GRB.MINIMIZE)

    # Constraints
    for m in range(M):
        constrain = gp.QuadExpr()
        for i in range(N):
            for j in range(N):
                if A[m][i][j] != 0:
                    constrain += A[m][i][j]*vars[i]*vars[j]
        for i in range(N):
            if C[m][i] != 0:
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
        for m in range(M):
            constrains = 0
            for i in range(N):
                for j in range(N):
                    if A[m][i][j] != 0:
                        constrains += A[m][i][j]*solution[i]*solution[j]
            for i in range(N):
                if C[m][i] != 0:
                    constrains += 2*C[m][i] * solution[i]
            print("constrain{0}: {1}".format(m, constrains))
            print("b{0}: {1}".format(m, b[m]))
        return solution
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return False
    else:
        # print("Objective: {0}".format(model.status))
        # print("Solution: {0}".format(model.status))
        return False


if __name__ == "__main__":
    with open("./QCQPData" + str(sys.argv[3]) + "/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["Q"], load_dict["c"], load_dict["A"], load_dict["C"], load_dict["b"], load_dict["u"])

    