# This file uses gurobi to solve the primary qcqp !(P)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, C, b, u):
    M = M + 1
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    model = gp.Model() 
    vars = model.addMVar(shape=N, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) 
    model.setObjective(np.ones(N) @ vars, GRB.MINIMIZE)
    
    # Constraints
    for m in range(1, M):
        model.addConstr( vars @ np.array(A[m]) @ vars + 2*np.array(C[m]) @ vars <= b[m] )

    # model.setParam('MIPGap', 0.1)


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
            solution.append(vars[i].X)
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
    with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["C"], load_dict["b"], load_dict["u"])

    