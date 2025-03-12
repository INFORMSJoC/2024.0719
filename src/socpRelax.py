# This file uses gurobi to solve the primary qcqp !(P)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, C, b, u): 

    model = gp.Model()
    # Variables
    x_vars = model.addMVar(shape=N, lb=-N, ub=N, vtype=GRB.CONTINUOUS) 
    X_vars = model.addMVar(shape=(N,N), vtype=GRB.CONTINUOUS)
     

    temp = gp.LinExpr()
    for i in range(N):
        temp += x_vars[i]
 
    model.setObjective(temp, GRB.MINIMIZE)
    
    # Constraints
    for m in range(1, M):
        eigenvalues, eigenvectors = np.linalg.eig(np.array(A[m]))
        if min(eigenvalues) >= 0: 
            constrain = gp.QuadExpr()
            for i in range(N):
                for j in range(N):
                    constrain += A[m][i][j]*x_vars[i]*x_vars[j]
            for i in range(N):
                constrain += 2*C[m][i] * x_vars[i]
            model.addConstr(constrain <= b[m])
        else:
            constrain = gp.LinExpr()
            for i in range(N):
                for j in range(N):
                    constrain += A[m][i][j]*X_vars[i][j]
            for i in range(N):
                constrain += 2*C[m][i] * x_vars[i]
            model.addConstr(constrain <= b[m])

            temp_A = np.zeros((N,N))
            for l in range(N):
                if eigenvalues[l] >= 0:
                    temp_A += eigenvalues[l] * (eigenvectors[:,l].reshape(N,1) @ eigenvectors[:,l].reshape(1,N))
                else:
                    ujujt = eigenvectors[:,l].reshape(N,1) @ eigenvectors[:,l].reshape(1,N)
                    constrain = gp.QuadExpr()
                    for i in range(N):
                        for j in range(N):
                            constrain += ujujt[i][j]*x_vars[i]*x_vars[j]
                    for i in range(N):
                        for j in range(N):
                            constrain -= ujujt[i][j]*X_vars[i][j]
                    model.addConstr(constrain <= 0)
            constrain = gp.QuadExpr()
            for i in range(N):
                for j in range(N):
                    constrain += temp_A[i][j]*x_vars[i]*x_vars[j]
            for i in range(N):
                for j in range(N):
                    constrain -= temp_A[i][j]*X_vars[i][j]
            model.addConstr(constrain <= 0)




        


    model.setParam('TimeLimit', 60*60*2) 
    # try:
    model.optimize()
    # except gp.GurobiError:
    #     print("Optimize failed due to non-convexity")
    #     model.params.NonConvex = 2
    #     model.optimize()

    solution = []
    if model.status == GRB.OPTIMAL:
        for i in range(N):
            solution.append(x_vars[i].X)
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

    