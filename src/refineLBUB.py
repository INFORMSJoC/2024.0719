from re import X
import sys 
import json

import numpy as np 
from scipy.optimize import minimize 
from scipy.optimize import Bounds
from scipy.optimize import LinearConstraint
from scipy.optimize import NonlinearConstraint

import gurobipy as gp
from gurobipy import GRB


with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
    load_dict = json.load(f)
    N, M, A, C, b = load_dict["N"], load_dict["M"], load_dict["A"], load_dict["C"], load_dict["b"]
M = M + 1
A.append(np.identity(N).tolist())
C.append([0 for i in range(N)])
b.append(N*N)

def objfun(x):
    return (x @ np.array(A[0]) @ x) + 2 * (np.array(C[0]) @ x)

def objfun_der(x):
    return 2*(x @ np.array(A[0])) + 2*np.array(C[0])

def objfun_hess(x):
    return 2*np.array(A[0])


# linear_constraint = LinearConstraint([[1, 2], [2, 1]], [-np.inf, 1], [1, 1])
bounds = Bounds([-N]*N, [N]*N)

def cons_f(x):
    cons_list = [] 
    for m in range(1,M):
        cons_list.append((x @ np.array(A[m]) @ x) + 2 * (np.array(C[m]) @ x))
    return cons_list
def cons_J(x):
    consJ_list = [] 
    for m in range(1,M):
        consJ_list.append(2*(x @ np.array(A[m])) + 2*np.array(C[m]))
    return consJ_list
def cons_H(x, v):
    result = np.zeros((N,N)) 
    for m in range(1,M):
        result = result + v[m-1] * 2*np.array(A[m])
    return result
nonlinear_constraint = NonlinearConstraint(cons_f, -np.inf, b[1:], jac=cons_J, hess=cons_H)

x0 = np.zeros(N) 
res = minimize(objfun, x0, method='trust-constr', jac=objfun_der, hess=objfun_hess,
               constraints=[nonlinear_constraint],
               options={'verbose': 3}, bounds=bounds)

print(res.fun)
print(res.x) 
# print(cons_f(res.x))
# print(b[1:])


def update_z(x, u, idx):
    """
    x: ndarray
    u: ndarray 
    """
    model = gp.Model()
    # Variables
    z = model.addMVar(shape=N, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)  
       
    # Objective  
    model.setObjective((z-x+u)@(z-x+u), GRB.MINIMIZE)
    
    # Constraints
    model.addConstr((z @ np.array(A[idx]) @ z) + 2 * (np.array(C[idx]) @ z) <= b[idx])

    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()

    if model.status == GRB.OPTIMAL:
        # print("Objective: {0}".format(model.objVal))
        # print("Solution: {0}".format(solution))
        return z.X
    elif model.status == GRB.TIME_LIMIT:
        # print("Objective: {0}".format("TIME_LIMIT"))
        # print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return

def update_y(x, u, rho):
    model = gp.Model()
    # Variables
    y = model.addMVar(shape=N, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)  
       
    # Objective  
    model.setObjective( (x @ np.array(A[0]) @ y) + 2 * (np.array(C[0])@x) + rho*(y-x+u)@(y-x+u), GRB.MINIMIZE)
    
    # Constraints
    for m in range(1, M): 
        model.addConstr(x @ np.array(A[m]) @ y + 2 * (np.array(C[m])@x) <= b[m])


    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()
 
    if model.status == GRB.OPTIMAL: 
        # print("Objective: {0}".format(model.objVal))
        # print("Solution: {0}".format(solution))
        return y.X
    elif model.status == GRB.TIME_LIMIT:
        # print("Objective: {0}".format("TIME_LIMIT"))
        # print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return

def update_x(y, u, rho):
    model = gp.Model()
    # Variables
    x = model.addMVar(shape=N, lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)  
       
    # Objective  
    model.setObjective( (x @ np.array(A[0]) @ y) + 2 * (np.array(C[0])@x) + rho*(y-x+u)@(y-x+u), GRB.MINIMIZE)
    
    # Constraints
    for m in range(1, M): 
        model.addConstr(x @ np.array(A[m]) @ y + 2 * (np.array(C[m])@x) <= b[m])


    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()
 
    if model.status == GRB.OPTIMAL: 
        # print("Objective: {0}".format(model.objVal))
        # print("Solution: {0}".format(solution))
        return x.X
    elif model.status == GRB.TIME_LIMIT:
        # print("Objective: {0}".format("TIME_LIMIT"))
        # print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return 

record = [] 
record_obj = []
record_cons = []

record.append(res.x)
record_obj.append(res.fun)
record_cons.append(cons_f(res.x)-np.array(b[1:]))

rho = 1000

admm_x = res.x  


# update y
# admm_u = np.zeros(N) 

# for i in range(200):
#     admm_y = update_y(admm_x, admm_u, rho)
#     admm_u = admm_u + admm_y - admm_x 
#     admm_x = update_x(admm_y, admm_u, rho)

#     record.append(admm_x)
#     obj = np.array(admm_x) @ np.array(A[0]) @ np.array(admm_x) + 2* (np.array(C[0]) @ np.array(admm_x))
#     constraints = [] 
#     for m in range(1, M):
#         lhs = np.array(admm_x) @ np.array(A[m]) @ np.array(admm_x) + 2* (np.array(C[m]) @ np.array(admm_x)) - b[m]
#         constraints.append(lhs)
    
#     record_obj.append(obj)
#     record_cons.append(constraints)

# print(record_obj)
# print(record_cons[-1])
# print( (admm_x @ np.array(A[0]) @ admm_y) + 2 * (np.array(C[0])@admm_x) + admm_u @ (admm_y-admm_x ) )
 
# print(admm_y - admm_x) 




# update z_1, ..., z_5 

admm_u1 = np.zeros(N) 
admm_u2 = np.zeros(N) 
admm_u3 = np.zeros(N) 
admm_u4 = np.zeros(N) 
admm_u5 = np.zeros(N) 

for i in range(5):
    admm_z1 = update_z(admm_x, admm_u1, 1)
    admm_z2 = update_z(admm_x, admm_u2, 2)
    admm_z3 = update_z(admm_x, admm_u3, 3)
    admm_z4 = update_z(admm_x, admm_u4, 4)
    admm_z5 = update_z(admm_x, admm_u5, 5)

    admm_u1 = admm_u1 + admm_z1 - admm_x 
    admm_u2 = admm_u2 + admm_z2 - admm_x 
    admm_u3 = admm_u3 + admm_z3 - admm_x 
    admm_u4 = admm_u4 + admm_z4 - admm_x 
    admm_u5 = admm_u5 + admm_z5 - admm_x 

    admm_x = np.linalg.inv( np.array(A[0]) + 5 * rho * np.eye(N) )  @  ( -np.array(C[0]) + rho * (admm_u1 + admm_z1 + admm_u2 + admm_z2 + admm_u3 + admm_z3 + admm_u4 + admm_z4 + admm_u5 + admm_z5) )
    record.append(admm_x)
    obj = np.array(admm_x) @ np.array(A[0]) @ np.array(admm_x) + 2* (np.array(C[0]) @ np.array(admm_x))
    constraints = [] 
    for m in range(1, M):
        lhs = np.array(admm_x) @ np.array(A[m]) @ np.array(admm_x) + 2* (np.array(C[m]) @ np.array(admm_x)) - b[m]
        constraints.append(lhs)
    
    record_obj.append(obj)
    record_cons.append(constraints)

print(record_obj)
print(record_cons[-1])
