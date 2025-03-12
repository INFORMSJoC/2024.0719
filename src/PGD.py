import sys
import json 
import numpy as np 
import gurobipy as gp
from gurobipy import GRB


import bental_poly_appro 
import primaryProblem 


def calculate_violation(N, M, A, C, b, solution):
    M = M + 1
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    obj = np.array(solution) @ np.array(A[0]) @ np.array(solution) + 2* (np.array(C[0]) @ np.array(solution))
    constraints = [] 
    for m in range(1, M):
        lhs = np.array(solution) @ np.array(A[m]) @ np.array(solution) + 2* (np.array(C[m]) @ np.array(solution)) - b[m]
        constraints.append(lhs)
    return (obj, constraints)  


# def projected_gradient_decent(N, M, A, C, b, cur_solution, gma, bta, alp, dec_rate, tigh_tol, t, T):
#     cur_obj, cur_violation = calculate_violation(N, M, A, C, b, cur_solution)

#     is_violation = False
#     for i in range(len(cur_violation)):
#         if (cur_violation[i] >= tigh_tol):
#             is_violation = True 

#     M = M + 1
#     A.append(np.identity(N).tolist())
#     C.append([0 for i in range(N)])
#     b.append(N*N)

#     model = gp.Model()
#     # Variables
#     d_vars = model.addMVar(shape=N, lb=-N, ub=N, vtype=GRB.CONTINUOUS)  
#     x_vars = model.addMVar(shape=N, lb=-N, ub=N, vtype=GRB.CONTINUOUS) 
#     lmd_vars = model.addMVar(shape=M-1, lb=0, vtype=GRB.CONTINUOUS)
       
#     # Objective 
#     dobj_dx = (2*np.array(A[0])) @ np.array(cur_solution) + 2*np.array(C[0])  
#     if is_violation:
#         model.setObjective(dobj_dx @ d_vars + bta * (d_vars @ d_vars), GRB.MINIMIZE) 
#     else:
#         model.setObjective(dobj_dx @ d_vars, GRB.MINIMIZE)
#         model.addConstr( (d_vars @ d_vars) <= alp * np.exp(-dec_rate*t/T) )
#     # model.setObjective(dobj_dx @ d_vars + gma*(lmd_vars@lmd_vars), GRB.MINIMIZE)
#     # model.addConstr( (d_vars @ d_vars) <= alp * np.exp(-dec_rate*t/T) )


#     model.addConstr(x_vars == np.array(cur_solution) + d_vars)

#     for m in range(1, M):
#         dcon_dx = (2*np.array(A[m])) @ np.array(cur_solution) + 2*np.array(C[m])
#         con_value = np.array(cur_solution) @ np.array(A[m]) @ np.array(cur_solution) + 2* (np.array(C[m]) @ np.array(cur_solution)) - b[m]
#         model.addConstr(con_value + dcon_dx @ d_vars <= tigh_tol )
#         # if con_value <= tigh_tol:
#         #     model.addConstr(con_value + dcon_dx @ d_vars <= tigh_tol )
#         # else:
#         #     model.addConstr(con_value + dcon_dx @ d_vars - lmd_vars[m-1] <= tigh_tol)


#     model.setParam('TimeLimit', 60*60*2) 
#     try:
#         model.optimize()
#     except gp.GurobiError:
#         print("Optimize failed due to non-convexity")
#         model.params.NonConvex = 2
#         model.optimize()

#     solution = []
#     if model.status == GRB.OPTIMAL:
#         for i in range(N):
#             solution.append(float(x_vars[i].X))
#         print("Objective: {0}".format(model.objVal))
#         print("Solution: {0}".format(solution))
#         return (model.objVal, solution)
#     elif model.status == GRB.TIME_LIMIT:
#         print("Objective: {0}".format("TIME_LIMIT"))
#         print("Solution: {0}".format("TIME_LIMIT"))
#         return
#     else:
#         return
    

def projected_gradient_decent_2(N, M, A, C, b, cur_solution, gma, bta, alp, dec_rate, tigh_tol, t, T):
    # cur_obj, cur_violation = calculate_violation(N, M, A, C, b, cur_solution)

    # is_violation = False
    # for i in range(len(cur_violation)):
    #     if (cur_violation[i] >= tigh_tol):
    #         is_violation = True 

    M = M + 1
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    model = gp.Model()
    # Variables
    d_vars = model.addMVar(shape=N, lb=-N, ub=N, vtype=GRB.CONTINUOUS)  
    x_vars = model.addMVar(shape=N, lb=-N, ub=N, vtype=GRB.CONTINUOUS) 
    lmd_vars = model.addMVar(shape=M-1, lb=0, vtype=GRB.CONTINUOUS)
       
    # Objective 
    dobj_dx = (2*np.array(A[0])) @ np.array(cur_solution) + 2*np.array(C[0])  
    # if is_violation:
    #     model.setObjective(dobj_dx @ d_vars + bta * (d_vars @ d_vars), GRB.MINIMIZE) 
    # else:
    #     model.setObjective(dobj_dx @ d_vars, GRB.MINIMIZE)
    #     model.addConstr( (d_vars @ d_vars) <= alp * np.exp(-dec_rate*t/T) )
    model.setObjective( (dobj_dx @ d_vars) + gma * lmd_vars.sum(), GRB.MINIMIZE)
    model.addConstr( (d_vars @ np.identity(N) @ d_vars) <= alp * np.exp(-dec_rate*t/T) )  


    model.addConstr(x_vars == np.array(cur_solution) + d_vars)

    for m in range(1, M):
        dcon_dx = (2*np.array(A[m])) @ np.array(cur_solution) + 2*np.array(C[m])
        con_value = np.array(cur_solution) @ np.array(A[m]) @ np.array(cur_solution) + 2* (np.array(C[m]) @ np.array(cur_solution)) - b[m]
        # model.addConstr(con_value + dcon_dx @ d_vars <= tigh_tol )
        # if con_value <= tigh_tol:
        #     model.addConstr(con_value + dcon_dx @ d_vars <= tigh_tol )
        # else:
        model.addConstr(con_value + dcon_dx @ d_vars  <=  lmd_vars[m-1] + tigh_tol)


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
            solution.append(float(x_vars[i].X))
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        return (model.objVal, solution)
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return


load_dict = 0
pri_obj = 0
pri_solution = 0
lb_obj = 0 
lb_solution = 0 
 
with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
    load_dict = json.load(f)
    # pri_obj, pri_solution = primaryProblem.qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["C"], load_dict["b"], load_dict["u"])
    # lb_obj, lb_solution = bental_poly_appro.qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"], 0.01)
N, M, A, C, b = load_dict["N"], load_dict["M"], load_dict["A"], load_dict["C"], load_dict["b"]


# print('Primary Problem Info: ===========================================')
# print(pri_obj)
# print(pri_solution) 
# _, pri_violation = calculate_violation(N, M, A, C, b, pri_solution) 
# print(pri_violation)
# print('Lower Bound Info: ===========================================')
# print(lb_obj)
# print(lb_solution) 
# _, lb_violation = calculate_violation(N, M, A, C, b, lb_solution) 
# print(lb_violation)

# gma = 1e3
# bta = 1e3
# alp = 1e-3
# dec_rate = 2
# # Tightness tolerance
# tigh_tol = 1e-5
# T = 100


gma = 1e3
bta = 1e3
alp = 1e-3
dec_rate = 2
# Tightness tolerance
tigh_tol = 1e-5
T = 1000


# TODO: from local solution /init solution --> iter,    fmincon
# TODO: dong hongbo

# lb_solution = np.zeros(N).tolist() 
lb_solution = (N * np.random.rand(N)).tolist()  
lb_obj = 0 
_, lb_violation = calculate_violation(N, M, A, C, b, lb_solution) 


cur_solution = lb_solution
obj_record = [lb_obj]
violation_record = [lb_violation]
for t in range(1, T+1):
    _, cur_solution = projected_gradient_decent_2(N, M, A, C, b, cur_solution, gma, bta, alp, dec_rate, tigh_tol, t, T)
    cur_obj, cur_violation = calculate_violation(N, M, A, C, b, cur_solution) 
    obj_record.append(cur_obj)
    violation_record.append(cur_violation)

print(obj_record)
print(f'Last Objective: {lb_obj} ==> {obj_record[-1]}, Primary Objective: {pri_obj}')
print(violation_record[-1])


 