import gurobipy as gp
from gurobipy import GRB
import numpy as np


# from bental_polyhedral_1 import coefficient_matrix_bental_polyhedral_1
from bental_polyhedral_2 import coefficient_matrix_bental_polyhedral_2


def optimize_primary_problem(N, T): 
    model = gp.Model()
    vars = [model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(N)]

    # Set Objective Function
    obj_coe = np.ones(N)
    model.setObjective(obj_coe @ vars, GRB.MAXIMIZE)  

    # Add Constraints 
    obj = gp.QuadExpr()
    for i in range(N):
        obj += vars[i]*vars[i]
    model.addConstr(obj <= T*T)


    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()

    result = 0
    solution = []
    if model.status == GRB.OPTIMAL:
        result = model.objVal
        for i in range(N):
            solution.append(vars[i].X)
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
    else:
        print("Model cannot be optimized")
    
    return (result, solution)



# def optimize_bental_polyhedral_1(N, T, epsilon): 
#     coe_matrix = coefficient_matrix_bental_polyhedral_1(N, epsilon)
#     P = coe_matrix[0]
#     Q = coe_matrix[1]
#     p = coe_matrix[2]

    

#     model = gp.Model()
#     vars = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
#     extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)

#     # Set Objective Function
#     obj = gp.QuadExpr()
#     for i in range(N):
#         obj += vars[i]*vars[i]
#     obj = np.ones(N)
#     model.setObjective(obj @ vars, GRB.MAXIMIZE)  

#     # Add Constraints 
#     model.addConstr(P @ vars + Q @ extra_vars + T*p >= 0)


#     model.setParam('TimeLimit', 60*60*2) 
#     try:
#         model.optimize()
#     except gp.GurobiError:
#         print("Optimize failed due to non-convexity")
#         model.params.NonConvex = 2
#         model.optimize()

#     result = 0
#     solution = vars.X
#     if model.status == GRB.OPTIMAL:
#         result = model.objVal 
#         print("Objective: {0}".format(model.objVal))
#         print("Solution: {0}".format(vars.X))
#     else:
#         print("Model cannot be optimized")
    
#     return (result, solution)


def optimize_bental_polyhedral_2(N, T, epsilon): 
    coe_matrix = coefficient_matrix_bental_polyhedral_2(N, epsilon)
    P = coe_matrix[0]
    Q = coe_matrix[1]
    p = coe_matrix[2]

    

    model = gp.Model()
    vars = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
    extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)

    # Set Objective Function
    # obj = gp.QuadExpr()
    # for i in range(N):
    #     obj += vars[i]*vars[i]
    obj = np.ones(N)
    model.setObjective(obj @ vars, GRB.MAXIMIZE)  

    # Add Constraints 
    model.addConstr(P @ vars + Q @ extra_vars + T*p >= 0)


    model.setParam('TimeLimit', 60*60*2) 
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        model.params.NonConvex = 2
        model.optimize()

    result = 0
    solution = vars.X
    if model.status == GRB.OPTIMAL:
        result = model.objVal 
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(vars.X))
    else:
        print("Model cannot be optimized")
    
    return (result, solution)




N = 100
T = 10 
epsilon = 0.001

pri_obj, pri_solution = optimize_primary_problem(N, T)

# poly1_obj, poly1_solution = optimize_bental_polyhedral_1(N, T, epsilon)

poly2_obj, poly2_solution = optimize_bental_polyhedral_2(N, T, epsilon)

print(f"""
Objective       Solution 
{pri_obj}         
{poly2_obj}      
""")