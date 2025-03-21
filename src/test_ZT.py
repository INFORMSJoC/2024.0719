import gurobipy as gp
from gurobipy import GRB
import numpy as np
# from scipy import sparse
import math

from bental_polyhedral_2 import coefficient_matrix_bental_polyhedral_2


def optimize_primary_problem(N, T, M): 
    model = gp.Model()
    vars = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(N)]

    # vars[0].start = 5 
    # for i in range(1, len(vars)):
    #     vars[i].start = 0 

    # Set Objective Function
    obj_coe = np.ones(N)
    model.setObjective(obj_coe @ vars, GRB.MINIMIZE)  

    # Add Constraints 
    obj = gp.QuadExpr()
    for i in range(N):
        obj += vars[i]*vars[i]
    model.addConstr(obj >= T*T)
    # model.addConstr(obj <= M*M)

    # constraint2 = gp.QuadExpr()
    # for i in range(N):
    #     constraint2 += i*vars[i]*vars[i]
    # model.addConstr(constraint2 >= T*T)
    



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
#     vars = model.addMVar(shape=N, lb=0, vtype=GRB.CONTINUOUS)
#     extra_w = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
#     extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)
#     extra_v = model.addMVar(shape=P.shape[0], lb=0, vtype=GRB.CONTINUOUS)
#     extra_bin = model.addMVar(shape=P.shape[0], vtype=GRB.BINARY)
#     M = 9999 
#     # Set Objective Function 
#     obj = np.ones(N)
#     model.setObjective(obj @ vars, GRB.MINIMIZE)  
 

#     P = sparse.csr_matrix(P)
#     Q = sparse.csr_matrix(Q)
#     # p = sparse.csr_matrix(p)
#     # Add Constraints 
#     model.addConstr(p.T @ extra_v >= T)
#     model.addConstr(P @ extra_w + Q @ extra_vars + p >= 0)
#     model.addConstr(P.T @ extra_v + vars ==0)
#     model.addConstr(Q.T @ extra_v == 0)
#     # for i in range(P.shape[0]):
#         # model.addConstr(extra_v.T @ (P @ extra_w + Q @ extra_vars + p) == 0)
#     model.addConstr(extra_v - (np.ones(P.shape[0])-extra_bin) * M <= 0)
#     model.addConstr(P @ extra_w + Q @ extra_vars + p - extra_bin * M <= 0)



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

# def bounded_v(N, T, M, epsilon, P, Q, p, i):
#     model = gp.Model()
#     vars = model.addMVar(shape=N, lb=0, vtype=GRB.CONTINUOUS) # vars >= 0 
#     extra_w = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
#     extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)
#     extra_v = model.addMVar(shape=P.shape[0], lb=0, vtype=GRB.CONTINUOUS) 

#     model.setObjective(extra_v[i], GRB.MAXIMIZE)  
    
#     model.addConstr(vars @ vars <= M*M)

#     P = sparse.csr_matrix(P)
#     Q = sparse.csr_matrix(Q)
#     # p = sparse.csr_matrix(p)
#     # Add Constraints 
#     model.addConstr(p.T @ extra_v >= T)
#     model.addConstr(P @ extra_w + Q @ extra_vars + p >= 0)
#     model.addConstr(P.T @ extra_v + vars ==0)
#     model.addConstr(Q.T @ extra_v == 0)

#     model.addConstr(p.T @ extra_v <= (1+epsilon) * (1+epsilon) * M) # TODO

#     try:
#         model.optimize()
#     except gp.GurobiError:
#         print("Optimize failed due to non-convexity")
#         model.params.NonConvex = 2
#         model.optimize()
#     return model.objVal

# def bounded_compl_v(N, T, M, epsilon, P, Q, p, i):
#     model = gp.Model()
    
#     extra_w = model.addMVar(shape=N, lb=-(1+epsilon), ub=(1+epsilon), vtype=GRB.CONTINUOUS)
#     extra_vars = model.addMVar(shape=Q.shape[1], lb=0, ub=(1+epsilon), vtype=GRB.CONTINUOUS) 

#     model.setObjective(P[i,:] @ extra_w + Q[i,:] @ extra_vars + p[i], GRB.MAXIMIZE)  
    

#     try:
#         model.optimize()
#     except gp.GurobiError:
#         print("Optimize failed due to non-convexity")
#         model.params.NonConvex = 2
#         model.optimize()
#     return model.objVal

def optimize_bental_polyhedral_2(N, T, M, epsilon): 
    coe_matrix = coefficient_matrix_bental_polyhedral_2(N, epsilon)
    P = coe_matrix[0]
    Q = coe_matrix[1]
    p = coe_matrix[2]
    

    model = gp.Model()
    vars = model.addMVar(shape=N, lb=0, ub=M, vtype=GRB.CONTINUOUS) # M >= vars >= 0 
    extra_w = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
    extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)
    extra_v = model.addMVar(shape=P.shape[0], lb=0, vtype=GRB.CONTINUOUS) 
    extra_bin = model.addMVar(shape=P.shape[0], vtype=GRB.BINARY)

    # UB_v = np.zeros(P.shape[0])
    # for i in range(P.shape[0]):
    #     UB_v[i] = bounded_v(N, T, M, epsilon, P, Q, p, i)

    # UB_compl_v = np.zeros(P.shape[0])
    # for i in range(P.shape[0]):
    #     UB_compl_v[i] = bounded_compl_v(N, T, M, epsilon, P, Q, p, i)
    
    

    # Set Objective Function 
    obj = np.ones(N)
    model.setObjective(obj @ vars, GRB.MINIMIZE)  
    
    # model.addConstr(vars @ vars <= M*M)

    # P = sparse.csr_matrix(P)
    # Q = sparse.csr_matrix(Q)
    # p = sparse.csr_matrix(p)
    # Add Constraints 
    # model.addConstr(np.ones(P.shape[0]) @ extra_bin <= 13)
    model.addConstr(p.T @ extra_v >= T)
    model.addConstr(P @ extra_w + Q @ extra_vars + p >= 0)
    model.addConstr(P.T @ extra_v + vars ==0)
    model.addConstr(Q.T @ extra_v == 0)

    # model.addConstr(extra_v - (np.ones(P.shape[0])-extra_bin) * UB_v <= 0)
    # model.addConstr(P @ extra_w + Q @ extra_vars + p - extra_bin * UB_compl_v <= 0)

    # model.addConstr(extra_v.T @ (P @ extra_w + Q @ extra_vars + p) == 0)
    for i in range(P.shape[0]):
        model.addConstr(extra_v[i] - extra_bin[i] * 9999 <= 0)
        model.addConstr(P[i,:] @ extra_w + Q[i,:] @ extra_vars + p[i] - (1-extra_bin[i]) * 9999 <= 0)
        
        



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



# def optimize_our_polyhedral(N, T, M, epsilon): 
#     # N shoud be the 2**

#     model = gp.Model()
#     z_i = model.addMVar(shape=N, lb=0, ub=M, vtype=GRB.CONTINUOUS) # M >= vars >= 0 
#     # extra_w = model.addMVar(shape=N, vtype=GRB.CONTINUOUS)
#     # extra_vars = model.addMVar(shape=Q.shape[1], vtype=GRB.CONTINUOUS)
#     # extra_v = model.addMVar(shape=P.shape[0], lb=0, vtype=GRB.CONTINUOUS) 
#     # extra_bin = model.addMVar(shape=P.shape[0], vtype=GRB.BINARY)
 
    
    

#     # Set Objective Function 
#     obj = np.ones(N)
#     model.setObjective(obj @ z_i, GRB.MINIMIZE)  
    

    
#     theta = math.ceil(math.log2(N))
#     vvv = []
#     xxx = list(range(1, 21))
#     yyy = [ 1/math.cos( math.pi/(2**(item+1)) ) for item in xxx ] 
#     total_accuracy = 1
#     for layer in range(theta):
#         for i in range(len(yyy)):
#             if total_accuracy * yyy[i] < 1+epsilon:
#                 vvv.append(i+1)
#                 total_accuracy = total_accuracy * yyy[i]
#                 break

#     # 先求层数：
#     layers = 0
#     layers_count = [] # 每一 成的数目
#     layers_variable = [] # 每一程的变量
#     current_layer_count = 2**theta
#     while current_layer_count != 1:
#         layers += 1
#         current_layer_count = current_layer_count // 2 + current_layer_count%2
#         layers_count.append(current_layer_count)
#         layers_variable.append([model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(current_layer_count)])
    
#     # print(layers)
#     # print(layers_count)
#     # print(layers_variable)
#     w_i = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(2**theta)] 
    
#     total_expr = [
#         [ [] for j in range(layers_count[i]) ] for i in range(layers)
#     ]
#     v = [
#         [ [] for j in range(layers_count[i]) ] for i in range(layers)
#     ]
#     i_i = [
#         [ [] for j in range(layers_count[i]) ] for i in range(layers)
#     ]

    
#     layers_variable[layers-1][0] = 1
#     # 第一层
#     for i in range(layers_count[0]):
#         temp_y1 = w_i[2*i]
#         temp_y2 = w_i[2*i+1]
#         temp_y3 = layers_variable[0][i]
#         # j 就是index_mul  v是multip_number
#         for index_mul in range(0, 2**vvv[0]+1):
#             model.addConstr( math.cos(math.pi*index_mul / (2**(vvv[0]+1))) * temp_y1 + math.sin(math.pi*index_mul / (2**(vvv[0]+1)))*temp_y2<= temp_y3)
#             total_expr[0][i].append(-math.cos(math.pi*index_mul / (2**(vvv[0]+1))) * temp_y1 - math.sin(math.pi*index_mul / (2**(vvv[0]+1)))*temp_y2 + temp_y3)
#             v[0][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
#             i_i[0][i].append(model.addVar(vtype=GRB.BINARY))
    
#     for lay in range(1, layers):
#         for i in range(layers_count[lay]):
#             temp_y1 = layers_variable[lay-1][2*i]
#             temp_y2 = layers_variable[lay-1][2*i+1]
#             temp_y3 = layers_variable[lay][i]
#             # j 就是index_mul  v是multip_number
#             for index_mul in range(0, 2**vvv[lay]+1):
#                 model.addConstr( math.cos(math.pi*index_mul / (2**(vvv[lay]+1))) * temp_y1 + math.sin(math.pi*index_mul / (2**(vvv[lay]+1)))*temp_y2<= temp_y3)
#                 total_expr[lay][i].append(-math.cos(math.pi*index_mul / (2**(vvv[lay]+1))) * temp_y1 - math.sin(math.pi*index_mul / (2**(vvv[lay]+1)))*temp_y2+ temp_y3)
#                 v[lay][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
#                 i_i[lay][i].append(model.addVar(vtype=GRB.BINARY))

#     # TODO:
#     M_bar = 9999 #9999 # GRB.INFINITY

#     for i in range(layers):
#         for j in range(layers_count[i]):
#             bound_constr = gp.LinExpr()
#             for k in range(2**vvv[i]+1):
#                 # bound_constr += i_i[i][j][k]
#                 model.addConstr(v[i][j][k] <= M_bar*i_i[i][j][k])
#                 model.addConstr(total_expr[i][j][k] <= M_bar*( 1-i_i[i][j][k] ))
#             # model.addConstr(bound_constr <= 2)

    
    
#     # 对于第一层
#     for index in range(layers_count[0]):
#         aT = []
#         for j in range(0, 2**vvv[0]+1):
#             aT.append(math.cos(math.pi*j / (2**(vvv[0]+1))))
#         bT = []
#         for j in range(0, 2**vvv[0]+1):
#             bT.append(math.sin( math.pi*j / (2**(vvv[0]+1)) ))
#         # constraints_number += 2
#         aTv = gp.LinExpr()
#         bTv = gp.LinExpr()
#         for temp_i in range(0, 2**vvv[0]+1):
#             aTv += aT[temp_i]*v[0][index][temp_i]
#             bTv += bT[temp_i]*v[0][index][temp_i]
        
#         model.addConstr(aTv >= -z_i[2*index])
#         model.addConstr(aTv >= z_i[2*index])

#         temp_var1 = model.addVar(vtype=GRB.BINARY)
#         abs_zi = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
#         model.addConstr(abs_zi + z_i[2*index] <= M_bar*temp_var1)
#         model.addConstr(abs_zi + z_i[2*index] >= -M_bar*temp_var1)
#         model.addConstr(abs_zi - z_i[2*index] <= M_bar*(1-temp_var1))
#         model.addConstr(abs_zi - z_i[2*index] >= -M_bar*(1-temp_var1))
#         temp_var2 = model.addVar(vtype=GRB.BINARY)
#         model.addConstr(aTv - abs_zi <= M_bar*temp_var2)
#         model.addConstr(w_i[2*index] <= M_bar*(1-temp_var2))

        
#         model.addConstr(bTv >= -z_i[2*index+1])
#         model.addConstr(bTv >= z_i[2*index+1])

#         temp_var1 = model.addVar(vtype=GRB.BINARY)
#         abs_zi = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
#         model.addConstr(abs_zi + z_i[2*index+1] <= M_bar*temp_var1)
#         model.addConstr(abs_zi + z_i[2*index+1] >= -M_bar*temp_var1)
#         model.addConstr(abs_zi - z_i[2*index+1] <= M_bar*(1-temp_var1))
#         model.addConstr(abs_zi - z_i[2*index+1] >= -M_bar*(1-temp_var1))
#         temp_var2 = model.addVar(vtype=GRB.BINARY)
#         model.addConstr(bTv - abs_zi <= M_bar*temp_var2)
#         model.addConstr(w_i[2*index+1] <= M_bar*(1-temp_var2))


    
#     for layer_index in range(layers-1):
#         for index in range(layers_count[layer_index]):
#             # constraints_number += 1
#             eTv = gp.LinExpr()
#             for temp_i in range(2**vvv[layer_index]+1):
#                 eTv += -v[layer_index][index][temp_i]
#             cTv = gp.LinExpr()
#             # 46 
#             if index % 2 == 0: # 0 2 4 6 
#                 if (index+1+1)//2 % 2 != 0: # 奇数 0  4 
#                     for temp_i in range(2**vvv[layer_index]+1):
#                         cTv += aT[temp_i]*v[layer_index+1][index//2][temp_i]
#                 else: # 2 6
#                     for temp_i in range(2**vvv[layer_index]+1):
#                         cTv += aT[temp_i]*v[layer_index+1][index//2][temp_i]
#             # 47
#             else: # 1 3 5 7
#                 if (index+1)//2 % 2 != 0: # 奇数 1 5
#                     for temp_i in range(2**vvv[layer_index]+1):
#                         cTv += bT[temp_i]*v[layer_index+1][(index-1)//2][temp_i]
#                 else:# 3 7
#                     for temp_i in range(2**vvv[layer_index]+1):
#                         cTv += bT[temp_i]*v[layer_index+1][(index-1)//2][temp_i]
#             # KKT 以及 本来的约束
#             model.addConstr(eTv+cTv >= 0)
#             temp_var = model.addVar(vtype=GRB.BINARY)
#             model.addConstr(eTv+cTv <= M_bar*temp_var)
#             model.addConstr(layers_variable[layer_index][index] <= M_bar*(1-temp_var))
            
            
    

#     # 最后一层
#     min_eTv = gp.LinExpr()
#     for i in range(2**vvv[theta-1]+1):
#         min_eTv += v[-1][0][i]
#     model.addConstr(min_eTv - T>=0)



#     model.setParam('TimeLimit', 60*60*2) 
#     try:
#         model.optimize()
#     except gp.GurobiError:
#         print("Optimize failed due to non-convexity")
#         model.params.NonConvex = 2
#         model.optimize()

#     result = 0
#     solution = z_i.X
#     if model.status == GRB.OPTIMAL:
#         result = model.objVal 
#         print("Objective: {0}".format(model.objVal))
#         print("Solution: {0}".format(z_i.X))
#     else:
#         print("Model cannot be optimized")
    
#     return (result, solution)


N = 10
T = 5 
M = 10000
epsilon = 0.001  

# N = 4                   7
# 4.809698831278216 # 4.678085115435202
# 4.951963201008076 # 4.9515034609403985
# 4.996739006679395 # 

# pri_obj, pri_solution = optimize_primary_problem(N, T, M)

# poly1_obj, poly1_solution = optimize_bental_polyhedral_1(N, T, epsilon)

poly2_obj, poly2_solution = optimize_bental_polyhedral_2(N, T, M, epsilon)

# poly3_obj, poly3_solution = optimize_our_polyhedral(N, T, M, epsilon)
# {poly1_obj}
# {poly2_obj}
print(f"""Objective    Solution 
  {poly2_obj}""")