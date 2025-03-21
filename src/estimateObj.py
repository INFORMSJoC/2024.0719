"""
# epslong = 0.01
    # multip_number = math.ceil(math.log(2/epslong))

    # Y = np.identity(N).tolist()
    # Z = [[]]
    # YT = np.transpose(Y)
    # ZT = np.transpose(Z)
    # k_upper = len(Y[0]) +1 
    # k_low = len(Z[0]) +1  
    
    # y_i = []
    # z_i = []
    # for i in range(k_upper):
    #     temp = gp.LinExpr()
    #     if i == k_upper-1:
    #         number = (N*N-1) / 2
    #         y_i.append(number)
    #     else:
    #         for j in range(N):
    #             temp += YT[i][j]*vars[j]
    #         y_i.append(temp)

    # for i in range(k_low):
    #     temp = gp.LinExpr()
    #     if i == k_low-1:
    #         number = (N*N+1) / 2
    #         z_i.append(number)
    #     else:
    #         for j in range(N):
    #             temp += ZT[i][j]*vars[j]
    #         z_i.append(temp)
    
    # t = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    
    # theta = 0
    # for i in range(999):
    #     if 2**i >= k_upper:
    #         theta = i
    #         break
    
    # if theta == 0:
    #     theta = 1
    
    
    
    # y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for l in range(theta+1)] for i in range(2**theta)]   
    # for i in range(k_upper):
    #     y_i_l[i][0] = y_i[i]
    
    # y_i_l[0][theta] = t
        
    # for fl in range(1,theta+1):
    #     v_l = math.floor(fl*multip_number)
    #     kexu_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl))] 
    #     eta_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl))] 
        
    #     for i in range(1, 2**(theta-fl) + 1):
    #         for j in range(v_l+1):

    #             if j == 0:
    #                 model.addLConstr(kexu_i_j[i-1][0] - y_i_l[2*i-2][fl-1], GRB.GREATER_EQUAL, 0)
    #                 model.addLConstr(kexu_i_j[i-1][0] + y_i_l[2*i-2][fl-1], GRB.GREATER_EQUAL, 0)
    #                 model.addLConstr(eta_i_j[i-1][0] - y_i_l[2*i-1][fl-1], GRB.GREATER_EQUAL, 0)
    #                 model.addLConstr(eta_i_j[i-1][0] + y_i_l[2*i-1][fl-1], GRB.GREATER_EQUAL, 0)

    #             else:
    #                 cos = math.cos(math.pi / (2**(j+1)))
    #                 sin = math.sin(math.pi / (2**(j+1)))
    #                 model.addLConstr(kexu_i_j[i-1][j] - cos*kexu_i_j[i-1][j-1] - sin*eta_i_j[i-1][j-1], GRB.EQUAL, 0)
    #                 model.addLConstr(eta_i_j[i-1][j] + sin*kexu_i_j[i-1][j-1] - cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
    #                 model.addLConstr(eta_i_j[i-1][j] - sin*kexu_i_j[i-1][j-1] + cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
                    
    #                 if j == v_l:
    #                     model.addLConstr(kexu_i_j[i-1][j] - y_i_l[i-1][fl], GRB.LESS_EQUAL, 0)
    #                     model.addLConstr(eta_i_j[i-1][j] - math.tan(math.pi / (2**(j+1))) * kexu_i_j[i-1][j], GRB.LESS_EQUAL, 0)

    # if k_low == 1:
    #     model.addLConstr(z_i[0]-t, GRB.GREATER_EQUAL, 0)
"""
import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json
import copy

#  -7.577894982720e+05 -757829.8294234548 -757789.4982720196
def estimate_obj(N, Q, c):
    model = gp.Model()

    # Variables
    vars = []
    for i in range(N): 
        # GRB.INFINITY
        vars.append(model.addVar(lb=-N, ub=N, vtype=GRB.CONTINUOUS))

    origin_Q = copy.deepcopy(Q)
    origin_c = copy.deepcopy(c)
    extra = 0
    e_vals,e_vecs = np.linalg.eig(Q)
    if min(e_vals.real) < 0:
        for i in range(N):
            Q[i][i] += -min(e_vals.real)
        extra = min(e_vals.real)
    
    # Objective
    obj = gp.QuadExpr()
    for i in range(N):
        for j in range(N):
            if Q[i][j] != 0:
                obj += Q[i][j]*vars[i]*vars[j]
    for i in range(N):
        if c[i] != 0:
            obj += 2*c[i]*vars[i]
    obj += extra
    

    # x 的有界  
    bounded = gp.QuadExpr()
    for i in range(N):
        bounded += vars[i] * vars[i]
    model.addConstr(bounded, GRB.LESS_EQUAL, N*N)

    


    model.setParam("OutputFlag", 0)
    model.setObjective(obj, GRB.MINIMIZE)
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        # Solve bilinear model
        model.params.NonConvex = 2
        model.optimize()

    if model.status == GRB.OPTIMAL:
        solution = []
        for i in range(N):
            solution.append(vars[i].X)
        # print("Solution:{}".format(solution))

        origin_obj = 0
        for i in range(N):
            for j in range(N):
                origin_obj += origin_Q[i][j]*solution[i]*solution[j]
        for i in range(N):
            origin_obj += 2*origin_c[i]*solution[i]

        # print("Origin_obj:{}".format(origin_obj))
        # print("model.objVal:{}".format(model.objVal))
        return origin_obj
    

        