#!/usr/bin/env python3
import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json



# N: varible count, M: constrain count, A: all martric: m*N*N, C: m* N, b m
def v_total_max(N, M, Q, c, A, all_Y, all_Z, C, b, u, epslong, m_index, obj_max):
    b[M-1] = obj_max

    model = gp.Model()
    # Variables
    vars = []
    for i in range(N): 
        vars.append(model.addVar(lb=-N, ub=N, vtype=GRB.CONTINUOUS))

    epslong = float(epslong)
    multip_number = math.ceil(math.log(2/epslong))

    m = m_index
    Y = all_Y[m]
    Z = all_Z[m]
    YT = np.transpose(Y)
    ZT = np.transpose(Z)
    k_upper = len(Y[0]) +1 
    k_low = len(Z[0]) +1  
    
    y_i = []
    z_i = []
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

    
    theta = 0
    for i in range(999):
        if 2**i >= k_low:
            theta = i
            break

    w_i = [model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(k_low)]
    
    y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for l in range(theta+1)] for i in range(2**theta)]   
    for i in range(k_low):
        y_i_l[i][0] = w_i[i]
    
    y_i_l[0][theta] = 1

    v_number = 0 
    for fl in range(1, theta+1):
        v_l = math.floor(fl*multip_number)
        for i in range(1, 2**(theta-fl) + 1):
            for j in range(v_l+1):
                if j == v_l:
                    v_number += 6
                else:
                    v_number += 4


    v = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(v_number)]
    total_expr = []
    P = []
    p = []
    max_Q = math.floor(2*theta* (2**theta) *(theta*multip_number+1) + (2**theta) * (theta+1))
    # max_Q = math.floor(2*theta* (2**(theta-1)) *(theta*multip_number+1) + (2**theta) * (theta+1))
    Q = []

    for fl in range(1,theta+1):
        v_l = math.floor(fl*multip_number)
        kexu_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl))] 
        eta_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl))] 
        
        for i in range(1, 2**(theta-fl) + 1):
            for j in range(v_l+1):

                if j == 0:
                    model.addLConstr(kexu_i_j[i-1][0] - y_i_l[2*i-2][fl-1], GRB.GREATER_EQUAL, 0)
                    model.addLConstr(kexu_i_j[i-1][0] + y_i_l[2*i-2][fl-1], GRB.GREATER_EQUAL, 0)
                    model.addLConstr(eta_i_j[i-1][0] - y_i_l[2*i-1][fl-1], GRB.GREATER_EQUAL, 0)
                    model.addLConstr(eta_i_j[i-1][0] + y_i_l[2*i-1][fl-1], GRB.GREATER_EQUAL, 0)
                    total_expr.append(kexu_i_j[i-1][0] - y_i_l[2*i-2][fl-1])
                    total_expr.append(kexu_i_j[i-1][0] + y_i_l[2*i-2][fl-1])
                    total_expr.append(eta_i_j[i-1][0] - y_i_l[2*i-1][fl-1])
                    total_expr.append(eta_i_j[i-1][0] + y_i_l[2*i-1][fl-1])
                    p.append(0)
                    p.append(0)
                    p.append(0)
                    p.append(0)
                    if 2*i-2 < k_low and fl == 1:
                        temp_list = [0 for i in range(k_low)]
                        temp_list[2*i-2] = -1
                        P.append(temp_list)
                        temp_list = [0 for i in range(k_low)]
                        temp_list[2*i-2] = 1
                        P.append(temp_list)
                        temp = [0 for i in range(max_Q)]
                        #  是否 + j
                        temp[(fl-1)*(2**theta) * (theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                        Q.append(temp)
                        Q.append(temp)
                    else:
                        P.append([0 for i in range(k_low)])
                        P.append([0 for i in range(k_low)])
                        temp = [0 for i in range(max_Q)]
                        temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                        temp[2*theta* (2**theta) *(theta*multip_number+1) + (2*i-2)* (theta+1) + fl -1] = -1
                        Q.append(temp)
                        temp = [0 for i in range(max_Q)]
                        temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                        temp[2*theta* (2**theta) *(theta*multip_number+1) + (2*i-2)* (theta+1) + fl -1] = 1
                        Q.append(temp)

                    if 2*i-1 < k_low and fl == 1:
                        temp_list = [0 for i in range(k_low)]
                        temp_list[2*i-1] = -1
                        P.append(temp_list)
                        temp_list = [0 for i in range(k_low)]
                        temp_list[2*i-1] = 1
                        P.append(temp_list)
                        temp = [0 for i in range(max_Q)]
                        temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                        Q.append(temp)
                        Q.append(temp)
                    else:
                        P.append([0 for i in range(k_low)])
                        P.append([0 for i in range(k_low)])
                        temp = [0 for i in range(max_Q)]
                        temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1)+ (i-1)*(theta*multip_number+1) + j] = 1
                        temp[2*theta* (2**theta) *(theta*multip_number+1) + (2*i-1)* (theta+1) + fl - 1] = -1
                        Q.append(temp)
                        temp = [0 for i in range(max_Q)]
                        temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                        temp[2*theta* (2**theta) *(theta*multip_number+1) + (2*i-1)* (theta+1) + fl -1] = 1
                        Q.append(temp)
                        

                else:
                    cos = math.cos(math.pi / (2**(j+1)))
                    sin = math.sin(math.pi / (2**(j+1)))
                    model.addLConstr(kexu_i_j[i-1][j] - cos*kexu_i_j[i-1][j-1] - sin*eta_i_j[i-1][j-1], GRB.EQUAL, 0)
                    model.addLConstr(eta_i_j[i-1][j] + sin*kexu_i_j[i-1][j-1] - cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
                    model.addLConstr(eta_i_j[i-1][j] - sin*kexu_i_j[i-1][j-1] + cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
                    total_expr.append(kexu_i_j[i-1][j] - cos*kexu_i_j[i-1][j-1] - sin*eta_i_j[i-1][j-1])
                    total_expr.append(-kexu_i_j[i-1][j] + cos*kexu_i_j[i-1][j-1] + sin*eta_i_j[i-1][j-1])
                    total_expr.append(eta_i_j[i-1][j] + sin*kexu_i_j[i-1][j-1] - cos*eta_i_j[i-1][j-1])
                    total_expr.append(eta_i_j[i-1][j] - sin*kexu_i_j[i-1][j-1] + cos*eta_i_j[i-1][j-1])
                    P.append([0 for i in range(k_low)])
                    P.append([0 for i in range(k_low)])
                    P.append([0 for i in range(k_low)])
                    P.append([0 for i in range(k_low)])
                    p.append(0)
                    p.append(0)
                    p.append(0)
                    p.append(0)
                    temp = [0 for i in range(max_Q)]
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = -cos
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = -sin
                    Q.append(temp)
                    temp = [0 for i in range(max_Q)]
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = -1
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = cos
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = sin
                    Q.append(temp)
                    temp = [0 for i in range(max_Q)]
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = sin
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = -cos
                    Q.append(temp)
                    temp = [0 for i in range(max_Q)]
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                    temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = -sin
                    temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j-1] = cos
                    Q.append(temp)
                    if j == v_l:
                        model.addLConstr(kexu_i_j[i-1][j] - y_i_l[i-1][fl], GRB.LESS_EQUAL, 0)
                        model.addLConstr(eta_i_j[i-1][j] - math.tan(math.pi / (2**(j+1))) * kexu_i_j[i-1][j], GRB.LESS_EQUAL, 0)
                        #  + -meiyou bian 
                        total_expr.append(-kexu_i_j[i-1][j] + y_i_l[i-1][fl])
                        total_expr.append(-eta_i_j[i-1][j] + math.tan(math.pi / (2**(j+1))) * kexu_i_j[i-1][j])
                        P.append([0 for i in range(k_low)])
                        P.append([0 for i in range(k_low)])
                        
                        if fl == theta and i == 1:
                            p.append(1)
                            p.append(0)
                            temp = [0 for i in range(max_Q)]
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = -1                                
                            Q.append(temp)
                            temp = [0 for i in range(max_Q)]
                            temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = -1
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = math.tan(math.pi / (2**(j+1)))
                            Q.append(temp)
                        else:
                            p.append(0)
                            p.append(0)
                            temp = [0 for i in range(max_Q)]
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = -1
                            temp[2*theta* (2**theta) *(theta*multip_number+1) + (i-1)* (theta+1) + fl ] = 1
                        
                            Q.append(temp)
                            temp = [0 for i in range(max_Q)]
                            temp[theta* (2**theta) *(theta*multip_number+1) + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = -1
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = math.tan(math.pi / (2**(j+1)))
                            Q.append(temp)

    for i in range(max_Q):
        temp = gp.LinExpr()
        for j in range(v_number):
            temp += v[j] * Q[j][i]
        model.addLConstr(temp, GRB.EQUAL, 0)

    for i in range(k_low):
        temp = gp.LinExpr()
        for j in range(v_number):
            temp += v[j] * P[j][i]
        model.addLConstr(temp+z_i[i], GRB.EQUAL, 0)



    e_vals,e_vecs = np.linalg.eig(A[m])
    M_bar = 0
    item_one = abs(-min(e_vals)) * N * N
    item_two = abs(b[m] + 1) / 2

    norm_C = 0
    for item in C[m]:
        norm_C += item*item
    
    item_three = math.sqrt(norm_C) * N

    M_bar = item_one + (item_two+item_three) * (item_two+item_three)
    M_bar = math.sqrt(M_bar)

    for i in range(v_number):
        temp += v[i] * p[i] 
    model.addConstr(temp, GRB.LESS_EQUAL, M_bar)

    total_expr_max = []
    model.setParam("OutputFlag", 0)
    for i in range(v_number):
        model.setObjective(total_expr[i], GRB.MAXIMIZE)
        try:
            model.optimize()
        except gp.GurobiError:
            print("Optimize failed due to non-convexity")
            # Solve bilinear model
            model.params.NonConvex = 2
            model.optimize()

        if model.status == GRB.OPTIMAL:
            total_expr_max.append(model.getObjective().getValue())
        else:
            total_expr_max.append(float("inf"))

    # print(total_expr_max)
    total_expr_max_later = []
    for item in total_expr_max:
        if abs(item) < 0.0000001:
            total_expr_max_later.append(0)
        else:
            total_expr_max_later.append(item)
    # print(total_expr_max_later)


    

    v_max = []
    for i in range(v_number):
        model.setObjective(v[i], GRB.MAXIMIZE)
        try:
            model.optimize()
        except gp.GurobiError:
            print("Optimize failed due to non-convexity")
            # Solve bilinear model
            model.params.NonConvex = 2
            model.optimize()

        if model.status == GRB.OPTIMAL:
            v_max.append(model.getObjective().getValue())
        else:
            v_max.append(float("inf"))
    # print(v_max)
    
    v_max_later = []
    for item in v_max:
        if abs(item) < 0.0000001:
            v_max_later.append(0)
        else:
            v_max_later.append(item)
    # print(v_max_later)
    
    mins_index = []
    for i in range(v_number):
        if v_max_later[i] == float("inf"):
            mins_index.append(i)
    # print(mins_index)
    
    # 两个v相减的 界限
    mins_v_upp = []
    # mins_v_low = []
    for i in range(0, len(mins_index), 2):
        first_index = mins_index[i]
        second_index = mins_index[i+1]
        model.setObjective(v[first_index]-v[second_index], GRB.MAXIMIZE)
        try:
            model.optimize()
        except gp.GurobiError:
            print("Optimize failed due to non-convexity")
            # Solve bilinear model
            model.params.NonConvex = 2
            model.optimize()

        if model.status == GRB.OPTIMAL:
            mins_v_upp.append(model.getObjective().getValue())
            v_max_later[first_index] = model.getObjective().getValue()
            v_max_later[second_index] = 0
        else:
            mins_v_upp.append(float("inf"))
            v_max_later[first_index] = float("inf")
            v_max_later[second_index] = 0
    
        # model.setObjective(v[first_index]-v[second_index], GRB.MINIMIZE)
        # try:
        #     model.optimize()
        # except gp.GurobiError:
        #     print("Optimize failed due to non-convexity")
        #     # Solve bilinear model
        #     model.params.NonConvex = 2
        #     model.optimize()

        # if model.status == GRB.OPTIMAL:
        #     mins_v_low.append(model.getObjective().getValue())
        #     v_max_later[second_index] = abs(-model.getObjective().getValue())
        # else:
        #     mins_v_low.append(float("-inf"))
        #     v_max_later[second_index] = -float("-inf")

    # print(mins_v_upp)
    # print(mins_v_low)
    for i in range(v_number):
        total_expr_max_later[i] = math.ceil(total_expr_max_later[i])
        v_max_later[i] = math.ceil(v_max_later[i])
    
    return total_expr_max_later, v_max_later


