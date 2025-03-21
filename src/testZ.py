#!/usr/bin/env python3
import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json

def qcqp_approximation():
    model = gp.Model()
    epslong = 0.01
    multip_number = math.ceil(math.log(2/epslong))


        
    z_i = [10,-16,13,-14,7,-12,-20]
    k_low = len(z_i)

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

    
    i_i = [model.addVar(vtype=GRB.BINARY) for i in range(v_number)]
    for i in range(v_number):
        # model.addConstr(v[i]*total_expr[i] == 0)
        model.addConstr(v[i]-99999*i_i[i], GRB.LESS_EQUAL, 0)
        model.addConstr(total_expr[i]-99999*(1-i_i[i]), GRB.LESS_EQUAL, 0)

    min_v = gp.LinExpr()
    for i in range(v_number):
        min_v += p[i] * v[i]
    
    max_w = gp.LinExpr()
    should_solution = 0
    for i in range(k_low):
        max_w += w_i[i]*z_i[i]
        should_solution += z_i[i]*z_i[i]
    # model.addLConstr(temp-t[m], GRB.GREATER_EQUAL, 0)   
    
    model.setObjective(0, GRB.MAXIMIZE)
    model.setParam('TimeLimit', 60*60*3) # s
    model.setParam('IntFeasTol', 0.000001)
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        # Solve bilinear model
        model.params.NonConvex = 2
        model.optimize()


    solution = []
    if model.status == GRB.OPTIMAL:
        s = 0
        for i in range(k_low):
            solution.append(w_i[i].X)
            s += w_i[i].X * z_i[i]
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        print(s)
        print("should_solution: {0}".format(math.sqrt(should_solution)))
        
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
    qcqp_approximation()




