#!/usr/bin/env python3
import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json
# from approximate_sub import v_total_max
# from estimateObj import estimate_obj


def qcqp_approximation(N, M, obj_Q, obj_c, A, all_Y, all_Z, C, b, u, epslong):
    # obj_min = estimate_obj(N, deepcopy(obj_Q), deepcopy(obj_c))
    # obj_max = abs(obj_min)

    # bounded
    M = M + 1
    all_Y.append(np.identity(N).tolist())
    all_Z.append([[] for i in range(N)])
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    # objective -> constraint
    M = M + 1
    A.append(obj_Q)
    C.append(obj_c)
    b.append(-6666)
    all_Y.append(u[0])
    all_Z.append(u[1])

    model = gp.Model()
    # Variables
    vars = []
    for i in range(N): 
        vars.append(model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))

    # Objective
    obj = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)

    epslong = float(epslong)
    multip_number = math.ceil(math.log(2/epslong))
    t = []
    for i in range(M):
        t.append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))

    for m in range(M):
        Y = all_Y[m]
        Z = all_Z[m]
        YT = np.transpose(Y)
        ZT = np.transpose(Z)
        k_upper = len(Y[0]) +1 
        k_low = len(Z[0]) +1  
        
        y_i = []
        z_i = []
        if m != M-1:
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
        else:
            for i in range(k_upper):
                temp = gp.LinExpr()
                if i == k_upper-1:
                    for j in range(N):
                        temp += C[m][j]*vars[j]
                    y_i.append(((obj-1) / 2) - temp)
                else:
                    for j in range(N):
                        temp += YT[i][j]*vars[j]
                    y_i.append(temp)

            for i in range(k_low):
                temp = gp.LinExpr()
                if i == k_low-1:
                    for j in range(N):
                        temp += C[m][j]*vars[j]
                    z_i.append(((obj+1) / 2) - temp)
                else:
                    for j in range(N):
                        temp += ZT[i][j]*vars[j]
                    z_i.append(temp)
        
        theta = 0
        for i in range(999):
            if 2**i >= k_upper:
                theta = i
                break
        if theta == 0:
            theta = 1
        
        y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for l in range(theta+1)] for i in range(2**theta)]   
        for i in range(k_upper):
            y_i_l[i][0] = y_i[i]
        
        y_i_l[0][theta] = t[m]
            
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

                    else:
                        cos = math.cos(math.pi / (2**(j+1)))
                        sin = math.sin(math.pi / (2**(j+1)))
                        model.addLConstr(kexu_i_j[i-1][j] - cos*kexu_i_j[i-1][j-1] - sin*eta_i_j[i-1][j-1], GRB.EQUAL, 0)
                        model.addLConstr(eta_i_j[i-1][j] + sin*kexu_i_j[i-1][j-1] - cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i-1][j] - sin*kexu_i_j[i-1][j-1] + cos*eta_i_j[i-1][j-1], GRB.GREATER_EQUAL, 0)
                        
                        if j == v_l:
                            model.addLConstr(kexu_i_j[i-1][j] - y_i_l[i-1][fl], GRB.LESS_EQUAL, 0)
                            model.addLConstr(eta_i_j[i-1][j] - math.tan(math.pi / (2**(j+1))) * kexu_i_j[i-1][j], GRB.LESS_EQUAL, 0)

        # temp = gp.QuadExpr()
        # for i in range(k_low):
        #     temp += z_i[i] * z_i[i]
        # model.addConstr(temp-t[m]*t[m], GRB.GREATER_EQUAL, 0)


        if k_low == 1:
            model.addLConstr(z_i[0]-t[m], GRB.GREATER_EQUAL, 0)
            # model.addLConstr(z_i[0]+t[m], GRB.GREATER_EQUAL, 0)
            continue

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

        temp = gp.LinExpr()
        for i in range(v_number):
            temp += p[i] * v[i]
        model.addLConstr(temp-t[m], GRB.GREATER_EQUAL, 0)   



        # total_expr_max_later, v_max_later = v_total_max(N, M, obj_Q, obj_c, A, all_Y, all_Z, C, b, u, epslong, m, obj_max)

        i_i = [model.addVar(vtype=GRB.BINARY) for i in range(v_number)]
        for i in range(v_number):
            model.addConstr(v[i]-99999*i_i[i], GRB.LESS_EQUAL, 0)
            model.addConstr(total_expr[i]-99999*(1-i_i[i]), GRB.LESS_EQUAL, 0)

    
    
    model.setObjective(obj, GRB.MINIMIZE)
    model.setParam('TimeLimit', 60*60*3) # s
    # model.setParam('IntFeasTol', 0.000001)
    try:
        model.optimize()
    except gp.GurobiError:
        print("Optimize failed due to non-convexity")
        # Solve bilinear model
        model.params.NonConvex = 2
        model.optimize()


    solution = []
    if model.status == GRB.OPTIMAL:
        for i in range(N):
            solution.append(vars[i].X)
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        real_obj = 0
        for i in range(N):
            for j in range(N):
                real_obj += obj_Q[i][j]*solution[i]*solution[j]
        for i in range(N):
            real_obj += 2*obj_c[i]*solution[i]
        print("RealObjective: {0}".format(real_obj))
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
        qcqp_approximation(load_dict["N"], load_dict["M"], load_dict["Q"], load_dict["c"], load_dict["A"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"], sys.argv[2])




