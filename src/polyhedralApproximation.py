# This file uses gurobi to solve the Y,Z decomposition qcqp !(2)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np
import math

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, all_Y, all_Z, C, b, u, epslong):
    M = M + 1
    A.append(np.identity(N).tolist())
    all_Y.append(np.identity(N).tolist())
    all_Z.append(np.zeros((N,0)).tolist())
    C.append([0 for i in range(N)])
    b.append(u)

    model = gp.Model()
    # Variables
    vars = [model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(N)]
       
    # Objective constraint
    # obj = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    

    temp = gp.QuadExpr()
    for i in range(N):
        for j in range(N):
            temp += A[0][i][j]*vars[i]*vars[j]
    for i in range(N):
        temp += 2*C[0][i]*vars[i]

    # model.addConstr(temp, GRB.LESS_EQUAL, obj)

    model.setObjective(temp, GRB.MINIMIZE)

    epslong = float(epslong)
    multip_number = math.ceil(math.log(2/epslong))
    t = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(M)]

    for m in range(1,M):
        Y = all_Y[m]
        Z = all_Z[m]
        YT = np.transpose(Y)
        ZT = np.transpose(Z)
        k_upper = len(Y[0]) +1 
        k_low = len(Z[0]) +1 

        y_i = []
        z_i = []
        if m != 0:
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

        # temp = gp.QuadExpr()
        # for i in range(k_upper):
        #     temp += y_i[i] * y_i[i]
        # model.addConstr(temp, GRB.LESS_EQUAL, t[m]*t[m])
        
        # temp = gp.QuadExpr()
        # for i in range(k_low):
        #     temp += z_i[i] * z_i[i]
        # model.addConstr(temp, GRB.GREATER_EQUAL, t[m]*t[m])

        theta = 0
        for i in range(999):
            if 2**i >= k_upper:
                theta = i
                break
        if theta == 0:
            theta = 1
        
        # y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for l in range(theta+1)] for i in range(2**theta)]   
        # for i in range(k_upper):
        #     y_i_l[i][0] = y_i[i]
        
        # y_i_l[0][theta] = t[m]
        y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(2**(theta-l)+1)] for l in range(theta+1)]   
        for i in range(k_upper):
            y_i_l[0][i+1] = y_i[i]
        
        y_i_l[theta][1] = t[m]
            
        for fl in range(1,theta+1):
            v_l = math.floor(fl*multip_number)
            kexu_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl)+1)] 
            eta_i_j =  [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl)+1)] 
            
            for i in range(1, 2**(theta-fl) + 1):
                for j in range(v_l+1):

                    if j == 0:
                        model.addLConstr(kexu_i_j[i][0] - y_i_l[fl-1][2*i-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(kexu_i_j[i][0] + y_i_l[fl-1][2*i-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][0] - y_i_l[fl-1][2*i], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][0] + y_i_l[fl-1][2*i], GRB.GREATER_EQUAL, 0)
                    else:
                        cos = math.cos(math.pi / (2**(j+1)))
                        sin = math.sin(math.pi / (2**(j+1)))
                        model.addLConstr(kexu_i_j[i][j] - cos*kexu_i_j[i][j-1] - sin*eta_i_j[i][j-1], GRB.EQUAL, 0)
                        model.addLConstr(eta_i_j[i][j] + sin*kexu_i_j[i][j-1] - cos*eta_i_j[i][j-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][j] - sin*kexu_i_j[i][j-1] + cos*eta_i_j[i][j-1], GRB.GREATER_EQUAL, 0)
                        if j == v_l:
                            model.addLConstr(kexu_i_j[i][j] - y_i_l[fl][i], GRB.LESS_EQUAL, 0)
                            model.addLConstr(eta_i_j[i][j] - math.tan(math.pi / (2**(j+1))) * kexu_i_j[i][j], GRB.LESS_EQUAL, 0)

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
        
        y_i_l = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for l in range(2**(theta-l)+1)] for l in range(theta+1)]   
        for i in range(k_low):
            y_i_l[0][i+1] = w_i[i]
        
        y_i_l[theta][1] = 1

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
        eta_begin = theta* (2**theta) *(theta*multip_number+1)
        y_begin = 2*theta* (2**theta) *(theta*multip_number+1)
        # max_Q = math.floor(2*theta* (2**(theta-1)) *(theta*multip_number+1) + (2**theta) * (theta+1))
        Q = []

        for fl in range(1,theta+1):
            v_l = math.floor(fl*multip_number)
            kexu_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl)+1)] 
            eta_i_j = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(v_l + 1)] for i in range(2**(theta-fl)+1)] 
            
            for i in range(1, 2**(theta-fl) + 1):
                for j in range(v_l+1):
                    if j == 0:
                        model.addLConstr(kexu_i_j[i][0] - y_i_l[fl-1][2*i-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(kexu_i_j[i][0] + y_i_l[fl-1][2*i-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][0] - y_i_l[fl-1][2*i], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][0] + y_i_l[fl-1][2*i], GRB.GREATER_EQUAL, 0)
                        total_expr.append(kexu_i_j[i][0] - y_i_l[fl-1][2*i-1])
                        total_expr.append(kexu_i_j[i][0] + y_i_l[fl-1][2*i-1])
                        total_expr.append(eta_i_j[i][0] - y_i_l[fl-1][2*i])
                        total_expr.append(eta_i_j[i][0] + y_i_l[fl-1][2*i])
                        p.append(0)
                        p.append(0)
                        p.append(0)
                        p.append(0)
                        if fl == 1:
                            temp_list = [0 for i in range(2**theta+1)]
                            temp_list[2*i-1] = -1
                            P.append(temp_list)
                            temp_list = [0 for i in range(2**theta+1)]
                            temp_list[2*i-1] = 1
                            P.append(temp_list)
                            temp_list = [0 for i in range(2**theta+1)]
                            temp_list[2*i] = -1
                            P.append(temp_list)
                            temp_list = [0 for i in range(2**theta+1)]
                            temp_list[2*i] = 1
                            P.append(temp_list)
 
                            temp = [0 for i in range(max_Q)]
                            temp[(fl-1)*(2**theta) * (theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                            Q.append(temp)
                            Q.append(temp)
                            temp = [0 for i in range(max_Q)]
                            temp[eta_begin+ (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                            Q.append(temp)
                            Q.append(temp)
                        else:
                            P.append([0 for i in range(2**theta+1)])
                            P.append([0 for i in range(2**theta+1)])
                            P.append([0 for i in range(2**theta+1)])
                            P.append([0 for i in range(2**theta+1)])
                            temp = [0 for i in range(max_Q)]
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                            temp[y_begin + (2*i-2)* (theta+1) + fl -1] = -1
                            Q.append(temp)
                            temp = [0 for i in range(max_Q)]
                            temp[(fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1)] = 1
                            temp[y_begin + (2*i-2)* (theta+1) + fl -1] = 1
                            Q.append(temp)

                            temp = [0 for i in range(max_Q)]
                            temp[eta_begin + (fl-1)*(2**theta) *(theta*multip_number+1)+ (i-1)*(theta*multip_number+1) + j] = 1
                            temp[y_begin + (2*i-1)* (theta+1) + fl - 1] = -1
                            Q.append(temp)
                            temp = [0 for i in range(max_Q)]
                            temp[eta_begin + (fl-1)*(2**theta) *(theta*multip_number+1) + (i-1)*(theta*multip_number+1) + j] = 1
                            temp[y_begin + (2*i-1)* (theta+1) + fl -1] = 1
                            Q.append(temp)
                            

                    else:
                        cos = math.cos(math.pi / (2**(j+1)))
                        sin = math.sin(math.pi / (2**(j+1)))
                        model.addLConstr(kexu_i_j[i][j] - cos*kexu_i_j[i][j-1] - sin*eta_i_j[i][j-1], GRB.EQUAL, 0)
                        model.addLConstr(eta_i_j[i][j] + sin*kexu_i_j[i][j-1] - cos*eta_i_j[i][j-1], GRB.GREATER_EQUAL, 0)
                        model.addLConstr(eta_i_j[i][j] - sin*kexu_i_j[i][j-1] + cos*eta_i_j[i][j-1], GRB.GREATER_EQUAL, 0)
                        total_expr.append(kexu_i_j[i][j] - cos*kexu_i_j[i][j-1] - sin*eta_i_j[i][j-1])
                        total_expr.append(-kexu_i_j[i][j] + cos*kexu_i_j[i][j-1] + sin*eta_i_j[i][j-1])
                        total_expr.append(eta_i_j[i][j] + sin*kexu_i_j[i][j-1] - cos*eta_i_j[i][j-1])
                        total_expr.append(eta_i_j[i][j] - sin*kexu_i_j[i][j-1] + cos*eta_i_j[i][j-1])
                        P.append([0 for i in range(2**theta+1)])
                        P.append([0 for i in range(2**theta+1)])
                        P.append([0 for i in range(2**theta+1)])
                        P.append([0 for i in range(2**theta+1)])
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
                            model.addLConstr(kexu_i_j[i][j] - y_i_l[fl][i], GRB.LESS_EQUAL, 0)
                            model.addLConstr(eta_i_j[i][j] - math.tan(math.pi / (2**(j+1))) * kexu_i_j[i][j], GRB.LESS_EQUAL, 0)
                            total_expr.append(-kexu_i_j[i][j] + y_i_l[fl][i])
                            total_expr.append(-eta_i_j[i][j] + math.tan(math.pi / (2**(j+1))) * kexu_i_j[i][j])
                            P.append([0 for i in range(2**theta+1)])
                            P.append([0 for i in range(2**theta+1)])
                            
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
                temp += v[j] * P[j][i+1]
            model.addLConstr(temp+z_i[i], GRB.EQUAL, 0)

        temp = gp.LinExpr()
        for i in range(v_number):
            temp += p[i] * v[i]
        model.addLConstr(temp-t[m], GRB.GREATER_EQUAL, 0)   



        # total_expr_max_later, v_max_later = v_total_max(N, M, obj_Q, obj_c, A, all_Y, all_Z, C, b, u, epslong, m, obj_max)

        i_i = [model.addVar(vtype=GRB.BINARY) for i in range(v_number)]
        for i in range(v_number):
            # model.addConstr(v[i] * total_expr[i], GRB.EQUAL, 0)
            model.addConstr(v[i]-9999*i_i[i], GRB.LESS_EQUAL, 0)
            model.addConstr(total_expr[i]-9999*(1-i_i[i]), GRB.LESS_EQUAL, 0)



    model.setParam('TimeLimit', 60*60*3) 
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
        return
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return
    



if __name__ == "__main__":
    with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["Y"], load_dict["Z"], load_dict["C"], load_dict["b"], load_dict["u"], sys.argv[2])

    