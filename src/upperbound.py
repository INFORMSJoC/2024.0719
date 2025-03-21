#!/usr/bin/env python3
import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json


# N: varible count, M: constrain count, A: all martric: m*N*N, C: m* N, b m
def qcqp_approximation(N, M, Q, c, A, all_Y, all_Z, C, b, u, epslong):
    M = M + 1
    all_Y.append(np.identity(N).tolist())
    all_Z.append([[] for i in range(N)])
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    M = M + 1
    A.append(Q)
    C.append(c)
    b.append(0)
    all_Y.append(u[0])
    all_Z.append(u[1])
    

    # N = N + 1
    # for i in range(M):
    #     for row in A[i]:
    #         row.append(0)
    #     A[i].append([0 for i in range(N)])

    #     C[i].append(0)
    # C[-1][-1] = -0.5

    # for i in range(M):
    #     all_Y[i].append([0 for k in range(len(all_Y[i][0]))])
    #     all_Z[i].append([0 for k in range(len(all_Z[i][0]))])

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
        
        # 先求层数：
        layers = 0
        layers_count = [] # 每一层的数目
        layers_variable = [] # 每一层的变量
        current_layer_count = k_upper
        while current_layer_count != 1:
            layers += 1
            current_layer_count = current_layer_count // 2 + current_layer_count%2
            layers_count.append(current_layer_count)
            layers_variable.append([model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(current_layer_count)])

        # print(layers)
        # print(layers_count)
        # print(layers_variable)
        

        layers_variable[layers-1][0] = t[m]
        for i in range(layers_count[0]):
            if 2*i+1 < k_upper:
                # model.addConstr(y_i[2*i]*y_i[2*i] + y_i[2*i+1]*y_i[2*i+1] - layers_variable[0][i]*layers_variable[0][i], GRB.LESS_EQUAL, 0)
                temp_y1 = y_i[2*i]
                temp_y2 = y_i[2*i+1]
                temp_y3 = layers_variable[0][i]
                xi_0 = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
                eta_0 = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
                model.addConstr(xi_0 >= -temp_y1)
                model.addConstr(xi_0 >= temp_y1)
                model.addConstr(eta_0 >= -temp_y2)
                model.addConstr(eta_0 >= temp_y2)
                # j 就是index_mul  v是multip_number
                for index_mul in range(0, 2**multip_number+1):
                    xi_0_const = math.sin(math.pi*(index_mul+1) / (2**(multip_number+1))) - math.sin(math.pi*index_mul / (2**(multip_number+1)))
                    eta_0_const = math.cos(math.pi*index_mul / (2**(multip_number+1))) - math.cos(math.pi*(index_mul+1) / (2**(multip_number+1)))
                    y3_const = math.sin(math.pi / (2**(multip_number+1)))
                    model.addConstr( (xi_0_const/y3_const) * xi_0 + (eta_0_const/y3_const)*eta_0<= temp_y3)
            else:
                layers_variable[0][i] = y_i[2*i]
        
        for lay in range(1, layers):
            for i in range(layers_count[lay]):
                if 2*i+1 < layers_count[lay-1]:
                    # model.addConstr(layers_variable[lay-1][2*i]*layers_variable[lay-1][2*i] + layers_variable[lay-1][2*i+1]*layers_variable[lay-1][2*i+1] - layers_variable[lay][i]*layers_variable[lay][i], GRB.LESS_EQUAL, 0)
                    temp_y1 = layers_variable[lay-1][2*i]
                    temp_y2 = layers_variable[lay-1][2*i+1]
                    temp_y3 = layers_variable[lay][i]
                    xi_0 = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
                    eta_0 = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
                    model.addConstr(xi_0 >= -temp_y1)
                    model.addConstr(xi_0 >= temp_y1)
                    model.addConstr(eta_0 >= -temp_y2)
                    model.addConstr(eta_0 >= temp_y2)
                    # j 就是index_mul  v是multip_number
                    for index_mul in range(0, 2**multip_number+1):
                        xi_0_const = math.sin(math.pi*(index_mul+1) / (2**(multip_number+1))) - math.sin(math.pi*index_mul / (2**(multip_number+1)))
                        eta_0_const = math.cos(math.pi*index_mul / (2**(multip_number+1))) - math.cos(math.pi*(index_mul+1) / (2**(multip_number+1)))
                        y3_const = math.sin(math.pi / (2**(multip_number+1)))
                        model.addConstr( (xi_0_const/y3_const) * xi_0 + (eta_0_const/y3_const)*eta_0<= temp_y3)
                    
                else:
                    layers_variable[lay][i] = layers_variable[lay-1][2*i]

        
        
        
        # temp = gp.QuadExpr()
        # for i in range(k_upper):
        #     temp += y_i[i] * y_i[i]
        # model.addConstr(temp-t[m]*t[m], GRB.LESS_EQUAL, 0)
        
        # temp = gp.QuadExpr()
        # for i in range(k_low):
        #     temp += z_i[i] * z_i[i]
        # model.addConstr(temp-t[m]*t[m], GRB.GREATER_EQUAL, 0)

        # Z部分
        if k_low == 1:
            model.addConstr( z_i[0] >= t[m] )
            continue

        theta = 0
        for i in range(999):
            if 2**i >= k_low:
                theta = i
                break

        w_i = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(2**theta)] # 每一项大于0
        for i in range(k_low, 2**theta):
            z_i.append(0)
        

        # 先求层数：
        layers = 0
        layers_count = [] # 每一 成的数目
        layers_variable = [] # 每一程的变量
        current_layer_count = 2**theta
        while current_layer_count != 1:
            layers += 1
            current_layer_count = current_layer_count // 2 + current_layer_count%2
            layers_count.append(current_layer_count)
            layers_variable.append([model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(current_layer_count)])
        # print(layers)
        # print(layers_count)
        # print(layers_variable)

        total_expr = [
            [ [] for j in range(layers_count[i]) ] for i in range(layers)
        ]
        v = [
            [ [] for j in range(layers_count[i]) ] for i in range(layers)
        ]
        i_i = [
            [ [] for j in range(layers_count[i]) ] for i in range(layers)
        ]

        layers_variable[layers-1][0] = 1
        # 第一层
        for i in range(layers_count[0]):
            # if 2*i+1 < k_low:
            temp_y1 = w_i[2*i]
            temp_y2 = w_i[2*i+1]
            temp_y3 = layers_variable[0][i]
            # j 就是index_mul  v是multip_number
            for index_mul in range(0, 2**multip_number+1):
                xi_0_const = math.sin(math.pi*(index_mul+1) / (2**(multip_number+1))) - math.sin(math.pi*index_mul / (2**(multip_number+1)))
                eta_0_const = math.cos(math.pi*index_mul / (2**(multip_number+1))) - math.cos(math.pi*(index_mul+1) / (2**(multip_number+1)))
                y3_const = math.sin(math.pi / (2**(multip_number+1)))
                model.addConstr( (xi_0_const/y3_const) * temp_y1 + (eta_0_const/y3_const)*temp_y2<= temp_y3)
                total_expr[0][i].append( -(xi_0_const/y3_const) * temp_y1 - (eta_0_const/y3_const)*temp_y2 + temp_y3)
                v[0][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
                i_i[0][i].append(model.addVar(vtype=GRB.BINARY))
            # else:
                # layers_variable[0][i] = w_i[2*i]
        
        for lay in range(1, layers):
            for i in range(layers_count[lay]):
                # if 2*i+1 < layers_count[lay-1]:
                temp_y1 = layers_variable[lay-1][2*i]
                temp_y2 = layers_variable[lay-1][2*i+1]
                temp_y3 = layers_variable[lay][i]
                # j 就是index_mul  v是multip_number
                for index_mul in range(0, 2**multip_number+1):
                    xi_0_const = math.sin(math.pi*(index_mul+1) / (2**(multip_number+1))) - math.sin(math.pi*index_mul / (2**(multip_number+1)))
                    eta_0_const = math.cos(math.pi*index_mul / (2**(multip_number+1))) - math.cos(math.pi*(index_mul+1) / (2**(multip_number+1)))
                    y3_const = math.sin(math.pi / (2**(multip_number+1)))
                    model.addConstr( (xi_0_const/y3_const) * temp_y1 + (eta_0_const/y3_const)*temp_y2<= temp_y3)
                    total_expr[lay][i].append(-(xi_0_const/y3_const) * temp_y1 - (eta_0_const/y3_const)*temp_y2+ temp_y3)
                    v[lay][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
                    i_i[lay][i].append(model.addVar(lb=0, ub=1, vtype=GRB.BINARY))
                # else:
                    # layers_variable[lay][i] = layers_variable[lay-1][2*i]

        # TODO:
        M_bar = 99999 #9999 # GRB.INFINITY
        for i in range(layers):
            for j in range(layers_count[i]):
                bound_constr = gp.LinExpr()
                for k in range(2**multip_number+1):
                    # bound_constr += i_i[i][j][k]
                    model.addConstr(v[i][j][k]-M_bar*i_i[i][j][k], GRB.LESS_EQUAL, 0)
                    model.addConstr(total_expr[i][j][k]-M_bar*(1-i_i[i][j][k]), GRB.LESS_EQUAL, 0)
                # model.addConstr(bound_constr <= 2)
        
        aT = []        
        for index_mul in range(0, 2**multip_number+1):
            xi_const = math.sin(math.pi*(index_mul+1) / (2**(multip_number+1))) - math.sin(math.pi*index_mul / (2**(multip_number+1)))
            y_const = math.sin(math.pi / (2**(multip_number+1)))
            aT.append(xi_const/y_const)

        bT = []
        for index_mul in range(0, 2**multip_number+1):
            eta_const = math.cos(math.pi*index_mul / (2**(multip_number+1))) - math.cos(math.pi*(index_mul+1) / (2**(multip_number+1)))
            y_const = math.sin(math.pi / (2**(multip_number+1)))
            bT.append(eta_const/y_const)

        
        
        # 对于第一层
        for index in range(layers_count[0]):
            # if(len(v[0][index]) > 0):
            aTv = gp.LinExpr()
            bTv = gp.LinExpr()
            for temp_i in range(0, 2**multip_number+1):
                aTv += aT[temp_i]*v[0][index][temp_i]
                bTv += bT[temp_i]*v[0][index][temp_i]
            # if 2*index < k_low:
            model.addConstr(aTv >= -z_i[2*index])
            model.addConstr(aTv >= z_i[2*index])
            temp_var1 = model.addVar(vtype=GRB.BINARY)
            abs_zi = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
            model.addConstr(abs_zi + z_i[2*index] <= M_bar*temp_var1)
            model.addConstr(abs_zi + z_i[2*index] >= -M_bar*temp_var1)
            model.addConstr(abs_zi - z_i[2*index] <= M_bar*(1-temp_var1))
            model.addConstr(abs_zi - z_i[2*index] >= -M_bar*(1-temp_var1))

            temp_var2 = model.addVar(vtype=GRB.BINARY)
            model.addConstr(aTv - abs_zi <= M_bar*temp_var2)
            model.addConstr(w_i[2*index] <= M_bar*(1-temp_var2))
            # if 2*index+1 < k_low:
            model.addConstr(bTv >= -z_i[2*index+1])
            model.addConstr(bTv >= z_i[2*index+1]) 
            temp_var1 = model.addVar(vtype=GRB.BINARY)
            abs_zi = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
            model.addConstr(abs_zi + z_i[2*index+1] <= M_bar*temp_var1)
            model.addConstr(abs_zi + z_i[2*index+1] >= -M_bar*temp_var1)
            model.addConstr(abs_zi - z_i[2*index+1] <= M_bar*(1-temp_var1))
            model.addConstr(abs_zi - z_i[2*index+1] >= -M_bar*(1-temp_var1))

            temp_var2 = model.addVar(vtype=GRB.BINARY)
            model.addConstr(bTv - abs_zi <= M_bar*temp_var2)
            model.addConstr(w_i[2*index+1] <= M_bar*(1-temp_var2))
            
            
        for layer_index in range(0, layers-1):
            for index in range(layers_count[layer_index]):
                # if len(v[layer_index][index]) > 0 and len(v[layer_index+1][index//2]) > 0:
                eTv = gp.LinExpr()
                for temp_i in range(2**multip_number+1):
                    eTv += -v[layer_index][index][temp_i]
                cTv = gp.LinExpr()
                # 46 
                if index % 2 == 0: # 0 2 4 6 
                    if (index+1+1)//2 % 2 != 0: # 奇数 0  4 
                        for temp_i in range(2**multip_number+1):
                            cTv += aT[temp_i]*v[layer_index+1][index//2][temp_i]
                    else: # 2 6
                        for temp_i in range(2**multip_number+1):
                            cTv += aT[temp_i]*v[layer_index+1][index//2][temp_i]
                # 47
                else: # 1 3 5 7
                    if (index+1)//2 % 2 != 0: # 奇数 1 5
                        for temp_i in range(2**multip_number+1):
                            cTv += bT[temp_i]*v[layer_index+1][(index-1)//2][temp_i]
                    else:# 3 7
                        for temp_i in range(2**multip_number+1):
                            cTv += bT[temp_i]*v[layer_index+1][(index-1)//2][temp_i]
                # KKT 以及 本来的约束
                model.addConstr(eTv+cTv >= 0)

                temp_var = model.addVar(vtype=GRB.BINARY)
                model.addConstr(eTv+cTv <= M_bar*temp_var)
                model.addConstr(layers_variable[layer_index][index] <= M_bar*(1-temp_var))
                
        # 最后一层
        min_eTv = gp.LinExpr()
        for i in range(2**multip_number+1):
            min_eTv += v[-1][0][i]        
        model.addConstr(min_eTv - t[m]>=0)


    model.setObjective(obj, GRB.MINIMIZE)
    model.setParam('TimeLimit', 60*60*3) # s
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
        for i in range(N-1):
            for j in range(N-1):
                real_obj += Q[i][j]*solution[i]*solution[j]
        for i in range(N-1):
            real_obj += 2*c[i]*solution[i]
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




