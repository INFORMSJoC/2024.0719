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
                    model.addConstr( math.cos(math.pi*index_mul / (2**(multip_number+1))) * xi_0 + math.sin(math.pi*index_mul / (2**(multip_number+1)))*eta_0<= temp_y3)
                    # tempConstr = 
                    # model.update()
                    # tempConstr.Lazy = 2
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
                        model.addConstr( math.cos(math.pi*index_mul / (2**(multip_number+1))) * xi_0 + math.sin(math.pi*index_mul / (2**(multip_number+1)))*eta_0<= temp_y3)
                    
                else:
                    layers_variable[lay][i] = layers_variable[lay-1][2*i]

        

        # Z部分
        if k_low == 1:
            model.addConstr( z_i[0] >= t[m] )
            continue

        theta = 0
        for i in range(999):
            if 2**i >= k_low:
                theta = i
                break

        w_i = [model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(2**theta)] 
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
            temp_y1 = w_i[2*i]
            temp_y2 = w_i[2*i+1]
            temp_y3 = layers_variable[0][i]
            # j 就是index_mul  v是multip_number
            for index_mul in range(0, 2**multip_number+1):
                model.addConstr( math.cos(math.pi*index_mul / (2**(multip_number+1))) * temp_y1 + math.sin(math.pi*index_mul / (2**(multip_number+1)))*temp_y2<= temp_y3)
                total_expr[0][i].append(-math.cos(math.pi*index_mul / (2**(multip_number+1))) * temp_y1 - math.sin(math.pi*index_mul / (2**(multip_number+1)))*temp_y2 + temp_y3)
                v[0][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
                i_i[0][i].append(model.addVar(vtype=GRB.BINARY))
        
        for lay in range(1, layers):
            for i in range(layers_count[lay]):
                temp_y1 = layers_variable[lay-1][2*i]
                temp_y2 = layers_variable[lay-1][2*i+1]
                temp_y3 = layers_variable[lay][i]
                # j 就是index_mul  v是multip_number
                for index_mul in range(0, 2**multip_number+1):
                    model.addConstr( math.cos(math.pi*index_mul / (2**(multip_number+1))) * temp_y1 + math.sin(math.pi*index_mul / (2**(multip_number+1)))*temp_y2<= temp_y3)
                    total_expr[lay][i].append(-math.cos(math.pi*index_mul / (2**(multip_number+1))) * temp_y1 - math.sin(math.pi*index_mul / (2**(multip_number+1)))*temp_y2+ temp_y3)
                    v[lay][i].append(model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS))
                    i_i[lay][i].append(model.addVar(vtype=GRB.BINARY))

        # TODO:
        M_bar = 9999 #9999 # GRB.INFINITY

        for i in range(layers):
            for j in range(layers_count[i]):
                bound_constr = gp.LinExpr()
                for k in range(2**multip_number+1):
                    # bound_constr += i_i[i][j][k]
                    model.addConstr(v[i][j][k] <= M_bar*i_i[i][j][k])
                    model.addConstr(total_expr[i][j][k] <= M_bar*( 1-i_i[i][j][k] ))
                # model.addConstr(bound_constr <= 2)

        aT = []
        for j in range(0, 2**multip_number+1):
            aT.append(math.cos(math.pi*j / (2**(multip_number+1))))
        bT = []
        for j in range(0, 2**multip_number+1):
            bT.append(math.sin( math.pi*j / (2**(multip_number+1)) ))
        
        # 对于第一层
        for index in range(layers_count[0]):
            # constraints_number += 2
            aTv = gp.LinExpr()
            bTv = gp.LinExpr()
            for temp_i in range(0, 2**multip_number+1):
                aTv += aT[temp_i]*v[0][index][temp_i]
                bTv += bT[temp_i]*v[0][index][temp_i]
            
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


        
        for layer_index in range(layers-1):
            for index in range(layers_count[layer_index]):
                # constraints_number += 1
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

    