# This file uses gurobi to solve the primary qcqp !(P)!
import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np
from numpy import array 

# N: varible count, M: constrain count, A: all martric, b
def qcqp_optimize(N, M, A, C, b, u):
    M = M + 1
    A.append(np.identity(N).tolist())
    C.append([0 for i in range(N)])
    b.append(N*N)

    model = gp.Model()
    # Variables
    vars = [model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for i in range(N)]
    lda = model.addVar(lb=0, ub=1, vtype=GRB.CONTINUOUS)
       
    # Objective constraint
    # obj = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
    x1 = [array(-0.20108637), array(-0.15743232), array(0.55310504), array(-0.09106077), array(0.05317563), array(0.40720906), array(0.20710049), array(-0.11720683), array(-0.05920927), array(0.16644911), array(0.18038183), array(-0.18484584), array(0.01446304), array(-0.26572161), array(0.27837232), array(0.04992307), array(0.36208691), array(0.18700689), array(0.18251305), array(0.47934106), array(0.15152914), array(-0.34924998), array(0.20789069), array(0.39502403), array(-0.10892558), array(-0.28415899), array(-0.35216895), array(-0.12207156), array(-0.4703647), array(-0.03735765), array(0.52662848), array(0.44721489), array(0.05918634), array(0.45886207), array(-0.25381869), array(0.56200077), array(0.06564324), array(-0.0113564), array(-0.00999413), array(-0.22222492), array(0.00808219), array(-0.34875321), array(0.33781307), array(0.38881774), array(0.1438025), array(0.3274466), array(0.0297796), array(0.10759786), array(-0.23154299), array(-0.02913569), array(0.29961866), array(-0.43757056), array(0.05754073), array(0.32683042), array(0.22813059), array(-0.13974309), array(0.23684328), array(-0.16117564), array(-0.0729356), array(0.46687476), array(0.40310423), array(0.17017103), array(0.05790432), array(0.25208665), array(0.24097284), array(-0.44744155), array(-0.06622656), array(0.34292872), array(-0.30952277), array(-0.23557343), array(-0.06175733), array(-0.30643633), array(-0.5526081), array(0.44024086), array(-0.47982429), array(-0.14103609), array(-0.21093597), array(-0.37597282), array(0.10216862), array(-0.18695763), array(0.22761717), array(-0.3693656), array(0.18432483), array(0.25086707), array(-0.27631159), array(0.16648743), array(0.22018566), array(0.17812744), array(0.02311157), array(-0.1049693), array(0.10351717), array(-0.1406968), array(0.43487982), array(0.22716258), array(-0.2672779), array(0.02404673), array(0.52588471), array(-0.24702702), array(-0.086285), array(-0.47119984)]
    x1 = [float(item) for item in x1]
    x2 = [array(-0.01871499), array(-0.08787407), array(0.23381703), array(-0.11891815), array(-0.0492816), array(0.16958659), array(0.16061666), array(-0.05510649), array(-0.06087119), array(-0.02948276), array(0.05591711), array(-0.06555785), array(-0.00907664), array(-0.05002822), array(0.11503377), array(-0.02225029), array(0.14770534), array(0.15228119), array(0.07149074), array(0.15623031), array(0.06252894), array(-0.20706358), array(0.05499442), array(0.20160275), array(0.0011385), array(-0.16218059), array(-0.10453103), array(-0.03929314), array(-0.19958158), array(-0.03306489), array(0.26085468), array(0.24295046), array(0.03183827), array(0.20127936), array(-0.15557991), array(0.19118081), array(0.1021196), array(0.00111797), array(-0.00981701), array(-0.15042031), array(0.08196842), array(-0.13011536), array(0.17145289), array(0.10192479), array(0.01567546), array(0.14593138), array(0.05090766), array(0.08094841), array(-0.07744873), array(-0.07031627), array(0.04629074), array(-0.16097637), array(0.13795225), array(0.11284573), array(0.0589699), array(-0.04918522), array(0.04101822), array(-0.04564054), array(-0.11261781), array(0.12002761), array(0.16577642), array(0.11443499), array(0.01708446), array(0.0687497), array(0.11939119), array(-0.21055353), array(0.05157515), array(0.14405996), array(-0.16751556), array(-0.11712361), array(-0.03658351), array(-0.19152567), array(-0.20359878), array(0.20696077), array(-0.13358546), array(-0.01627022), array(-0.08281993), array(-0.04454049), array(-0.01605468), array(-0.06955726), array(0.0676535), array(-0.16620214), array(0.05095096), array(0.09061472), array(-0.16524318), array(0.11765692), array(0.06948661), array(0.07555307), array(0.00246678), array(-0.09140098), array(0.0071634), array(-0.04459941), array(0.13919143), array(-0.02212618), array(-0.06075185), array(-0.06037191), array(0.26175211), array(-0.08992471), array(-0.04239453), array(-0.1306687)]
    x2 = [float(item) for item in x2]

    constraints = [] 
    for m in range(1, M):
        lhs = np.array(x2) @ np.array(A[m]) @ np.array(x2) + 2* (np.array(C[m]) @ np.array(x2)) - b[m]
        constraints.append(lhs)
    print(constraints)

    temp = gp.QuadExpr()
    for i in range(N):
        for j in range(N):
            temp += A[0][i][j]*vars[i]*vars[j]
    for i in range(N):
        temp += 2*C[0][i]*vars[i]

    # model.addConstr(temp, GRB.LESS_EQUAL, obj)

    model.setObjective(temp, GRB.MINIMIZE)
    

    # vars = lda*np.array(x1) + (1-lda)*np.array(x2)
    for n in range(N): # lda*x1[n] + (1-lda)*
        model.addConstr(x2[n] == vars[n])
    # Constraints
    for m in range(1, M):
        constrain = gp.QuadExpr()
        for i in range(N):
            for j in range(N):
                constrain += A[m][i][j]*vars[i]*vars[j]
        for i in range(N):
            constrain += 2*C[m][i] * vars[i]
        model.addConstr(constrain <= b[m])


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
            solution.append(vars[i].X)
        print("Objective: {0}".format(model.objVal))
        print("Solution: {0}".format(solution))
        return (model.objVal, solution)
    elif model.status == GRB.TIME_LIMIT:
        print("Objective: {0}".format("TIME_LIMIT"))
        print("Solution: {0}".format("TIME_LIMIT"))
        return
    else:
        return
    



if __name__ == "__main__":
    with open("./QCQPData/" + str(sys.argv[1]) + ".json",'r') as f:
        load_dict = json.load(f)
        qcqp_optimize(load_dict["N"], load_dict["M"], load_dict["A"], load_dict["C"], load_dict["b"], load_dict["u"])

    