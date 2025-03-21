import sys
import math
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import json

from scipy import linalg
model = gp.Model()

A = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(2)] for i in range(3)]
B = [[model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS) for j in range(3)] for i in range(2)]

AB = [[8, 2, -2],[2,5,4],[-2,4,5]]


e_vals,e_vecs = linalg.eig(AB)
print(e_vecs)
print(e_vals)

for i in range(3):
    for j in range(3):
        model.addConstr(A[i][0]*B[0][j] + A[i][1]*B[1][j] == AB[i][j])

# BA = [[9,0],[0,9]]
# for i in range(2):
#     for j in range(2):
#         model.addConstr(B[i][0]*A[0][j] + B[i][1]*A[1][j] + B[i][2]*A[2][j] == BA[i][j])

model.addConstr(A[0][0] >= 0)

model.setObjective(0, GRB.MINIMIZE)
try:
    model.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")
    # Solve bilinear model
    model.params.NonConvex = 2
    model.optimize()


A_X = [[0 for i in range(2)] for i in range(3)]
B_X = [[0 for i in range(3)] for i in range(2)]
for i in range(3):
    for j in range(2):
        A_X[i][j] = A[i][j].X

for i in range(2):
    for j in range(3):
        B_X[i][j] = B[i][j].X

print(A_X)
print(B_X)
# for i in range(3):
#     for j in range(3):
#         print(A[i][0].X*B[0][j].X + A[i][1].X*B[1][j].X)

for i in range(2):
    for j in range(2):
        print(B[i][0].X*A[0][j].X + B[i][1].X*A[1][j].X + B[i][2].X*A[2][j].X)