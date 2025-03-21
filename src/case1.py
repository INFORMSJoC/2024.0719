import sys
import gurobipy as gp
from gurobipy import GRB
import json
import numpy as np



model = gp.Model()
x1 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)
x2 = model.addVar(lb=0, ub=GRB.INFINITY, vtype=GRB.CONTINUOUS)

model.setObjective(x1+x2, GRB.MINIMIZE)

model.addConstr(7-x2, GRB.LESS_EQUAL, x1)
model.addConstr(3*(4-x2), GRB.LESS_EQUAL, x1)


try:
    model.optimize()
except gp.GurobiError:
    print("Optimize failed due to non-convexity")
    model.params.NonConvex = 2
    model.optimize()

print("Objective: {0}".format(model.objVal))
print("x1: {0}".format(x1.X))
print("x2: {0}".format(x2.X))