import gurobipy as gp
from gurobipy import GRB

try:
    # Create a new model
    m = gp.Model("mip1")

    # Create variables
    x = m.addVar(vtype=GRB.CONTINUOUS, name="x")
    y = m.addVar(vtype=GRB.CONTINUOUS, name="y")
    z = m.addVar(vtype=GRB.CONTINUOUS, name="z")

    # Set objective
    m.setObjective(x + y + 2 * z, GRB.MAXIMIZE)

    # Add constraint: x + 2 y + 3 z <= 4
    m.addConstr(x + 2 * y + 3 * z <= 4, "c0")

    # Add constraint: x + y >= 1
    m.addConstr(x + y >= 1, "c1")

    

    # Optimize model
    m.optimize()

    for v in m.getVars():
        print(f"{v.VarName} {v.X}")

    print(f"Obj: {m.ObjVal}")

    print(m.getCoeff(m.getConstrs()[0],x))
    print(m.getCol(y).size())
    print(m.getCol(z).size())

except gp.GurobiError as e:
    print(f"Error code {e.errno}: {e}")

except AttributeError:
    print("Encountered an attribute error")