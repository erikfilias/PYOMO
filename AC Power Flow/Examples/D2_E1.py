from pyomo.environ import SolverFactory
import pyomo.environ as pe
import math

model = pe.AbstractModel()

model.L = pe.Set()    #Set of models
model.B = pe.Set(within=model.L*model.L)    #Set of branches


model.cap     = pe.Param(model.L) #Capacity of manufacture
model.D     = pe.Param(model.L) #Demand
model.dist  = pe.Param(model.B) #Distance
model.pos     = pe.Param(model.L) #Demand

model.cant = pe.Var(model.B)
model.inv = pe.Var(model.L, domain=pe.Binary)


def obj_expression(model):
    return sum(10*model.dist[i,j]*model.cant[i,j] for i,j in model.B)


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.minimize)


def ax_constraint_ruleA(model,k):
    return sum(model.cant[i,j] for i,j in model.B if j==k and i!=j) == model.D[k]


model.AxbConstraintA = pe.Constraint(model.L, rule=ax_constraint_ruleA)
#


def ax_constraint_ruleB(model):
    return sum(model.inv[k] for k in model.L) == 1


model.AxbConstraintB = pe.Constraint(rule=ax_constraint_ruleB)


def ax_constraint_ruleC(model,k):
    return sum(model.cant[i,j] for i,j in model.B if i==k) <= model.cap[k]  + model.inv[k]*1000


model.AxbConstraintC = pe.Constraint(model.L, rule=ax_constraint_ruleC)


def ax_constraint_ruleD(model,i,j):
    if i==j:
        return model.cant[i,j] == 0
    else:
        return model.cant[i,j] >= 0


model.AxbConstraintD = pe.Constraint(model.B, rule=ax_constraint_ruleD)


def ax_constraint_ruleE(model,k):
    if model.pos[k] == 0:
        return model.inv[k] == 0
    else:
        return model.inv[k] >= 0


model.AxbConstraintE = pe.Constraint(model.L, rule=ax_constraint_ruleE)


model.pprint()

instance = model.create_instance('D2_E1.dat')  # To choose instance
instance.pprint()

opt = SolverFactory('gurobi')

results = opt.solve(instance)

results.write()


# instance.dist.pprint()
instance.inv.pprint()
instance.cant.pprint()

print("Profit = ", sum(10*instance.dist[i,j]*instance.cant[i,j].value for i,j in instance.B))


