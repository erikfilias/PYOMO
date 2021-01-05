from pyomo.environ import SolverFactory
import pyomo.environ as pe
import math

model = pe.AbstractModel()

model.N = pe.Set()    #Set of dials
model.J = pe.Set()    #Set of options

model.n = pe.Param(model.N, model.J) #values of options on each dial

model.x = pe.Var(model.N, model.J, domain=pe.Binary)
model.alpha1 = pe.Var(domain=pe.NonNegativeReals)


def obj_expression(model):
    return model.alpha1


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.minimize)
model.c = pe.ConstraintList()

def ax_constraint_ruleA(model):
    return sum(model.n[i,j]*model.x[i,j] for j in model.J for i in model.N) == 419


model.AxbConstraintA = pe.Constraint(rule=ax_constraint_ruleA)
#

def ax_constraint_ruleB(model,i):
    return sum(model.x[i,j] for j in model.J) == 1


model.AxbConstraintB = pe.Constraint(model.N, rule=ax_constraint_ruleB)
#

model.pprint()

instance = model.create_instance('D2_E3.dat')  # To choose instance
instance.pprint()

opt = SolverFactory('gurobi')

results = opt.solve(instance)

results.write()


# instance.dist.pprint()
instance.x.pprint()
# instance.cant.pprint()

# print("Profit = ", sum(instance.c[i]*instance.x[i].value for i in instance.N))

# Iterate to eliminate the previously found solution
# for i in range(5):
#     expr = 0
#     for j in instance.x:
#         if pe.value(instance.x[j]) == 0:
#             expr += instance.x[j]
#         else:
#             expr += (1-instance.x[j])
#     instance.c.add( expr >= 1 )
#     results = opt.solve(instance)
#     print ("\n===== iteration",i)
#     instance.display()

