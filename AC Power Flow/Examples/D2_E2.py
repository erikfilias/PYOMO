from pyomo.environ import SolverFactory
import pyomo.environ as pe
import math

model = pe.AbstractModel()

model.N = pe.Set()    #Set of crew members


model.f     = pe.Param(model.N) #fishing skill
model.s     = pe.Param(model.N) #sailing skill
model.n     = pe.Param(model.N) #navigating skill
model.c     = pe.Param(model.N) #salary

model.x = pe.Var(model.N, domain=pe.Binary)


def obj_expression(model):
    return sum(model.c[i]*model.x[i] for i in model.N)


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.minimize)


def ax_constraint_ruleA(model):
    return sum(model.f[i]*model.x[i] for i in model.N) >= 15


model.AxbConstraintA = pe.Constraint(rule=ax_constraint_ruleA)
#

def ax_constraint_ruleB(model):
    return sum(model.s[i]*model.x[i] for i in model.N) >= 15


model.AxbConstraintB = pe.Constraint(rule=ax_constraint_ruleB)
#

def ax_constraint_ruleC(model):
    return sum(model.n[i]*model.x[i] for i in model.N) >= 15


model.AxbConstraintC = pe.Constraint(rule=ax_constraint_ruleC)
#

model.pprint()

instance = model.create_instance('D2_E2.dat')  # To choose instance
instance.pprint()

opt = SolverFactory('gurobi')

results = opt.solve(instance)

results.write()


# instance.dist.pprint()
instance.x.pprint()
# instance.cant.pprint()

print("Profit = ", sum(instance.c[i]*instance.x[i].value for i in instance.N))


