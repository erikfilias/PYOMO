from pyomo.environ import SolverFactory
import pyomo.environ as pe

model = pe.AbstractModel()

model.I = pe.Set()
model.J = pe.Set()
model.K = pe.Set()

model.s = pe.Param(model.K)
model.p = pe.Param(model.I, model.J)
model.cx = pe.Param(model.J, model.K)
model.cy = pe.Param(model.J, model.K)
# model.c = pe.Param(model.I, model.J, model.K)
model.d = pe.Param(model.I)
model.m = pe.Param(model.K)

# the next line declares a variable indexed by the set J
model.x = pe.Var(model.I, model.J, domain=pe.NonNegativeIntegers)
model.f = pe.Var(model.K, domain=pe.NonNegativeReals)


def obj_expression(model):
    return sum(model.p[i, j]*model.x[i, j] for i in model.I for j in model.J)-sum(model.m[k]*model.f[k] for k in model.K)


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.maximize)


def ax_constraint_ruleA(model, k):
    # return the expression for the constraint for k
    return model.s[k]-sum(model.cx[j, k]*model.x['X', j]for j in model.J)-sum(model.cy[j, k]*model.x['Y', j]for j in model.J) == model.f[k]
    # return model.s[k]-sum(model.c[i, j, k]*model.x[i, j] for i in model.I for j in model.J) == model.f[k]


# the next line creates one constraint for each member of the set model.K
model.AxbConstraintA = pe.Constraint(model.K, rule=ax_constraint_ruleA)

def ax_constraint_ruleB(model,i):
    # return the expression for the constraint for i
    return sum(model.x[i, j] for j in model.J) >= model.d[i]


# the next line creates one constraint for each member of the set model.I
model.AxbConstraintB = pe.Constraint(model.I, rule=ax_constraint_ruleB)

model.pprint()

instance = model.create_instance('D1_E3.dat')  # To choose instance
instance.pprint()

opt = SolverFactory('gurobi')

results = opt.solve(instance)

results.write()

# instance.obj_expression.pprint()
instance.x.pprint()
for v in instance.component_objects(pe.Var, active=True):
    print("Variable", v)
    varobject = getattr(instance, str(v))
    for index in varobject:
        print("   ", index, varobject[index].value)
# instance.OBJ.value
print("Profit = ", sum(instance.p[i, j]*instance.x[i, j].value for i in instance.I for j in instance.J)-sum(instance.m[k]*instance.f[k].value for k in instance.K))
# print("Profit = ", instance.bus_e[i].value)


# for j in instance.J:
#     print("x[", j, "] = ", '%10.4f' % (instance.x[j].value))

