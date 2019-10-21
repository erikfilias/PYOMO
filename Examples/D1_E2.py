from pyomo.environ import SolverFactory
import pyomo.environ as pe

model = pe.AbstractModel()

model.I = pe.Set()
model.J = pe.Set()

model.c = pe.Param(model.I, model.J)
model.f = pe.Param(model.I)
model.cf = pe.Param()
model.d = pe.Param()

# the next line declares a variable indexed by the set J
model.x = pe.Var(model.I, model.J, domain=pe.Binary)


def obj_expression(model):
    return sum((model.c[i, j]+model.cf)*model.x[i, j] for i in model.I for j in model.J)


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.minimize)


def ax_constraint_ruleA(model, i):
    # return the expression for the constraint for i
    return sum(model.x[i, j] for j in model.J) <= 1


# the next line creates one constraint for each member of the set model.I
model.AxbConstraintA = pe.Constraint(model.I, rule=ax_constraint_ruleA)

def ax_constraint_ruleB(model):
    # return the expression for the constraint for i
    return sum(j*model.f[i]*model.x[i, j] for i in model.I for j in model.J) >= model.d


# the next line creates one constraint for each member of the set model.I
model.AxbConstraintB = pe.Constraint(rule=ax_constraint_ruleB)

model.pprint()

instance = model.create_instance('D1_E2.dat')  # To choose instance
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
print("Profit = ", sum((instance.c[i, j]+instance.cf)*instance.x[i, j].value for i in instance.I for j in instance.J))
# print("Profit = ", instance.bus_e[i].value)


# for j in instance.J:
#     print("x[", j, "] = ", '%10.4f' % (instance.x[j].value))

