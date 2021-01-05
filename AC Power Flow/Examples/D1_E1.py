from pyomo.dataportal import DataPortal
from pyomo.environ import SolverFactory
import pyomo.environ as pe
import time

StartTime   = time.time()
model = pe.AbstractModel()
data = DataPortal()

model.I = pe.Set()
model.J = pe.Set()

model.c = pe.Param(model.I, model.J)
model.f = pe.Param(model.J)
model.p = pe.Param(model.I)
model.d = pe.Param(model.I)

# the next line declares a variable indexed by the set J
model.x = pe.Var(model.I, model.J, domain=pe.NonNegativeReals)


def obj_expression(model):
    return sum(model.p[i]*model.x[i, j] for i in model.I for j in model.J) - pe.summation(model.c, model.x)


model.OBJ = pe.Objective(rule=obj_expression, sense=pe.maximize)


def ax_constraint_ruleA(model, j):
    # return the expression for the constraint for j
    return sum(model.x[i, j] for i in model.I) <= model.f[j]


# the next line creates one constraint for each member of the set model.J
model.AxbConstraintA = pe.Constraint(model.J, rule=ax_constraint_ruleA)

def ax_constraint_ruleB(model, i):
    # return the expression for the constraint for i
    return sum(model.x[i, j] for j in model.J) >= model.d[i]


# the next line creates one constraint for each member of the set model.I
model.AxbConstraintB = pe.Constraint(model.I, rule=ax_constraint_ruleB)

model.pprint()

ModelingTime = time.time() - StartTime
StartTime   = time.time()


instance = model.create_instance('D1_E1.dat')  # To choose instance
instance.pprint()

ReadingTime = time.time() - StartTime
StartTime   = time.time()

opt = SolverFactory('gurobi')

results = opt.solve(instance)

SolvingTime = time.time() - StartTime
StartTime   = time.time()

results.write()

for v in instance.component_objects(pe.Var, active=True):
    print("Variable", v)
    varobject = getattr(instance, str(v))
    for index in varobject:
        print("   ", index, varobject[index].value)

print("Profit = ", sum(instance.p[i]*instance.x[i, j].value for i in instance.I for j in instance.J) - sum(instance.c[i, j]*instance.x[i, j].value for i in instance.I for j in instance.J))

WritingTime = time.time() - StartTime
StartTime   = time.time()

print('Modeling               time', ModelingTime       )
print('Reading DATA           time', ReadingTime  )
print('Solving                time', SolvingTime )
print('Writing                time', WritingTime )

