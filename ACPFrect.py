import pyomo.environ as pe
from pyomo.environ import SolverFactory

# import ACPFR
# from PYOMO import ACPFR
# include ACPFR
# include sistema14nos.dat
# load "sistema14nos"
# from pyomo.core import Param
# ... more constraints

model = pe.AbstractModel()

# SETS
model.BAR = pe.Set()
model.RAM = pe.Set(within=model.BAR * model.BAR, ordered=True)

# create PARAM BAR
model.tipo = pe.Param(model.BAR)
model.Pd = pe.Param(model.BAR)
model.Qd = pe.Param(model.BAR)
model.gshb = pe.Param(model.BAR)
model.bshb = pe.Param(model.BAR)
model.Pg0 = pe.Param(model.BAR)
model.Qg0 = pe.Param(model.BAR)
model.V0 = pe.Param(model.BAR)
model.th0 = pe.Param(model.BAR)

model.Pd_pu = pe.Param(model.BAR, initialize=model.Pd)
model.Qd_pu = pe.Param(model.BAR, initialize=model.Qd)
model.Pg0_pu = pe.Param(model.BAR, initialize=model.Pg0)
model.Qg0_pu = pe.Param(model.BAR, initialize=model.Qg0)

# create PARAM RAM
model.a = pe.Param(model.RAM)
model.r = pe.Param(model.RAM)
model.x = pe.Param(model.RAM)
model.g = pe.Param(model.RAM, initialize=0)
model.b = pe.Param(model.RAM, initialize=0)
model.fi = pe.Param(model.RAM)
model.smax = pe.Param(model.RAM)
model.bsh = pe.Param(model.RAM)

model.bsh_half = pe.Param(model.RAM, initialize=model.bsh)
model.fi_rad = pe.Param(model.RAM, initialize=model.fi)
model.smax_pu = pe.Param(model.RAM, initialize=model.smax)
model.a_init = pe.Param(model.RAM, initialize=model.a)

# create PARAM SYSTEM
model.Sbase = pe.Param()
model.Vnom = pe.Param()
model.Vmin = pe.Param()
model.Vmax = pe.Param()

# VARIABLES
# model.bus_voltage = Var(range(nbus),bounds = lambda model,i : (bus_voltage_min[bus[i].bustype], bus_voltage_max[bus[i].bustype]), initialize=1)
model.bus_e = pe.Var(model.BAR)  # Real part of voltage
model.bus_f = pe.Var(model.BAR, initialize=0)  # Imaginary part of voltage
model.branch_Ppara = pe.Var(model.RAM, initialize=0)
model.branch_Qpara = pe.Var(model.RAM, initialize=0)
model.branch_Pde = pe.Var(model.RAM, initialize=0)
model.branch_Qde = pe.Var(model.RAM, initialize=0)
model.bus_Pg = pe.Var(model.BAR)
model.bus_Qg = pe.Var(model.BAR)


# CALCULUS

def g_init(model, i, j):
    a = model.r[i, j] / (model.r[i, j] ** 2 + model.x[i, j] ** 2)
    return a


model.g = pe.Param(model.RAM, initialize=g_init)


def b_init(model, i, j):
    a = -model.x[i, j] / (model.r[i, j] ** 2 + model.x[i, j] ** 2)
    return a


model.b = pe.Param(model.RAM, initialize=b_init)


#
def bsh_init(model, i, j):
    a = (model.bsh[i, j]) / 2
    return a


model.bsh_half = pe.Param(model.RAM, initialize=bsh_init)


def fi_init(model, i, j):
    return (model.fi[i, j] * 3.14159265359) / 180


model.fi_rad = pe.Param(model.RAM, initialize=fi_init)


def smax_init(model, i, j):
    return (model.smax[i, j]) / model.Sbase


model.smax_pu = pe.Param(model.RAM, initialize=smax_init)


def a_init(model, i, j):
    if model.a[i, j] == 0:
        a = 1
    else:
        a = 1 / (model.a[i, j])
    return a


model.a_init = pe.Param(model.RAM, initialize=a_init)


def Pd_init(model, i):
    a = (model.Pd[i]) / model.Sbase
    return a


model.Pd_pu = pe.Param(model.BAR, initialize=Pd_init)


def Qd_init(model, i):
    a = (model.Qd[i]) / model.Sbase
    return a


model.Qd_pu = pe.Param(model.BAR, initialize=Qd_init)


def Pg0_init(model, i):
    a = (model.Pg0[i]) / model.Sbase
    return a


model.Pg0_pu = pe.Param(model.BAR, initialize=Pg0_init)


def Qg0_init(model, i):
    a = (model.Qg0[i]) / model.Sbase
    return a


model.Qg0_pu = pe.Param(model.BAR, initialize=Qg0_init)


def e_init(model, i):
    a2 = model.Vnom
    return a2


model.bus_e = pe.Var(model.BAR, initialize=e_init)


def f_init(model, i):
    a = 0
    return a


model.bus_f = pe.Var(model.BAR, initialize=f_init)


def Pg_init(model, i):
    if model.tipo[i] != 3:
        a2 = model.Pg0_pu[i]
        model.bus_Pg[i].fixed = True
    else:
        a2 = 0
    return a2


model.bus_Pg = pe.Var(model.BAR, initialize=Pg_init)


def Qg_init(model, i):
    if model.tipo[i] == 0:
        a2 = model.Qg0_pu[i]
        model.bus_Qg[i].fixed = True
    else:
        a2 = 0
    return a2


model.bus_Qg = pe.Var(model.BAR, initialize=Qg_init)


###########################################################################################
# OBJECTIVE
def TotalLoss(model):
    return sum(model.g[i, j] * (model.a_init[i, j] ** 2 * (model.bus_e[i] ** 2 + model.bus_f[i] ** 2) + \
                                (model.bus_e[j] ** 2 + model.bus_f[j] ** 2) - 2 * model.a_init[i, j] * \
                                (model.bus_e[i] * model.bus_e[j] + model.bus_f[i] * model.bus_f[j])) \
               for i, j in model.RAM)


model.Loss = pe.Objective(rule=TotalLoss)


##CONSTRAINTS
# Active Power Balance
def A_p_balance_rule(model, k):
    return model.bus_Pg[k] - model.Pd_pu[k] - \
           sum(model.branch_Ppara[j, i]
               for j, i in model.RAM if k == i) - \
           sum(model.branch_Pde[i, j]
               for i, j in model.RAM if k == i) + \
           (model.bus_e[k] ** 2 + model.bus_f[k] ** 2) * model.gshb[k] == 0


model.A_p_balance_rule = pe.Constraint(model.BAR, rule=A_p_balance_rule)


# Reactive Power Balance
def B_q_balance_rule(model, k):
    return model.bus_Qg[k] - model.Qd_pu[k] - \
           sum(model.branch_Qpara[j, i]
               for j, i in model.RAM if k == i) - \
           sum(model.branch_Qde[i, j]
               for i, j in model.RAM if k == i) + \
           (model.bus_e[k] ** 2 + model.bus_f[k] ** 2) * model.bshb[k] == 0


model.B_q_balance_constr = pe.Constraint(model.BAR, rule=B_q_balance_rule)


# Pde Constraint
def C_Pde_rule(model, i, j):
    return (model.branch_Pde[i, j] == model.a_init[i, j] ** 2 * model.g[i, j] * (
            model.bus_e[i] ** 2 + model.bus_f[i] ** 2) - \
            model.a_init[i, j] * model.g[i, j] * (model.bus_e[i] * model.bus_e[j] + model.bus_f[i] * model.bus_f[j]) + \
            model.a_init[i, j] * model.b[i, j] * (model.bus_e[i] * model.bus_f[j] - model.bus_e[j] * model.bus_f[i]))


model.C_Pde_rule = pe.Constraint(model.RAM, rule=C_Pde_rule)


# Ppara Constraint
def D_Ppara_rule(model, i, j):
    return (model.branch_Ppara[i, j] == model.g[i, j] * (model.bus_e[j] ** 2 + model.bus_f[j] ** 2) - \
            model.a_init[i, j] * model.g[i, j] * (model.bus_e[i] * model.bus_e[j] + model.bus_f[i] * model.bus_f[j]) - \
            model.a_init[i, j] * model.b[i, j] * (model.bus_e[i] * model.bus_f[j] - model.bus_e[j] * model.bus_f[i]))


model.D_Ppara_rule = pe.Constraint(model.RAM, rule=D_Ppara_rule)


# Qpara Constraint
def E_Qde_rule(model, i, j):
    return (model.branch_Qde[i, j] == -model.a_init[i, j] ** 2 * (model.b[i, j] + model.bsh_half[i, j]) * (
            model.bus_e[i] ** 2 + model.bus_f[i] ** 2) + \
            model.a_init[i, j] * model.g[i, j] * (model.bus_e[i] * model.bus_f[j] - model.bus_e[j] * model.bus_f[i]) + \
            model.a_init[i, j] * model.b[i, j] * (model.bus_e[i] * model.bus_e[j] + model.bus_f[i] * model.bus_f[j]))


model.E_Qde_rule = pe.Constraint(model.RAM, rule=E_Qde_rule)


# Qde Constraint
def F_Qpara_rule(model, i, j):
    return (model.branch_Qpara[i, j] == -(model.b[i, j] + model.bsh_half[i, j]) * (
            model.bus_e[j] ** 2 + model.bus_f[j] ** 2) - \
            model.a_init[i, j] * model.g[i, j] * (model.bus_e[i] * model.bus_f[j] - model.bus_e[j] * model.bus_f[i]) + \
            model.a_init[i, j] * model.b[i, j] * (model.bus_e[i] * model.bus_e[j] + model.bus_f[i] * model.bus_f[j]))


model.F_Qpara_rule = pe.Constraint(model.RAM, rule=F_Qpara_rule)


# Generation Voltage
def G_Vg_rule(model, i):
    switch = model.tipo[i]
    if switch != 0:
        return (model.bus_e[i] ** 2 + model.bus_f[i] ** 2 == model.V0[i] ** 2)
    else:
        return pe.Constraint.Skip


model.G_Vg_rule = pe.Constraint(model.BAR, rule=G_Vg_rule)


# Angle
def H_angle_rule(model, i):
    switch = model.tipo[i]
    if switch == 3:
        return (model.bus_f[i] == model.bus_e[i] * pe.tan(model.th0[i] * 3.14159265359 / 180))
    else:
        return pe.Constraint.Skip


model.H_angle_rule = pe.Constraint(model.BAR, rule=H_angle_rule)

##Limit Current Constraint
# def F_limit_current_rule(model,i,j):
#    return (0,model.branch_Iij_sqr[i,j],None)
# model.F_limit_current_rule = pe.Constraint(model.RAM, rule = F_limit_current_rule)
#
##Limit Voltage Constraint
# def G_limit_voltage_rule(model,i):
#    return (0,model.bus_voltage_sqr[i],None)
# model.G_limit_voltage_rule = pe.Constraint(model.BAR, rule = G_limit_voltage_rule)

# model.pprint()


# model.pprint()

instance = model.create_instance('sistema14nos.dat')
instance.pprint()

opt = SolverFactory('ipopt')

results = opt.solve(instance)

results.write()

print("------------------------------------------------------------------------------------------------------------")
print("-------------------------------SUMMARY----------------------------------------------------------------------")
print("------------------------------------------------------------------------------------------------------------")
print("BUS RESULTS")
print("------------------------------------------------------------------------------------------------------------")
print("   Bus      V[pu]  Theta[Degree]     Pg[MW]     Qg[MVAr]       Pd[MW]     Qd[MVAr]      Gsh[MW]    Bsh[MVAr]")
# for i in range(len(model.bus_angle)): print("Bus_Angle "+str(i)+":" + str(model.bus_angle[i].value*180/3.14159265359))
for i in instance.BAR:
    a = (instance.bus_e[i].value ** 2 + instance.bus_f[i].value ** 2)
    b = 180 / 3.14159265359 * pe.atan((instance.bus_f[i].value) / instance.bus_e[i].value)
    print('%5d  %10.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f'
          % (i, a ** 0.5, b,
             instance.Sbase * instance.bus_Pg[i].value, instance.Sbase * instance.bus_Qg[i].value,
             instance.Sbase * instance.Pd_pu[i], instance.Sbase * instance.Qd_pu[i],
             instance.Sbase * instance.gshb[i] * a,
             instance.Sbase * instance.bshb[i] * a))

a = 0
b = 0
c = 0
d = 0
e = 0
f = 0
for i in instance.BAR:
    a1 = (instance.bus_e[i].value ** 2 + instance.bus_f[i].value ** 2)
    a = a + instance.bus_Pg[i].value
    b = b + instance.bus_Qg[i].value
    c = c + instance.Pd[i]
    d = d + instance.Qd[i]
    e = e + instance.gshb[i] * a1
    f = f + instance.bshb[i] * a1
print(
        "------------------------------------------------------------------------------------------------------------")
print('TOTAL %37.4f %12.4f %12.4f %12.4f %12.4f %12.4f'
      % (instance.Sbase * a, instance.Sbase * b, c, d, instance.Sbase * e, instance.Sbase * f))
#
print("------------------------------------------------------------------------------------------------------------")
print("BRANCH RESULTS")
print("------------------------------------------------------------------------------------------------------------")
print("    i     j        Pij[MW]         Pji[MW]       Qij[MVAr]       Qji[MVAr]         Pls[MW]       Qls[MVAr]")
for i, j in instance.RAM:
    print('%5d %5d %15.4f %15.4f %15.4f %15.4f %15.4f %15.4f'
          % (i, j,
             instance.Sbase * (instance.branch_Pde[i, j].value),
             instance.Sbase * instance.branch_Ppara[i, j].value,
             instance.Sbase * (instance.branch_Qde[i, j].value),
             instance.Sbase * (instance.branch_Qpara[i, j].value),
             ((instance.Sbase * (instance.branch_Pde[i, j].value + instance.branch_Ppara[i, j].value)) ** 2) ** 0.5,
             ((instance.Sbase * (instance.branch_Qde[i, j].value + instance.branch_Qpara[i, j].value)) ** 2) ** 0.5))

a = 0
b = 0
for i, j in instance.RAM:
    a = a + instance.branch_Ppara[i, j].value + instance.branch_Pde[i, j].value
    b = b + instance.branch_Qpara[i, j].value + instance.branch_Qde[i, j].value
print("------------------------------------------------------------------------------------------------------------")
print('TOTAL %85.4f %15.4f'
      % (instance.Sbase * a, instance.Sbase * b))
print("------------------------------------------------------------------------------------------------------------")