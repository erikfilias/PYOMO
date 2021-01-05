# -*- coding: utf-8 -*-
"""
Created on Wed Dec 27 22:23:41 2017

@author: Erik
"""
import pyomo.environ as pe

model = pe.AbstractModel()

#SETS
model.BAR = pe.Set()
model.RAM = pe.Set(within=model.BAR*model.BAR, ordered=True)
model.L = pe.RangeSet(20)

#create PARAM BAR
model.tipo = pe.Param(model.BAR)
model.Pd = pe.Param(model.BAR)
model.Qd = pe.Param(model.BAR)
model.gshb = pe.Param(model.BAR)
model.bshb = pe.Param(model.BAR)
model.Pg0 = pe.Param(model.BAR)
model.Qg0 = pe.Param(model.BAR)
model.V0 = pe.Param(model.BAR)
model.th0 = pe.Param(model.BAR)

model.Pd_pu = pe.Param(model.BAR, initialize = model.Pd)
model.Qd_pu = pe.Param(model.BAR, initialize = model.Qd)
model.Pg0_pu = pe.Param(model.BAR, initialize = model.Pg0)
model.Qg0_pu = pe.Param(model.BAR, initialize = model.Qg0)

#create PARAM RAM
model.a = pe.Param(model.RAM)
model.r = pe.Param(model.RAM)
model.x = pe.Param(model.RAM)
model.z2 = pe.Param(model.RAM, initialize = 0)
model.fi = pe.Param(model.RAM)
model.smax = pe.Param(model.RAM)
model.bsh = pe.Param(model.RAM)

model.bsh_half = pe.Param(model.RAM, initialize = model.bsh)
model.fi_rad = pe.Param(model.RAM, initialize = model.fi)
model.smax_pu = pe.Param(model.RAM, initialize = model.smax)
model.a_init = pe.Param(model.RAM, initialize = model.a)

#create PARAM SYSTEM
model.Sbase = pe.Param()
model.Vnom = pe.Param()
model.Vmin = pe.Param()
model.Vmax = pe.Param()

#create PARAM LINEARIZATION
model.delta_S = pe.Param(model.RAM)
model.m = pe.Param(within=model.RAM*model.L)

#create VAR LINEARIZATION
model.delta_P = pe.Var(model.RAM,model.L)
model.delta_Q = pe.Var(model.RAM,model.L)
model.P_max = pe.Var(model.RAM)
model.P_min = pe.Var(model.RAM)
model.Q_max = pe.Var(model.RAM)
model.Q_min = pe.Var(model.RAM)

#VARIABLES
#model.bus_voltage = Var(range(nbus),bounds = lambda model,i : (bus_voltage_min[bus[i].bustype], bus_voltage_max[bus[i].bustype]), initialize=1)
model.bus_voltage_sqr = pe.Var(model.BAR)
model.bus_angle = pe.Var(model.BAR, initialize=0)
model.branch_Iij_sqr = pe.Var(model.RAM)
model.branch_P = pe.Var(model.RAM)
model.branch_Q = pe.Var(model.RAM)
model.bus_Pg = pe.Var(model.BAR)
model.bus_Qg = pe.Var(model.BAR)

#CALCULUS

def z2_init(model,i,j):
    a = model.r[i,j]**2 + model.x[i,j]**2
    return a
model.z2 = pe.Param(model.RAM, initialize = z2_init)
#
def bsh_init(model,i,j):
    a = (model.bsh[i,j])/2
    return a
model.bsh_half = pe.Param(model.RAM, initialize = bsh_init)

def fi_init(model,i,j):
    return (model.fi[i,j]*3.14159265359)/180
model.fi_rad = pe.Param(model.RAM, initialize = fi_init)

def smax_init(model,i,j):
    return (model.smax[i,j])/model.Sbase
model.smax_pu = pe.Param(model.RAM, initialize = smax_init)

def a_init(model,i,j):
    if model.a[i,j] == 0:
        a = 1
    else:
        a = 1/(model.a[i,j])
    return a
model.a_init = pe.Param(model.RAM, initialize = a_init)

def Pd_init(model,i):
    a = (model.Pd[i])/model.Sbase
    return a
model.Pd_pu = pe.Param(model.BAR, initialize = Pd_init)

def Qd_init(model,i):
    a = (model.Qd[i])/model.Sbase
    return a
model.Qd_pu = pe.Param(model.BAR, initialize = Qd_init)

def Pg0_init(model,i):
    a = (model.Pg0[i])/model.Sbase
    return a
model.Pg0_pu = pe.Param(model.BAR, initialize = Pg0_init)

def Qg0_init(model,i):
    a = (model.Qg0[i])/model.Sbase
    return a
model.Qg0_pu = pe.Param(model.BAR, initialize = Qg0_init)

def Theta_init(model,i):
    switch = model.tipo[i]
    if switch == 0:
        a1 = 0
    elif switch == 2:
        a1 = 0
    else:
        a1 = model.th0[i]*3.14159265359/180;
        model.bus_angle[i].fixed = True
    return a1 
model.bus_angle = pe.Var(model.BAR, initialize = Theta_init)

def Voltage_init(model,i):
    if model.tipo[i] == 0:
        a2 = model.Vnom**2
        model.bus_voltage_sqr[i].fixed = False
    else:
        a2 = model.V0[i]**2
        model.bus_voltage_sqr[i].fixed = True
    return a2
model.bus_voltage_sqr = pe.Var(model.BAR, initialize = Voltage_init)

def Isqr_init(model,i,j):
    a2 = 0
    return a2
model.branch_Iij_sqr = pe.Var(model.RAM, initialize = Isqr_init)

def Pij_init(model,i,j):
    a2 = 0
    return a2
model.branch_P = pe.Var(model.RAM, initialize = Pij_init)

def Qij_init(model,i,j):
    a2 = 0
    return a2
model.branch_Q = pe.Var(model.RAM, initialize = Qij_init)

def Pg_init(model,i):
    if model.tipo[i] == 0 or model.tipo[i] == 2:
        a2 = model.Pg0_pu[i]
        model.bus_Pg[i].fixed = True
    else:
        a2 = 0
    return a2
model.bus_Pg = pe.Var(model.BAR, initialize = Pg_init)

def Qg_init(model,i):
    if model.tipo[i] == 0:
        a2 = model.Qg0_pu[i]
        model.bus_Qg[i].fixed = True
    else:
        a2 = 0
    return a2
model.bus_Qg = pe.Var(model.BAR, initialize = Qg_init)

def delta_S_init(model,i,j):
    a = (model.smax_pu[i,j])/len(model.L)
    return a
model.delta_S = pe.Param(model.RAM, initialize = delta_S_init)

def m_init(model,i,j,l):
    a = (2*l-1)*(model.delta_S[i,j])
    return a
model.m = pe.Param(model.RAM,model.L, initialize = m_init)

#####################################################################################################################################
#OBJECTIVE
#####################################################################################################################################
def TotalLoss(model):
    return sum(model.r[i,j]*(model.branch_Iij_sqr[i,j]) \
               for i,j in model.RAM)
model.Loss = pe.Objective(rule=TotalLoss)

##CONSTRAINTS
#Active Power Balance
def A_p_balance_rule(model,k):
    return model.bus_Pg[k] - model.Pd_pu[k] + \
            sum(model.branch_P[j,i] 
                for j,i in model.RAM if k == i) - \
            sum(model.branch_P[i,j] + model.r[i,j]*model.branch_Iij_sqr[i,j] 
                for i,j in model.RAM if k == i) + \
            model.bus_voltage_sqr[k]*model.gshb[k] == 0 
model.A_p_balance_rule = pe.Constraint(model.BAR, rule = A_p_balance_rule)

#Reactive Power Balance
def B_q_balance_rule(model, k):
    return model.bus_Qg[k] - model.Qd_pu[k] + \
            sum(model.branch_Q[j,i] + model.bsh_half[j,i]*model.bus_voltage_sqr[k] 
            for j,i in model.RAM if k == i) - \
            sum(model.branch_Q[i,j] - model.bsh_half[i,j]*model.bus_voltage_sqr[k] + \
                model.x[i,j]*(model.branch_Iij_sqr[i,j]) 
                for i,j in model.RAM if k == i) + \
            model.bus_voltage_sqr[k]*model.bshb[k] == 0 

model.B_q_balance_constr = pe.Constraint(model.BAR, rule = B_q_balance_rule)

#Voltage Constraint
def C_voltage_rule(model,i,j):   
    return (model.bus_voltage_sqr[i]*model.a_init[i,j]**2 - \
            model.bus_voltage_sqr[j] == 2*(model.r[i,j]*model.branch_P[i,j] + \
            model.x[i,j]*model.branch_Q[i,j]) + model.z2[i,j]*model.branch_Iij_sqr[i,j])
model.C_voltage_rule = pe.Constraint(model.RAM, rule = C_voltage_rule)

#Angle Constraint
def D_angle_rule(model,i,j):    
    return ((model.V0[i]*model.V0[j])*model.a_init[i,j] * \
            (model.bus_angle[i]-model.bus_angle[j] + \
                model.fi_rad[i,j]) == model.x[i,j]*model.branch_P[i,j] - \
            model.r[i,j]*model.branch_Q[i,j])
model.D_angle_rule = pe.Constraint(model.RAM, rule = D_angle_rule)

#Current Constraint
def E_current_rule(model,i,j):
    return ((model.V0[j]**2)*model.branch_Iij_sqr[i,j] == sum(model.m[i,j,l]*model.delta_P[i,j,l] for l in model.L) + \
            sum(model.m[i,j,l]*model.delta_Q[i,j,l] for l in model.L))
model.E_current_rule = pe.Constraint(model.RAM, rule = E_current_rule)

#Pij_Limit
def F_Pij_limit_rule(model,i,j):
    return model.P_max[i,j]-model.P_min[i,j] == model.branch_P[i,j]
model.F_Pij_limit_rule = pe.Constraint(model.RAM, rule = F_Pij_limit_rule)

#Delta_Pij_Limit
def G_delta_Pij_limit_rule(model,i,j):
    return model.P_max[i,j]+model.P_min[i,j] == sum(model.delta_P[i,j,l] for l in model.L)
model.G_delta_Pij_limit_rule = pe.Constraint(model.RAM, rule = G_delta_Pij_limit_rule)

#Pij_Limit
def H_Qij_limit_rule(model,i,j):
    return model.Q_max[i,j]-model.Q_min[i,j] == model.branch_Q[i,j]
model.H_Qij_limit_rule = pe.Constraint(model.RAM, rule = H_Qij_limit_rule)

#Delta_Pij_Limit
def I_delta_Qij_limit_rule(model,i,j):
    return model.Q_max[i,j]+model.Q_min[i,j] == sum(model.delta_Q[i,j,l] for l in model.L)
model.I_delta_Qij_limit_rule = pe.Constraint(model.RAM, rule = I_delta_Qij_limit_rule)

#Limit Delta_Sij_limit_P1
def J_delta_Sij_limit_P1_rule(model,i,j,l):
    return (0,model.delta_P[i,j,l],model.delta_S[i,j])
model.J_delta_Sij_limit_P1_rule = pe.Constraint(model.RAM, model.L, rule = J_delta_Sij_limit_P1_rule)

#Limit Delta_Sij_limit_Q1
def K_delta_Sij_limit_Q1_rule(model,i,j,l):
    return (0,model.delta_Q[i,j,l],model.delta_S[i,j])
model.K_delta_Sij_limit_Q1_rule = pe.Constraint(model.RAM, model.L, rule = K_delta_Sij_limit_Q1_rule)

#Limit Pmax_limit
def L_Pmax_limit_rule(model,i,j):
    return (0,model.P_max[i,j],None)
model.L_Pmax_limit_rule = pe.Constraint(model.RAM, rule = L_Pmax_limit_rule)

#Limit Pmin_limit
def M_Pmin_limit_rule(model,i,j):
    return (0,model.P_min[i,j],None)
model.M_Pmin_limit_rule = pe.Constraint(model.RAM, rule = M_Pmin_limit_rule)

#Limit Qmax_limit
def N_Qmax_limit_rule(model,i,j):
    return (0,model.Q_max[i,j],None)
model.N_Qmax_limit_rule = pe.Constraint(model.RAM, rule = N_Qmax_limit_rule)

#Limit Pmin_limit
def O_Qmin_limit_rule(model,i,j):
    return (0,model.Q_min[i,j],None)
model.O_Qmin_limit_rule = pe.Constraint(model.RAM, rule = O_Qmin_limit_rule)

#Limit Current Constraint
def P_limit_current_rule(model,i,j):
    return (0,model.branch_Iij_sqr[i,j],None)
model.P_limit_current_rule = pe.Constraint(model.RAM, rule = P_limit_current_rule)

#Limit Voltage Constraint
def Q_limit_voltage_rule(model,i):
    return (0,model.bus_voltage_sqr[i],None)
model.Q_limit_voltage_rule = pe.Constraint(model.BAR, rule = Q_limit_voltage_rule)

#model.pprint()