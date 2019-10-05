# StarGen Lite Stochastic Daily Unit Commitment of Thermal and Hydro Units (SDUC)

# Developed by

#    Andres Ramos
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu

#    MIT Energy Initiative
#    Massachusetts Institute of Technology
#    arght@mit.edu

#    January 22, 2019

# from __future__ import division         #Division returns float
import numpy    as np
import pandas   as pd
import openpyxl as op
import numpy.matlib

from pyomo.environ import *
from pyomo.core import Var  # To fix binary variables
from pyomo.core import Constraint  # To get dual   variables
from pyomo.core.kernel import *

# system dimensions
N = 24  # hours
SC = 3  # scenarios
G = 13  # thermal and hydro generating units

# reading data from Excel
InputFile = 'StarGenLite_SDUC.xlsm'
dfDemand = pd.read_excel(InputFile, sheet_name='DemandReserveIG', header=0, usecols='D', nrows=N)
dfOperReserve = pd.read_excel(InputFile, sheet_name='DemandReserveIG', header=0, usecols='I', nrows=N)
dfOperReserveUp = pd.read_excel(InputFile, sheet_name='DemandReserveIG', header=0, usecols='N', nrows=N)
dfOperReserveDw = pd.read_excel(InputFile, sheet_name='DemandReserveIG', header=0, usecols='S', nrows=N)
dfIntermGen = pd.read_excel(InputFile, sheet_name='DemandReserveIG', header=4, usecols='X:Z', nrows=N)
dfThermalGen = pd.read_excel(InputFile, sheet_name='Generation', header=5, usecols='C:R', nrows=G, index_col=0)

# Parameter
pScenProb = [0.3, 0.5, 0.2]  # probabilities of scenarios
pENSCost = 10  # cost of energy not served      [MEur per GWh]
pCO2Cost = 5  # cost of CO2 emission           [Eur per tCO2]

# scaling of parameters to GW and MEur
pDemand = dfDemand.values * 1e-3  # hourly load                    [GW]
pOperReserve = dfOperReserve.values * 1e-3  # hourly operating reserve       [GW]
pOperReserveUp = dfOperReserveUp.values * 1e-3  # hourly operating reserve up    [GW]
pOperReserveDw = dfOperReserveDw.values * 1e-3  # hourly operating reserve down  [GW]
pIntermGen = dfIntermGen.values * 1e-3  # stochastic IG generation       [GW]
pMaxProd = dfThermalGen.iloc[:, 0].values * 1e-3  # maximum output                 [GW]
pMinProd = dfThermalGen.iloc[:, 1].values * 1e-3  # minimum output                 [GW]
pIniOut = dfThermalGen.iloc[:, 2].values * 1e-3  # initial output                 [GW]
pRampUp = dfThermalGen.iloc[:, 3].values * 1e-3  # ramp up                        [GW/h]
pRampDw = dfThermalGen.iloc[:, 4].values * 1e-3  # ramp down                      [GW/h]
pSlopeVarCost = dfThermalGen.iloc[:, 6].values * 1e-3 * dfThermalGen.iloc[:, 5].values + dfThermalGen.iloc[:,
                                                                                         8].values * 1e-3  # slope variable cost            [MEur/GWh]
pInterVarCost = dfThermalGen.iloc[:, 7].values * 1e-6 * dfThermalGen.iloc[:,
                                                        5].values  # intercept variable cost        [MEur/GWh]
pMinTU0 = dfThermalGen.iloc[:, 9].values  # minimum time up                [h]
pMinTD0 = dfThermalGen.iloc[:, 10].values  # minimum time down              [h]
pEmissionCost = dfThermalGen.iloc[:, 11].values * 1e-3 * pCO2Cost  # emission cost                  [MEur/GWh]
pStartupCost = dfThermalGen.iloc[:, 12].values * 1e-6 * dfThermalGen.iloc[:,
                                                        5].values  # startup cost                   [MEur]
pShutDownCost = dfThermalGen.iloc[:, 13].values * 1e-6 * dfThermalGen.iloc[:,
                                                         5].values  # shutdown cost                  [MEur]

pIniUC = np.zeros(G)
pMinTU = np.zeros(G)
pMinTD = np.zeros(G)

for g in range(G):
    if pIniOut[g] >= pMinProd[g]:
        pIniUC[g] = 1
    else:
        pIniUC[g] = 0
    if pMinTU0[g] == 0:
        pMinTU[g] = 1
    else:
        pMinTU[g] = round(pMinTU0[g])
    if pMinTD0[g] == 0:
        pMinTD[g] = 1
    else:
        pMinTD[g] = round(pMinTD0[g])

# Stochastic Daily Unit Commitment (UC) SDUC
mSDUC = ConcreteModel()

# sets
mSDUC.N = Set(initialize=RangeSet(0, N - 1), doc='hours')
mSDUC.G = Set(initialize=RangeSet(0, G - 1), doc='units')
mSDUC.SC = Set(initialize=RangeSet(0, SC - 1), doc='scenarios')

pProduct1_up = np.zeros((SC, N, G))
pProduct_up = np.zeros((SC, N, G))
pIG_up = np.zeros((SC, N))
pENS_up = np.zeros((SC, N))

# definition of variable bounds
for sc in range(SC):
    for n in range(N):
        pIG_up[sc, n] = pIntermGen[n, sc]
        pENS_up[sc, n] = pDemand[n]
        for g in range(G):
            pProduct1_up[sc, n, g] = pMaxProd[g] - pMinProd[g]
            pProduct_up[sc, n, g] = pMaxProd[g]


# variable bounds
def vProduct_bd(mSDUC, sc, n, g):
    return ((0.0, pProduct_up[sc, n, g]))


def vProduct1_bd(mSDUC, sc, n, g):
    return ((0.0, pProduct1_up[sc, n, g]))


def vIG_bd(mSDUC, sc, n):
    return ((0.0, pIG_up[sc, n]))


def vENS_bd(mSDUC, sc, n):
    return ((0.0, pENS_up[sc, n]))


# variables
mSDUC.vProduct = Var(mSDUC.SC, mSDUC.N, mSDUC.G, within=NonNegativeReals, bounds=vProduct_bd,
                     doc='output of the unit            [GW]')
mSDUC.vProduct1 = Var(mSDUC.SC, mSDUC.N, mSDUC.G, within=NonNegativeReals, bounds=vProduct1_bd,
                      doc='output of the unit > min load [GW]')
mSDUC.vIG = Var(mSDUC.SC, mSDUC.N, within=NonNegativeReals, bounds=vIG_bd, doc='intermittent generation       [GW]')
mSDUC.vENS = Var(mSDUC.SC, mSDUC.N, within=NonNegativeReals, bounds=vENS_bd, doc='energy not served             [GW]')
mSDUC.vCommitt = Var(mSDUC.N, mSDUC.G, within=Binary, doc='commitment of the unit       {0,1}')
mSDUC.vStartup = Var(mSDUC.N, mSDUC.G, within=Binary, doc='startup    of the unit       {0,1}')
mSDUC.vShutdown = Var(mSDUC.N, mSDUC.G, within=Binary, doc='shutdown   of the unit       {0,1}')

mSDUC.obj = Objective(expr=sum(pENSCost * mSDUC.vENS[sc, n] * pScenProb[sc] for sc in mSDUC.SC for n in mSDUC.N) +
                           sum(pSlopeVarCost[g] * mSDUC.vProduct[sc, n, g] * pScenProb[sc] for sc in mSDUC.SC for n in
                               mSDUC.N for g in mSDUC.G) +
                           sum(pEmissionCost[g] * mSDUC.vProduct[sc, n, g] * pScenProb[sc] for sc in mSDUC.SC for n in
                               mSDUC.N for g in mSDUC.G) +
                           sum(pInterVarCost[g] * mSDUC.vCommitt[n, g] for n in mSDUC.N for g in mSDUC.G) +
                           sum(pStartupCost[g] * mSDUC.vStartup[n, g] for n in mSDUC.N for g in mSDUC.G) +
                           sum(pShutDownCost[g] * mSDUC.vShutdown[n, g] for n in mSDUC.N for g in mSDUC.G),
                      sense=minimize, doc='total system variable cost     [Meur]')


# constraints
def rule_eBalance(mSDUC, sc, n):
    return sum(mSDUC.vProduct[sc, n, g] for g in mSDUC.G) + mSDUC.vIG[sc, n] + mSDUC.vENS[sc, n] == float(pDemand[n])


mSDUC.eBalance = Constraint(mSDUC.SC, mSDUC.N, rule=rule_eBalance, doc='load generation balance        [GW]')


def rule_eOpReserve(mSDUC, n):
    return sum(pMaxProd[g] * mSDUC.vCommitt[n, g] for g in mSDUC.G) >= float(pOperReserve[n]) + float(pDemand[n])


mSDUC.eOpReserve = Constraint(mSDUC.N, rule=rule_eOpReserve, doc='operating reserve              [GW]')


def rule_eReserveUp(mSDUC, sc, n):
    return sum(pMaxProd[g] * mSDUC.vCommitt[n, g] - mSDUC.vProduct[sc, n, g] for g in mSDUC.G) >= float(
        pOperReserveUp[n])


mSDUC.eReserveUp = Constraint(mSDUC.SC, mSDUC.N, rule=rule_eReserveUp, doc='operating reserve upwards      [GW]')


def rule_eReserveDw(mSDUC, sc, n):
    return sum(pMinProd[g] * mSDUC.vCommitt[n, g] - mSDUC.vProduct[sc, n, g] for g in mSDUC.G) <= -float(
        pOperReserveDw[n])


mSDUC.eReserveDw = Constraint(mSDUC.SC, mSDUC.N, rule=rule_eReserveDw, doc='operating reserve downwards    [GW]')


def rule_eMaxOutput(mSDUC, sc, n, g):
    return mSDUC.vProduct[sc, n, g] / float(pMaxProd[g]) <= mSDUC.vCommitt[n, g]


mSDUC.eMaxOutput = Constraint(mSDUC.SC, mSDUC.N, mSDUC.G, rule=rule_eMaxOutput,
                              doc='max output of a committed unit [GW]')


def rule_eMinOutput(mSDUC, sc, n, g):
    return mSDUC.vProduct[sc, n, g] / float(pMinProd[g]) >= mSDUC.vCommitt[n, g]


mSDUC.eMinOutput = Constraint(mSDUC.SC, mSDUC.N, mSDUC.G, rule=rule_eMinOutput,
                              doc='min output of a committed unit [GW]')


def rule_eTotOutput(mSDUC, sc, n, g):
    return mSDUC.vProduct[sc, n, g] == float(pMinProd[g]) * mSDUC.vCommitt[n, g] + mSDUC.vProduct1[sc, n, g]


mSDUC.eTotOutput = Constraint(mSDUC.SC, mSDUC.N, mSDUC.G, rule=rule_eTotOutput,
                              doc='tot output of a committed unit [GW]')


def rule_eRampUp(mSDUC, sc, n, g):
    return mSDUC.vProduct1[sc, n, g] - max(float(pIniOut[g] - pMinProd[g]), 0) <= float(pRampUp[g]) if n == 0 else \
    mSDUC.vProduct1[sc, n, g] - mSDUC.vProduct1[sc, n - 1, g] <= float(pRampUp[g])


mSDUC.eRampUp = Constraint(mSDUC.SC, mSDUC.N, mSDUC.G, rule=rule_eRampUp, doc='bound on ramp up               [GW]')


def rule_eRampDw(mSDUC, sc, n, g):
    return mSDUC.vProduct1[sc, n, g] - max(float(pIniOut[g] - pMinProd[g]), 0) >= -float(pRampDw[g]) if n == 0 else \
    mSDUC.vProduct1[sc, n, g] - mSDUC.vProduct1[sc, n - 1, g] >= -float(pRampDw[g])


mSDUC.eRampDw = Constraint(mSDUC.SC, mSDUC.N, mSDUC.G, rule=rule_eRampDw, doc='bound on ramp down             [GW]')


def rule_eUCStrShut(mSDUC, n, g):
    return mSDUC.vCommitt[n, g] - float(pIniUC[g]) == mSDUC.vStartup[n, g] - mSDUC.vShutdown[n, g] if n == 0 else \
    mSDUC.vCommitt[n, g] - mSDUC.vCommitt[n - 1, g] == mSDUC.vStartup[n, g] - mSDUC.vShutdown[n, g]


mSDUC.eUCStrShut = Constraint(mSDUC.N, mSDUC.G, rule=rule_eUCStrShut,
                              doc='relation among commitment startup and shutdown')


def rule_eMinTUp(mSDUC, n, g):
    if n >= int(pMinTU[g]) & n <= N:
        return sum(mSDUC.vStartup[nn, g] for nn in range(n + 1 - int(pMinTU[g]), n + 1)) <= mSDUC.vCommitt[n, g]
    else:
        Constraint.Skip()


mSDUC.eMinTUp = Constraint(mSDUC.N, mSDUC.G, rule=rule_eMinTUp, doc='minimum up   time (    committed)')


def rule_eMinTDw(mSDUC, n, g):
    if n >= int(pMinTD[g]) & n <= N:
        return sum(mSDUC.vShutdown[nn, g] for nn in range(n + 1 - int(pMinTD[g]), n + 1)) <= 1 - mSDUC.vCommitt[n, g]
    else:
        Constraint.Skip()


mSDUC.eMinTDw = Constraint(mSDUC.N, mSDUC.G, rule=rule_eMinTDw, doc='minimum down time (not committed)')

mSDUC.write('modelPyomo.lp', io_options={'symbolic_solver_labels': True})  # create file for lp of the mSDUC
solver = SolverFactory('gurobi')  # Select solver
# solver = SolverFactory('glpk')  # Select solver
solver.options['MIPGap'] = 0.0  # Change Mipgap
results = solver.solve(mSDUC, tee=True)  # tee=True displays the output of the solver
mSDUC.solutions.load_from(results)  # necessary for fixing the values of the binary variables

# In case dual variables are required
mSDUC.dual = Suffix(direction=Suffix.IMPORT)
# pSRMC = np.zeros((SC,N))
# for sc in mSDUC.SC:
#    for n in mSDUC.N:
#        print('ppp', sc, n, mSDUC.dual[mSDUC.eBalance[sc,n]])
