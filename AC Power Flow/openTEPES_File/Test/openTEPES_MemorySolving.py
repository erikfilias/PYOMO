# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import psutil
from   pyomo.opt     import SolverFactory
from   pyomo.environ import Param

import pandas as pd
dual_eBalance = pd.Series([0.]*len(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd), index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd)))

print('Memory solving   ****')

StartTime = time.time()

#%% solving the problem
Solver = SolverFactory('gurobi_persistent')                                            # select solver
Solver.set_instance(mTEPES)
Solver.set_gurobi_param('LogFile', 'openTEPES_'+CaseName+'.log')
#Solver.set_gurobi_param('IISFile', 'openTEPES_'+CaseName+'.ilp')                      # introduced to show results of IIS
Solver.set_gurobi_param('Kappa'  , 1   )                                               # barrier method
Solver.set_gurobi_param('Method' , 2   )                                               # barrier method
Solver.set_gurobi_param('MIPGap' , 0.03)
Solver.set_gurobi_param('Threads', round((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2))
#Solver.set_gurobi_param('TimeLimit'     , 40000)
#Solver.set_gurobi_param('IterationLimit', 10000000)
if mTEPES.pIndBinGenInvest*len(mTEPES.gc) + mTEPES.pIndBinNetInvest*len(mTEPES.lc) + mTEPES.pIndBinGenOperat() == 0:
    Solver.set_gurobi_param('MIPGap' , 1)
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        dual_eBalance[sc,p,n,nd] = Solver.get_linear_constraint_attr(mTEPES.eBalance[sc,p,n,nd], attr='Pi')
SolverResults = Solver.solve(mTEPES, tee=True, report_timing=True)                     # tee=True displays the log of the solver
SolverResults.write()

#%% fix values of binary variables to get dual variables and solve it again
if mTEPES.pIndBinGenInvest*len(mTEPES.gc) + mTEPES.pIndBinNetInvest*len(mTEPES.lc) + mTEPES.pIndBinGenOperat():
    if mTEPES.pIndBinGenInvest*len(mTEPES.gc):
        for gc in mTEPES.gc:
            mTEPES.vGenerationInvest[gc].setlb(mTEPES.vGenerationInvest[gc]())
            mTEPES.vGenerationInvest[gc].setub(mTEPES.vGenerationInvest[gc]())
            Solver.update_var(mTEPES.vGenerationInvest[gc])
    if mTEPES.pIndBinNetInvest*len(mTEPES.lc):
        for lc in mTEPES.lc:
            mTEPES.vNetworkInvest   [lc].setlb(mTEPES.vNetworkInvest   [lc]())
            mTEPES.vNetworkInvest   [lc].setub(mTEPES.vNetworkInvest   [lc]())
            Solver.update_var(mTEPES.vNetworkInvest   [lc])
    if mTEPES.pIndBinGenOperat:
        for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
            mTEPES.vStartUp [sc,p,n,t].setlb(mTEPES.vStartUp [sc,p,n,t]())
            mTEPES.vStartUp [sc,p,n,t].setub(mTEPES.vStartUp [sc,p,n,t]())
            mTEPES.vShutDown[sc,p,n,t].setlb(mTEPES.vShutDown[sc,p,n,t]())
            mTEPES.vShutDown[sc,p,n,t].setub(mTEPES.vShutDown[sc,p,n,t]())
            Solver.update_var(mTEPES.vCommitment[sc,p,n,t])
            Solver.update_var(mTEPES.vStartUp   [sc,p,n,t])
            Solver.update_var(mTEPES.vShutDown  [sc,p,n,t])
    SolverResults = Solver.solve(mTEPES, tee=True, report_timing=True)                 # tee=True displays the log of the solver
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        dual_eBalance[sc,p,n,nd] = Solver.get_linear_constraint_attr(mTEPES.eBalance[sc,p,n,nd], attr='Pi')
    SolverResults.write()                                                              # summary of the solver results

SolvingTime = time.time() - StartTime
StartTime   = time.time()
print('Solving                               ... ', round(SolvingTime), 's')

print('Total system cost [MEUR]                  ', mTEPES.eTotalTCost.expr())
