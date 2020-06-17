# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import psutil
from   pyomo.opt     import SolverFactory
from   pyomo.environ import Suffix

print('Problem solving         ****')

StartTime = time.time()

#%% solving the problem
Solver = SolverFactory('gurobi')                                                       # select solver
Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
#Solver.options['IISFile'      ] = 'openTEPES_'+CaseName+'.ilp'                        # should be uncommented to show results of IIS
Solver.options['Method'        ] = 2                                                   # barrier method
Solver.options['Presolve'      ] = 2
Solver.options['MIPGap'        ] = 0.03
Solver.options['Threads'       ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)
#Solver.options['TimeLimit'     ] =    40000
#Solver.options['IterationLimit'] = 10000000
if mTEPES.pIndBinGenInvest*len(mTEPES.gc) + mTEPES.pIndBinNetInvest*len(mTEPES.lc) + mTEPES.pIndBinGenOperat() == 0:
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual = Suffix(direction=Suffix.IMPORT)
SolverResults = Solver.solve(mTEPES, tee=True, report_timing=True)                     # tee=True displays the log of the solver
SolverResults.write()                                                                  # summary of the solver results

#%% fix values of binary variables to get dual variables and solve it again
if mTEPES.pIndBinGenInvest*len(mTEPES.gc) + mTEPES.pIndBinNetInvest*len(mTEPES.lc) + mTEPES.pIndBinGenOperat():
    if mTEPES.pIndBinGenInvest*len(mTEPES.gc):
        for gc in mTEPES.gc:
            mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
    if mTEPES.pIndBinNetInvest*len(mTEPES.lc):
        for lc in mTEPES.lc:
            mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
    if mTEPES.pIndBinGenOperat:
        for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
            mTEPES.vCommitment[sc,p,n,t].fix(mTEPES.vCommitment[sc,p,n,t]())
            mTEPES.vStartUp   [sc,p,n,t].fix(mTEPES.vStartUp   [sc,p,n,t]())
            mTEPES.vShutDown  [sc,p,n,t].fix(mTEPES.vShutDown  [sc,p,n,t]())
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual   = Suffix(direction=Suffix.IMPORT)
    SolverResults = Solver.solve(mTEPES, warmstart=True, tee=False, report_timing=True)                                     # tee=True displays the log of the solver
    SolverResults.write()                                                              # summary of the solver results

SolvingTime = time.time() - StartTime
StartTime   = time.time()
print('Solving                               ... ', round(SolvingTime), 's')

print('Total system cost [MEUR]                  ', mTEPES.eTotalTCost.expr())
