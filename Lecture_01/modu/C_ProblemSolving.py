# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.7.26 - December 28, 2020

import time
import psutil
from   pyomo.opt     import SolverFactory
from   pyomo.environ import Suffix

def ProblemSolving(CaseName, SolverName, model):
    print('Problem solving             ****')

    StartTime = time.time()

    #%% solving the problem
    Solver = SolverFactory(SolverName)                                                       # select solver
    if SolverName == 'gurobi':
        Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
        # Solver.options['IISFile'     ] = CaseName+'/openTEPES_'+CaseName+'.ilp'              # should be uncommented to show results of IIS
        Solver.options['Method'        ] = 2                                                   # barrier method
        Solver.options['Presolve'      ] = 2
        Solver.options['Crossover'     ] = 0
        Solver.options['MIPGap'        ] = 0.025
        Solver.options['Threads'       ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)
        Solver.options['TimeLimit'     ] =    7200
        Solver.options['IterationLimit'] = 7200000
        # Solver.options['relax_integrality'] =  1                                       # introduced to show results of the dual variables
        # Solver.options['Crossover'        ] = -1
        model.dual = Suffix(direction=Suffix.IMPORT)
    SolverResults = Solver.solve(model, tee=True, report_timing=True)                     # tee=True displays the log of the solver
    SolverResults.write()                                                                  # summary of the solver results


    SolvingTime = time.time() - StartTime
    StartTime   = time.time()
    print('Solving                               ... ', round(SolvingTime), 's')

    print('Total system cost [M$]                  ', model.eTotalTCost.expr())
