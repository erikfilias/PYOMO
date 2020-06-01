# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.0.0 - April 26, 2020

#%% solving the problem
Solver = SolverFactory('gurobi_persistent')                                            # select solver
Solver.set_instance(mTEPES)
Solver.set_gurobi_param('LogFile', 'openTEPES_'+CaseName+'.log')
#Solver.set_gurobi_param('IISFile', 'openTEPES_'+CaseName+'.ilp')                      # introduced to show results of IIS
Solver.set_gurobi_param('Method' , 2   )                                               # barrier method
Solver.set_gurobi_param('MIPGap' , 0.03)
Solver.set_gurobi_param('Threads', round((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2))
#Solver.set_gurobi_param('TimeLimit'     , 40000)
#Solver.set_gurobi_param('IterationLimit', 10000000)
SolverResults = Solver.solve(mTEPES, tee=True)                                         # tee=True displays the log of the solver
SolverResults.write()

#%% fix values of binary variables to get dual variables and solve it again
if pIndBinGenInvest*len(mTEPES.gc) + pIndBinNetInvest*len(mTEPES.lc) + pIndBinGenOperat:
    if pIndBinGenInvest*len(mTEPES.gc):
        for gc in mTEPES.gc:
            mTEPES.vGenerationInvest[gc].setlb(mTEPES.vGenerationInvest[gc]())
            mTEPES.vGenerationInvest[gc].setub(mTEPES.vGenerationInvest[gc]())
            Solver.update_var(mTEPES.vGenerationInvest[gc])
    if pIndBinNetInvest*len(mTEPES.lc):
        for lc in mTEPES.lc:
            mTEPES.vNetworkInvest   [lc].setlb(mTEPES.vNetworkInvest   [lc]())
            mTEPES.vNetworkInvest   [lc].setub(mTEPES.vNetworkInvest   [lc]())
            Solver.update_var(mTEPES.vNetworkInvest   [lc])
    if pIndBinGenOperat:
        for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
            mTEPES.vStartUp [sc,p,n,t].setlb(mTEPES.vStartUp [sc,p,n,t]())
            mTEPES.vStartUp [sc,p,n,t].setub(mTEPES.vStartUp [sc,p,n,t]())
            mTEPES.vShutDown[sc,p,n,t].setlb(mTEPES.vShutDown[sc,p,n,t]())
            mTEPES.vShutDown[sc,p,n,t].setub(mTEPES.vShutDown[sc,p,n,t]())
            Solver.update_var(mTEPES.vCommitment[sc,p,n,t])
            Solver.update_var(mTEPES.vStartUp   [sc,p,n,t])
            Solver.update_var(mTEPES.vShutDown  [sc,p,n,t])
    SolverResults = Solver.solve(mTEPES, tee=True)                                     # tee=True displays the log of the solver
    SolverResults.write()                                                              # summary of the solver results

    # this latest solve is exclusively for obtaining the dual variables
    Solver = SolverFactory('gurobi')                                                   # select solver
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual   = Suffix(direction=Suffix.IMPORT)
    SolverResults = Solver.solve(mTEPES, tee=True)                                     # tee=True displays the log of the solver
    SolverResults.write()                                                              # summary of the solver results
