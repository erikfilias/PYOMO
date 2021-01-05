Solver = SolverFactory('gurobi_persistent')
Solver.set_instance(mTEPES)
Solver.set_gurobi_param('LazyConstraints', 1)

for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
    if   pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n]:
        Solver.set_linear_constraint_attr(mTEPES.eRampUp     [sc,p,n,t], 'Lazy', 1)
    if   pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n]:
        Solver.set_linear_constraint_attr(mTEPES.eRampDw     [sc,p,n,t], 'Lazy', 1)
    if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
        Solver.set_linear_constraint_attr(mTEPES.eMinUpTime  [sc,p,n,t], 'Lazy', 1)
    if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
        Solver.set_linear_constraint_attr(mTEPES.eMinDownTime[sc,p,n,t], 'Lazy', 1)

Solver.set_gurobi_param('LogFile', 'openTEPES_'+CaseName+'.log')
Solver.set_gurobi_param('Method' , 2   )                                               # barrier method
Solver.set_gurobi_param('Threads', round((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2))
Solver.set_gurobi_param('MIPGap' , 0.03)
SolverResults = Solver.solve(mTEPES, tee=True)                                         # tee=True displays the output of the solver
SolverResults.write()

Solver.set_gurobi_param('MIPGap', 0.1)
SolverResults = Solver.solve(mTEPES, tee=True)
SolverResults.write()
