# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.0.0 - May 7, 2020

pStageDuration = 168

# load the problem into memory
Solver = SolverFactory('gurobi_persistent', report_timing=True, profile_memory=10)
Solver.set_instance(mTEPES, symbolic_solver_labels=True)
Solver.set_gurobi_param('Threads', int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2))
Solver.set_gurobi_param('MIPGap' , 0.03)

# investment decisions are fixed and their costs not taken into account
for gc in mTEPES.gc:
    mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
for lc in mTEPES.lc:
    mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
Solver.remove_constraint(mTEPES.eTotalFCost)
print('Total system fixed cost [MEUR] ... ', mTEPES.vTotalFCost[:]())

MemoryLoadTime = time.time() - StartTime
StartTime      = time.time()
print('Loading problem into memory           ... ', round(SolvingTime), 's')

#%% iterative solve for each stage of a year
for st in range(1,int(pDuration.sum()/pStageDuration+1)):

    # remove all the constraints of the problem that depend on n, all the constraints in the first stage and the constraints from the previous stage after the first one
    Solver.remove_constraint(mTEPES.eTotalVCost)
    Solver.remove_constraint(mTEPES.eTotalECost)
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        Solver.remove_constraint(mTEPES.eInstalGenCap[sc,p,n,gc])
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        Solver.remove_constraint(mTEPES.eInstalGenESS[sc,p,n,ec])
        Solver.remove_constraint(mTEPES.eInstalConESS[sc,p,n,ec])
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if pOperReserveUp[sc,p,n]:
            Solver.remove_constraint(mTEPES.eOperReserveUp[sc,p,n])
        if pOperReserveDw[sc,p,n]:
            Solver.remove_constraint(mTEPES.eOperReserveDw[sc,p,n])
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        Solver.remove_constraint(mTEPES.eBalance[sc,p,n,nd])
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % pESSTimeStep[es] == 0:
            Solver.remove_constraint(mTEPES.eESSInventory[sc,p,n,es])
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]: #and (st == 1 or (n != mTEPES.n.last() and st > 1))
            Solver.remove_constraint(mTEPES.eMaxOutput2ndBlock[sc,p,n,nr])
        if pOperReserveDw[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            Solver.remove_constraint(mTEPES.eMinOutput2ndBlock[sc,p,n,nr])
        Solver.remove_constraint(mTEPES.eTotalOutput[sc,p,n,nr])
    if st == 1 or (n != mTEPES.n.first() and st > 1):
        Solver.remove_constraint(mTEPES.eUCStrShut[sc,p,n,nr])
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.remove_constraint(mTEPES.eRampUp[sc,p,n,t])
        if pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.remove_constraint(mTEPES.eRampDw[sc,p,n,t])
        if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
            Solver.remove_constraint(mTEPES.eMinUpTime[sc,p,n,t])
        if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
            Solver.remove_constraint(mTEPES.eMinDownTime[sc,p,n,t])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        Solver.remove_constraint(mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc])
        Solver.remove_constraint(mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        Solver.remove_constraint(mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        Solver.remove_constraint(mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc])
        Solver.remove_constraint(mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if pIndNetLosses:
            Solver.remove_constraint(mTEPES.eLineLosses1[sc,p,n,ni,nf,cc])
            Solver.remove_constraint(mTEPES.eLineLosses2[sc,p,n,ni,nf,cc])

    # must-run units are fixed to 0 for all the load levels
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pMustRun[nr] == 1:
            mTEPES.vCommitment[sc,p,n,nr].fix(0)

    # remove all the load levels and add only load levels of this stage
    mTEPES.del_component(mTEPES.n )
    mTEPES.del_component(mTEPES.n2)
    mTEPES.n  = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in mTEPES.nn and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration and pDuration[nn] > 0)
    mTEPES.n2 = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in mTEPES.nn and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration and pDuration[nn] > 0)

    # fix the inventory level at the end of this stage
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]())

    # must-run units are fixed to 1 for all the load levels of this stage
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pMustRun[nr] == 1:
            mTEPES.vCommitment[sc,p,n,nr].fix(1)

    # add all the constraints of the problem that depend on n but not those that tie with previous or next stage
    Solver.add_constraint(mTEPES.eTotalVCost)
    Solver.add_constraint(mTEPES.eTotalECost)
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        Solver.add_constraint(mTEPES.eInstalGenCap[sc,p,n,gc])
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        Solver.add_constraint(mTEPES.eInstalGenESS[sc,p,n,ec])
        Solver.add_constraint(mTEPES.eInstalConESS[sc,p,n,ec])
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if pOperReserveUp[sc,p,n]:
            Solver.add_constraint(mTEPES.eOperReserveUp[sc,p,n])
        if pOperReserveDw[sc,p,n]:
            Solver.add_constraint(mTEPES.eOperReserveDw[sc,p,n])
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        Solver.add_constraint(mTEPES.eBalance[sc,p,n,nd])
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % pESSTimeStep[es] == 0:
            Solver.add_constraint(mTEPES.eESSInventory[sc,p,n,es])
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]: #and (st == int(pDuration.sum()/pStageDuration) or (n != mTEPES.n.last() and st < int(pDuration.sum()/pStageDuration)))
            Solver.add_constraint(mTEPES.eMaxOutput2ndBlock[sc,p,n,nr])
        if pOperReserveDw[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            Solver.add_constraint(mTEPES.eMinOutput2ndBlock[sc,p,n,nr])
        Solver.add_constraint(mTEPES.eTotalOutput[sc,p,n,nr])
        if st == 1 or (n != mTEPES.n.first() and st > 1):
            Solver.add_constraint(mTEPES.eUCStrShut[sc,p,n,nr])
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.add_constraint(mTEPES.eRampUp[sc,p,n,t])
        if pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.add_constraint(mTEPES.eRampDw[sc,p,n,t])
        if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
            Solver.add_constraint(mTEPES.eMinUpTime[sc,p,n,t])
        if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
            Solver.add_constraint(mTEPES.eMinDownTime[sc,p,n,t])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        Solver.add_constraint(mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc])
        Solver.add_constraint(mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        Solver.add_constraint(mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        Solver.add_constraint(mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc])
        Solver.add_constraint(mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if pIndNetLosses:
            Solver.add_constraint(mTEPES.eLineLosses1[sc,p,n,ni,nf,cc])
            Solver.add_constraint(mTEPES.eLineLosses2[sc,p,n,ni,nf,cc])

    SolverResults = Solver.solve(mTEPES, warmstart=True, save_results=False, load_solutions=False)
    if SolverResults.solver.termination_condition == TerminationCondition.optimal:
        SolverResults = Solver.load_vars()
        ObjVal   = Solver.get_model_attr('ObjVal'  )
        ObjBound = Solver.get_model_attr('ObjBound')

    StageSolutionTime = time.time() - StartTime
    StartTime         = time.time()
    print('Solving stage    ', st, ' ... ', round(StageSolutionTime), 's      Total system variable cost [MEUR] ... ', ObjVal)
