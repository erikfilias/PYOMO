# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.0.0 - May 7, 2020

# run these commands in a python terminal, each one individually, to set up the pyro solver server
# pyomo_ns
# dispatch_srvr
# run this command as many times as the number of cores to allow parallel solves
# pyro_mip_server

pStageDuration = 168

# investment decisions are fixed and their costs not taken into account
for gc in mTEPES.gc:
    mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
for lc in mTEPES.lc:
    mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
mTEPES.eTotalFCost.deactivate()
print('Total system fixed cost [MEUR] ... ', sum(mTEPES.vTotalFCost[:]()))

from pyomo.opt.parallel import SolverManagerFactory
solver_manager    = SolverManagerFactory('pyro')
action_handle_map = {} # maps action handles to instances

# load the problem
Solver = SolverFactory('gurobi')
Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
#Solver.options['IISFile'      ] = 'openTEPES_'+CaseName+'.ilp'                        # should be uncommented to show results of IIS
Solver.options['MIPGap'        ] = 0.03
Solver.options['Threads'       ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)

#%% iterative solve for each stage of a year
for st in range(1,int(pDuration.sum()/pStageDuration+1)):

    # deactivate all the constraints of the problem that depend on n
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        mTEPES.eInstalGenCap[sc,p,n,gc].deactivate()
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        mTEPES.eInstalGenESS[sc,p,n,ec].deactivate()
        mTEPES.eInstalConESS[sc,p,n,ec].deactivate()
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if pOperReserveUp[sc,p,n]:
            mTEPES.eOperReserveUp[sc,p,n].deactivate()
        if pOperReserveDw[sc,p,n]:
            mTEPES.eOperReserveDw[sc,p,n].deactivate()
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        mTEPES.eBalance[sc,p,n,nd].deactivate()
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % pESSTimeStep[es] == 0:
            mTEPES.eESSInventory[sc,p,n,es].deactivate()
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            mTEPES.eMaxOutput2ndBlock[sc,p,n,nr].deactivate()
        if pOperReserveDw[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            mTEPES.eMinOutput2ndBlock[sc,p,n,nr].deactivate()
        mTEPES.eTotalOutput[sc,p,n,nr].deactivate()
        mTEPES.eUCStrShut[sc,p,n,nr].deactivate()
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n]:
            mTEPES.eRampUp[sc,p,n,t].deactivate()
        if pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n]:
            mTEPES.eRampDw[sc,p,n,t].deactivate()
        if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
            mTEPES.eMinUpTime[sc,p,n,t].deactivate()
        if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
            mTEPES.eMinDownTime[sc,p,n,t].deactivate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc].deactivate()
        mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc].deactivate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc].deactivate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc].deactivate()
        mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc].deactivate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if pIndNetLosses:
            mTEPES.eLineLosses1[sc,p,n,ni,nf,cc].deactivate()
            mTEPES.eLineLosses2[sc,p,n,ni,nf,cc].deactivate()

    # remove all the load levels and add only load levels of this stage
    mTEPES.del_component(mTEPES.n )
    mTEPES.del_component(mTEPES.n2)
    mTEPES.n  = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in mTEPES.nn and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration and pDuration[nn] > 0)
    mTEPES.n2 = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in mTEPES.nn and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration and pDuration[nn] > 0)

    # fix the inventory level at the end of this stage
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]())

    # activate all the constraints of the problem that depend on n but not those that tie with previous or next stage
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        mTEPES.eInstalGenCap[sc,p,n,gc].activate()
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        mTEPES.eInstalGenESS[sc,p,n,ec].activate()
        mTEPES.eInstalConESS[sc,p,n,ec].activate()
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if pOperReserveUp[sc,p,n]:
            mTEPES.eOperReserveUp[sc,p,n].activate()
        if pOperReserveDw[sc,p,n]:
            mTEPES.eOperReserveDw[sc,p,n].activate()
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        mTEPES.eBalance[sc,p,n,nd].activate()
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % pESSTimeStep[es] == 0:
            mTEPES.eESSInventory[sc,p,n,es].activate()
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            mTEPES.eMaxOutput2ndBlock[sc,p,n,nr].activate()
        if pOperReserveDw[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
            mTEPES.eMinOutput2ndBlock[sc,p,n,nr].activate()
        mTEPES.eTotalOutput[sc,p,n,nr].activate()
        if st == 1 or (n != mTEPES.n.first() and st > 1):
            mTEPES.eUCStrShut[sc,p,n,nr].activate()
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            mTEPES.eRampUp[sc,p,n,t].activate()
        if pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            mTEPES.eRampDw[sc,p,n,t].activate()
        if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
            mTEPES.eMinUpTime[sc,p,n,t].activate()
        if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
            mTEPES.eMinDownTime[sc,p,n,t].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc].activate()
        mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc].activate()
        mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if pIndNetLosses:
            mTEPES.eLineLosses1[sc,p,n,ni,nf,cc].activate()
            mTEPES.eLineLosses2[sc,p,n,ni,nf,cc].activate()

    # mTEPES.write('openTEPES_HY_'+str(st)+'.lp', io_options={'symbolic_solver_labels': True})

    # send the stage problem to the solver
    action_handle = solver_manager.queue(mTEPES, opt=Solver, warmstart=False, tee=True)
    action_handle_map[action_handle] = 19591121 + st

    # save the inventory level at the end of this stage for the next one
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        pESSInitialInventory[es] = mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]()

    StageSolutionTime = time.time() - StartTime
    StartTime         = time.time()
    print('Sending stage    ', st, ' ... ', round(StageSolutionTime), 's')

# retrieve the stage solution from the solver
for st in range(1,int(pDuration.sum()/pStageDuration+1)):

    this_action_handle = solver_manager.wait_any()
    solved_name        = action_handle_map[this_action_handle]
    SolverResults      = solver_manager.get_results(this_action_handle)

    StageSolutionTime = time.time() - StartTime
    StartTime         = time.time()
    print('Retrieving stage ', st, ' ... ', round(StageSolutionTime), 's      Total system variable cost [MEUR] ', mTEPES.eTotalTCost.expr())
