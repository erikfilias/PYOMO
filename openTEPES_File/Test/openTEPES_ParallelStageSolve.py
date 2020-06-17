# Open Generation and Transmission Operation and Expansion mTEPES.planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import psutil
import pandas as pd
import pyomo.environ as pyo
from   pyomo.environ      import Set, Param, Suffix, TerminationCondition
from   pyomo.opt          import SolverFactory
from   pyomo.opt.parallel import SolverManagerFactory

print('Solving stage in parallel   ****')

StartTime = time.time()

pTotalOCost = 0

if len(mTEPES.gc):
    pGenerationInvestMarginalG = pd.Series([0.]*len(mTEPES.gc ), index=pd.Index                 (mTEPES.gc ))
if len(mTEPES.ec):
    pGenerationInvestMarginalE = pd.Series([0.]*len(mTEPES.ec ), index=pd.Index                 (mTEPES.ec ))
if len(mTEPES.lc):
    pNetworkInvestMarginalCap  = pd.Series([0.]*len(mTEPES.lc ), index=pd.MultiIndex.from_tuples(mTEPES.lc ))
if len(mTEPES.lca):
    pNetworkInvestMarginal2nd  = pd.Series([0.]*len(mTEPES.lca), index=pd.MultiIndex.from_tuples(mTEPES.lca))

# deactivate the original complete o.f. and activate the subproblem o.f.
mTEPES.eTotalTCost.deactivate()
mTEPES.eTotalOCost.activate()

# solver_manager    = SolverManagerFactory('pyro')
# action_handle_map = {}

Solver = SolverFactory('gurobi', symbolic_solver_labels=True)
# Solver.options['LogFile' ] = 'openTEPES_Sbp_'+CaseName+'.log'
Solver.options['Method'  ] = 2
Solver.options['Presolve'] = 2
Solver.options['Threads' ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)

#%% iterative solve for each stage of a year
for st in range(1,int(sum(mTEPES.pDuration.values())/pStageDuration+1)):

    # deactivate all the constraints of the problem that depend on n
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        mTEPES.eInstalGenCap[sc,p,n,gc].deactivate()
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        mTEPES.eInstalGenESS[sc,p,n,ec].deactivate()
        mTEPES.eInstalConESS[sc,p,n,ec].deactivate()
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if mTEPES.pOperReserveUp[sc,p,n]:
            mTEPES.eOperReserveUp[sc,p,n].deactivate()
        if mTEPES.pOperReserveDw[sc,p,n]:
            mTEPES.eOperReserveDw[sc,p,n].deactivate()
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        mTEPES.eBalance[sc,p,n,nd].deactivate()
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % mTEPES.pCycleTimeStep[es] == 0:
            mTEPES.eESSInventory[sc,p,n,es].deactivate()
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if mTEPES.pOperReserveUp[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            mTEPES.eMaxOutput2ndBlock[sc,p,n,nr].deactivate()
        if mTEPES.pOperReserveDw[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            mTEPES.eMinOutput2ndBlock[sc,p,n,nr].deactivate()
        mTEPES.eTotalOutput[sc,p,n,nr].deactivate()
        mTEPES.eUCStrShut[sc,p,n,nr].deactivate()
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t]:
            mTEPES.eRampUp[sc,p,n,t].deactivate()
        if mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t]:
            mTEPES.eRampDw[sc,p,n,t].deactivate()
        if mTEPES.pUpTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pUpTime[t]:
            mTEPES.eMinUpTime[sc,p,n,t].deactivate()
        if mTEPES.pDwTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pDwTime[t]:
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
        if mTEPES.pIndNetLosses:
            mTEPES.eLineLosses1[sc,p,n,ni,nf,cc].deactivate()
            mTEPES.eLineLosses2[sc,p,n,ni,nf,cc].deactivate()

    # remove all the load levels and add only load levels of this stage
    mTEPES.del_component(mTEPES.n )
    mTEPES.del_component(mTEPES.n2)
    mTEPES.n  = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration) and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration)
    mTEPES.n2 = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration) and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration)

    # fix the inventory level at the end of this stage
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]())

    # activate all the constraints of the problem that depend on n but not those that tie with previous or next stage
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        if pyo.value(mTEPES.vGenerationInvest[gc]) > 0.:
            mTEPES.eInstalGenCap[sc,p,n,gc].activate()
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        if pyo.value(mTEPES.vGenerationInvest[ec]) > 0.:
            mTEPES.eInstalGenESS[sc,p,n,ec].activate()
            mTEPES.eInstalConESS[sc,p,n,ec].activate()
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if mTEPES.pOperReserveUp[sc,p,n]:
            mTEPES.eOperReserveUp[sc,p,n].activate()
        if mTEPES.pOperReserveDw[sc,p,n]:
            mTEPES.eOperReserveDw[sc,p,n].activate()
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        mTEPES.eBalance[sc,p,n,nd].activate()
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % mTEPES.pCycleTimeStep[es] == 0:
            mTEPES.eESSInventory[sc,p,n,es].activate()
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if mTEPES.pOperReserveUp[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            mTEPES.eMaxOutput2ndBlock[sc,p,n,nr].activate()
        if mTEPES.pOperReserveDw[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            mTEPES.eMinOutput2ndBlock[sc,p,n,nr].activate()
        mTEPES.eTotalOutput[sc,p,n,nr].activate()
        if st == 1 or (n != mTEPES.n.first() and st > 1):
            mTEPES.eUCStrShut[sc,p,n,nr].activate()
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            mTEPES.eRampUp[sc,p,n,t].activate()
        if mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            mTEPES.eRampDw[sc,p,n,t].activate()
        if mTEPES.pUpTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pUpTime[t]:
            mTEPES.eMinUpTime[sc,p,n,t].activate()
        if mTEPES.pDwTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pDwTime[t]:
            mTEPES.eMinDownTime[sc,p,n,t].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc].activate()
        mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
            mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc].activate()
            mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc].activate()
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if mTEPES.pIndNetLosses:
            mTEPES.eLineLosses1[sc,p,n,ni,nf,cc].activate()
            mTEPES.eLineLosses2[sc,p,n,ni,nf,cc].activate()

    # mTEPES.write('openTEPES_Sbp_'+str(st)+'.lp', io_options={'symbolic_solver_labels': True})

    # send the stage problem to the solver
    # action_handle = solver_manager.queue(mTEPES, opt=Solver, warmstart=False, keepfiles=False, tee=True, load_solutions=False, report_timing=True)
    # action_handle_map[action_handle] = 19591121 + st

    # save the inventory level at the end of this stage for the next one
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        mTEPES.pInitialInventory[es] = mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]()

    StageSolutionTime = time.time() - StartTime
    StartTime         = time.time()
    print('Sending stage    ', st, ' ... ', round(StageSolutionTime), 's')

# retrieve the stage solution from the solver
# for st in range(1,int(sum(mTEPES.pDuration.values())/pStageDuration+1)):
#
#     this_action_handle = solver_manager.wait_any()
#     solved_name        = action_handle_map[this_action_handle]
#     SolverResults      = solver_manager.get_results(this_action_handle)
    # print(solved_name, action_handle_map)
    # SolverResults.write(filename='pp.lp', io_options={'symbolic_solver_labels': True})
    # ObjVal        = SolverResults.Solution.Objective.get('__default_objective__').get('Value')
    # mTEPES.dual = SolverResults.Solution.Constraint()
    # if SolverResults.solver.termination_condition == TerminationCondition.optimal:
    #     SolverResults = Solver.load_vars()
    #     SolverResults = Solver.load_duals()
    #

    SolverResults = Solver.solve(mTEPES, warmstart=False, tee=False, report_timing=False)

    # if SolverResults.solver.termination_condition == TerminationCondition.optimal:
    ObjVal = mTEPES.eTotalOCost.expr()

    # print(mTEPES.eTotalOCost.expr())

    pTotalOCost += ObjVal

    for gc in mTEPES.gc:
        pGenerationInvestMarginalG   [gc      ] = sum(mTEPES.dual[mTEPES.eInstalGenCap       [sc,p,n,gc      ]]                                                             for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n) + pGenerationInvestMarginalG[gc      ]
    for ec in mTEPES.ec:
        pGenerationInvestMarginalE   [ec      ] = sum(mTEPES.dual[mTEPES.eInstalGenESS       [sc,p,n,ec      ]] + mTEPES.rc[mTEPES.eInstalConESS         [sc,p,n,ec      ]] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n) + pGenerationInvestMarginalE[ec      ]
    for ni,nf,cc in mTEPES.lc:
        pNetworkInvestMarginalCap    [ni,nf,cc] = sum(mTEPES.dual[mTEPES.eInstalNetCap1      [sc,p,n,ni,nf,cc]] + mTEPES.dual[mTEPES.eInstalNetCap2      [sc,p,n,ni,nf,cc]] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n) + pNetworkInvestMarginalCap [ni,nf,cc]
    for ni,nf,cc in mTEPES.lca:
        if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
            pNetworkInvestMarginal2nd[ni,nf,cc] = sum(mTEPES.dual[mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc]] + mTEPES.dual[mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc]] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n) + pNetworkInvestMarginal2nd [ni,nf,cc]

    StageSolutionTime = time.time() - StartTime
    StartTime         = time.time()
    print('Retrieving stage ', st, ' ... ', round(StageSolutionTime), 's      Total system variable cost [MEUR] ... ', ObjVal)

mTEPES.pTotalOCost = pTotalOCost

for gc in mTEPES.gc:
    mTEPES.pGenerationInvestMarginalG   [gc      ] = pGenerationInvestMarginalG[gc      ]
for ec in mTEPES.ec:
    mTEPES.pGenerationInvestMarginalE   [ec      ] = pGenerationInvestMarginalE[ec      ]
for ni,nf,cc in mTEPES.lc:
    mTEPES.pNetworkInvestMarginalCap    [ni,nf,cc] = pNetworkInvestMarginalCap [ni,nf,cc]
for ni,nf,cc in mTEPES.lca:
    if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
        mTEPES.pNetworkInvestMarginal2nd[ni,nf,cc] = pNetworkInvestMarginal2nd [ni,nf,cc]
