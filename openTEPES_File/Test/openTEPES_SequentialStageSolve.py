# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import psutil        # access the number of CPUs
import pandas as pd
from   pyomo.environ import Set, TerminationCondition
from   pyomo.opt     import SolverFactory

print('Solving stages sequentially ****')

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

StartTime = time.time()

# load the problem into memory
Solver = SolverFactory('gurobi_persistent', report_timing=True, profile_memory=10)
Solver.set_instance(mTEPES, symbolic_solver_labels=True)
# Solver.set_gurobi_param('LogFile' , 'openTEPES_Sbp_'+CaseName+'.log')
Solver.set_gurobi_param('Method'  , 2)
Solver.set_gurobi_param('Presolve', 2)
Solver.set_gurobi_param('Threads' , int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2))

MemoryLoadTime = time.time() - StartTime
StartTime      = time.time()
print('Loading problem into memory           ... ', round(MemoryLoadTime), 's')

Solver.remove_constraint(mTEPES.eTotalFCost)

#%% iterative solve for each stage of a year
for st in range(1,int(sum(mTEPES.pDuration.values())/pStageDuration+1)):

    # remove all the constraints of the problem that depend on n, all the constraints in the first stage and the constraints from the previous stage after the first one
    Solver.remove_constraint(mTEPES.eTotalVCost)
    Solver.remove_constraint(mTEPES.eTotalECost)
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        Solver.remove_constraint(mTEPES.eInstalGenCap[sc,p,n,gc])
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        Solver.remove_constraint(mTEPES.eInstalGenESS[sc,p,n,ec])
        Solver.remove_constraint(mTEPES.eInstalConESS[sc,p,n,ec])
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if mTEPES.pOperReserveUp[sc,p,n]:
            Solver.remove_constraint(mTEPES.eOperReserveUp[sc,p,n])
        if mTEPES.pOperReserveDw[sc,p,n]:
            Solver.remove_constraint(mTEPES.eOperReserveDw[sc,p,n])
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        Solver.remove_constraint(mTEPES.eBalance[sc,p,n,nd])
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % mTEPES.pCycleTimeStep[es] == 0:
            Solver.remove_constraint(mTEPES.eESSInventory[sc,p,n,es])
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if mTEPES.pOperReserveUp[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]: #and (st == 1 or (n != mTEPES.n.last() and st > 1))
            Solver.remove_constraint(mTEPES.eMaxOutput2ndBlock[sc,p,n,nr])
        if mTEPES.pOperReserveDw[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            Solver.remove_constraint(mTEPES.eMinOutput2ndBlock[sc,p,n,nr])
        Solver.remove_constraint(mTEPES.eTotalOutput[sc,p,n,nr])
    if st == 1 or (n != mTEPES.n.first() and st > 1):
        Solver.remove_constraint(mTEPES.eUCStrShut[sc,p,n,nr])
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.remove_constraint(mTEPES.eRampUp[sc,p,n,t])
        if mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.remove_constraint(mTEPES.eRampDw[sc,p,n,t])
        if mTEPES.pUpTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pUpTime[t]:
            Solver.remove_constraint(mTEPES.eMinUpTime[sc,p,n,t])
        if mTEPES.pDwTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pDwTime[t]:
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
        if mTEPES.pIndNetLosses:
            Solver.remove_constraint(mTEPES.eLineLosses1[sc,p,n,ni,nf,cc])
            Solver.remove_constraint(mTEPES.eLineLosses2[sc,p,n,ni,nf,cc])

    # remove all the load levels and add only load levels of this stage
    mTEPES.del_component(mTEPES.n )
    mTEPES.del_component(mTEPES.n2)
    mTEPES.n  = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration) and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration)
    mTEPES.n2 = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration) and mTEPES.nn.ord(nn) > (st-1)*pStageDuration and mTEPES.nn.ord(nn) <= st*pStageDuration)

    # fix the inventory level at the end of this stage
    for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
        mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es]())

    # add all the constraints of the problem that depend on n but not those that tie with previous or next stage
    Solver.add_constraint(mTEPES.eTotalVCost)
    Solver.add_constraint(mTEPES.eTotalECost)
    for sc,p,n,gc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gc:
        Solver.add_constraint(mTEPES.eInstalGenCap[sc,p,n,gc])
    for sc,p,n,ec in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ec:
        Solver.add_constraint(mTEPES.eInstalGenESS[sc,p,n,ec])
        Solver.add_constraint(mTEPES.eInstalConESS[sc,p,n,ec])
    for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
        if mTEPES.pOperReserveUp[sc,p,n]:
            Solver.add_constraint(mTEPES.eOperReserveUp[sc,p,n])
        if mTEPES.pOperReserveDw[sc,p,n]:
            Solver.add_constraint(mTEPES.eOperReserveDw[sc,p,n])
    for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd:
        Solver.add_constraint(mTEPES.eBalance[sc,p,n,nd])
    for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
        if mTEPES.n.ord(n) % mTEPES.pCycleTimeStep[es] == 0:
            Solver.add_constraint(mTEPES.eESSInventory[sc,p,n,es])
    for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
        if mTEPES.pOperReserveUp[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]: #and (st == int(pDuration.sum()/pStageDuration) or (n != mTEPES.n.last() and st < int(pDuration.sum()/pStageDuration)))
            Solver.add_constraint(mTEPES.eMaxOutput2ndBlock[sc,p,n,nr])
        if mTEPES.pOperReserveDw[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
            Solver.add_constraint(mTEPES.eMinOutput2ndBlock[sc,p,n,nr])
        Solver.add_constraint(mTEPES.eTotalOutput[sc,p,n,nr])
        if st == 1 or (n != mTEPES.n.first() and st > 1):
            Solver.add_constraint(mTEPES.eUCStrShut[sc,p,n,nr])
    for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
        if mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.add_constraint(mTEPES.eRampUp[sc,p,n,t])
        if mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and (st == 1 or (n != mTEPES.n.first() and st > 1)):
            Solver.add_constraint(mTEPES.eRampDw[sc,p,n,t])
        if mTEPES.pUpTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pUpTime[t]:
            Solver.add_constraint(mTEPES.eMinUpTime[sc,p,n,t])
        if mTEPES.pDwTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pDwTime[t]:
            Solver.add_constraint(mTEPES.eMinDownTime[sc,p,n,t])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lc:
        Solver.add_constraint(mTEPES.eInstalNetCap1[sc,p,n,ni,nf,cc])
        Solver.add_constraint(mTEPES.eInstalNetCap2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lea:
        Solver.add_constraint(mTEPES.eKirchhoff2ndLawExst[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.lca:
        if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
            mTEPES.vNetworkInvest.display()
            Solver.add_constraint(mTEPES.eKirchhoff2ndLawCnd1[sc,p,n,ni,nf,cc])
            Solver.add_constraint(mTEPES.eKirchhoff2ndLawCnd2[sc,p,n,ni,nf,cc])
    for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.ll:
        if mTEPES.pIndNetLosses:
            Solver.add_constraint(mTEPES.eLineLosses1[sc,p,n,ni,nf,cc])
            Solver.add_constraint(mTEPES.eLineLosses2[sc,p,n,ni,nf,cc])

    SolverResults = Solver.solve(mTEPES, warmstart=False, keepfiles=False, tee=False, load_solutions=False, save_results=False, report_timing=False)
    if SolverResults.solver.termination_condition == TerminationCondition.optimal:
        SolverResults = Solver.load_vars()
        SolverResults = Solver.load_duals()
        ObjVal        = Solver.get_model_attr('ObjVal')
    else:
        print('Subproblem infeasible')

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
    print('Solving stage    ', st, ' ... ', round(StageSolutionTime), 's      Total system variable cost [MEUR] ... ', ObjVal)

mTEPES.pTotalOCost = pTotalOCost

mTEPES.del_component(mTEPES.n )
mTEPES.del_component(mTEPES.n2)
mTEPES.n  = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration))
mTEPES.n2 = Set(initialize=mTEPES.nn, ordered=True, doc='load levels', filter=lambda mTEPES,nn: nn in list(mTEPES.pDuration))

for gc in mTEPES.gc:
    mTEPES.pGenerationInvestMarginalG   [gc      ] = pGenerationInvestMarginalG[gc      ]
for ec in mTEPES.ec:
    mTEPES.pGenerationInvestMarginalE   [ec      ] = pGenerationInvestMarginalE[ec      ]
for ni,nf,cc in mTEPES.lc:
    mTEPES.pNetworkInvestMarginalCap    [ni,nf,cc] = pNetworkInvestMarginalCap [ni,nf,cc]
for ni,nf,cc in mTEPES.lca:
    if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
        mTEPES.pNetworkInvestMarginal2nd[ni,nf,cc] = pNetworkInvestMarginal2nd [ni,nf,cc]
