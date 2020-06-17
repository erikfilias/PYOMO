# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import builtins
import psutil
import importlib
import pandas as pd
import pyomo.environ as pyo
from   pyomo.environ      import ConcreteModel, Set, RangeSet, Param, Var, NonNegativeReals, Binary, UnitInterval, RealSet, Constraint, Objective, minimize, Suffix
from   pyomo.opt          import SolverFactory

print('Benders decomposition       ****')

StartTime = time.time()

# solve the stage problem in memory (sequential) or writing the lp file (parallel)
pIndSequentialSolving = 1
# duration of the stage (weekly -168- or monthly -672- or trimonthly -2184- is what makes sense from a system operation point of view
# this value must be larger or equal than the shortest duration of any storage type (e.g., weekly)
builtins.pStageDuration = 4368

# Benders convergence tolerance
BdTol   = 1e-6

mMaster    = ConcreteModel('Master problem')
# maximum number of Benders iterations
mMaster.l  = RangeSet(12)
mMaster.ll = Set(doc='iterations')

mMaster.vTotalFCost = Var(initialize=0., within=NonNegativeReals, doc='total system fixed    cost                     [MEUR]')

if mTEPES.pIndBinGenInvest == 0:
    mMaster.vGenerationInvest = Var(mTEPES.gc, within=UnitInterval, doc='generation investment decision exists in a year [0,1]')
else:
    mMaster.vGenerationInvest = Var(mTEPES.gc, within=Binary,       doc='generation investment decision exists in a year {0,1}')

if mTEPES.pIndBinNetInvest == 0:
    mMaster.vNetworkInvest    = Var(mTEPES.lc, within=UnitInterval, doc='network    investment decision exists in a year [0,1]')
else:
    mMaster.vNetworkInvest    = Var(mTEPES.lc, within=Binary,       doc='network    investment decision exists in a year {0,1}')

mMaster.vTheta  = Var(initialize=0., doc='total system variable cost                     [MEUR]', within=RealSet)

def eTotalFCost(mMaster):
   return mMaster.vTotalFCost == sum(mTEPES.pGenFixedCost[gc] * mMaster.vGenerationInvest[gc] for gc in mTEPES.gc) + sum(mTEPES.pNetFixedCost[lc] * mMaster.vNetworkInvest[lc] for lc in mTEPES.lc)
mMaster.eTotalFCost = Constraint(rule=eTotalFCost, doc='total system fixed    cost [MEUR]')

def eTotalTCostMst(mMaster):
    return mMaster.vTotalFCost + mMaster.vTheta
mMaster.eTotalTCostMst = Objective(rule=eTotalTCostMst, sense=minimize, doc='total system cost [MEUR]')

def eBd_Cuts(mMaster, ll):
    return mMaster.vTheta - pTotalOCost_L[ll] >= - sum(pGenerationInvestMarginalG_L[ll,gc      ] * (pGenerationInvest[ll,gc      ] - mMaster.vGenerationInvest[gc      ]) for gc       in mTEPES.gc ) \
                                                 - sum(pGenerationInvestMarginalE_L[ll,ec      ] * (pGenerationInvest[ll,ec      ] - mMaster.vGenerationInvest[ec      ]) for ec       in mTEPES.ec ) \
                                                 - sum(pNetworkInvestMarginalCap_L [ll,ni,nf,cc] * (pNetworkInvest   [ll,ni,nf,cc] - mMaster.vNetworkInvest   [ni,nf,cc]) for ni,nf,cc in mTEPES.lc ) \
                                                 + sum(pNetworkInvestMarginal2nd_L [ll,ni,nf,cc] * (pNetworkInvest   [ll,ni,nf,cc] - mMaster.vNetworkInvest   [ni,nf,cc]) for ni,nf,cc in mTEPES.lca) \

def eTotalOCost(mTEPES):
    return mTEPES.vTotalVCost + mTEPES.vTotalECost
mTEPES.eTotalOCost = Objective(rule=eTotalOCost, sense=minimize, doc='total system cost [MEUR]')

mTEPES.dual                       = Suffix(direction=Suffix.IMPORT)
mTEPES.pTotalOCost                = Param(initialize=0.,  doc='total system operation cost                    [MEUR]', mutable=True)
mTEPES.pGenerationInvestMarginalG = Param(mTEPES.gc ,     doc='marginal of generation investment decision     [MEUR]', mutable=True)
mTEPES.pGenerationInvestMarginalE = Param(mTEPES.ec ,     doc='marginal of generation investment decision     [MEUR]', mutable=True)
mTEPES.pNetworkInvestMarginalCap  = Param(mTEPES.lc ,     doc='marginal of network    investment decision     [MEUR]', mutable=True)
mTEPES.pNetworkInvestMarginal2nd  = Param(mTEPES.lca,     doc='marginal of network    investment decision     [MEUR]', mutable=True)

Solver = SolverFactory('gurobi')
# Solver.options['LogFile' ] = 'openTEPES_Mst_'+CaseName+'.log'
Solver.options['MIPGap'  ] = 0.03
Solver.options['Threads' ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)

# initialization
Z_Lower = float('-inf')
Z_Upper = float(' inf')

pTotalOCost_L = pd.Series([0.]*len(mMaster.l), index=mMaster.l)

if len(mTEPES.gc):
    pGenerationInvest            = pd.Series([0.]*len(mMaster.l*mTEPES.gc ), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.gc ))
    pGenerationInvestMarginalG_L = pd.Series([0.]*len(mMaster.l*mTEPES.gc ), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.gc ))
if len(mTEPES.ec):
    pGenerationInvestMarginalE_L = pd.Series([0.]*len(mMaster.l*mTEPES.ec ), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.ec ))
if len(mTEPES.lc):
    pNetworkInvest               = pd.Series([0.]*len(mMaster.l*mTEPES.lc ), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.lc ))
    pNetworkInvestMarginalCap_L  = pd.Series([0.]*len(mMaster.l*mTEPES.lc ), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.lc ))
if len(mTEPES.lca):
    pNetworkInvestMarginal2nd_L  = pd.Series([0.]*len(mMaster.l*mTEPES.lca), index=pd.MultiIndex.from_tuples(mMaster.l*mTEPES.lca))

# Benders algorithm
mMaster.vTheta.fix(0)
for l in mMaster.l:
    if (abs(1-Z_Lower/Z_Upper) > BdTol or l == 1) and Z_Lower < Z_Upper:

        # solving master problem
        # mMaster.write('openTEPES_Mst_'+str(l)+'_'+CaseName+'.lp', io_options={'symbolic_solver_labels': True})  # create lp-format file
        SolverResultsMst = Solver.solve(mMaster)
        pTotalTCost      = mMaster.eTotalTCostMst.expr()

        if l == 1:
            mMaster.vTheta.free()

        # storing the master solution and fixing investment decision for the subproblem
        for gc in mTEPES.gc:
            pGenerationInvest[l,gc] = pyo.value(mMaster.vGenerationInvest[gc])
            mTEPES.vGenerationInvest[gc].fix(pGenerationInvest_L[l,gc])
        for gc in mTEPES.ec:
            pGenerationInvest[l,ec] = pyo.value(mMaster.vGenerationInvest[ec])
            mTEPES.vGenerationInvest[ec].fix(pGenerationInvest_L[l,ec])
        for ni,nf,cc in mTEPES.lc:
            pNetworkInvest[l,ni,nf,cc] = pyo.value(mMaster.vNetworkInvest[ni,nf,cc])
            mTEPES.vNetworkInvest[ni,nf,cc].fix(pNetworkInvest[l,ni,nf,cc])

        # solving subproblem
        if pIndSequentialSolving == 1:
            if l == 1:
                import           openTEPES_SequentialStageSolve
            else:
                importlib.reload(openTEPES_SequentialStageSolve)
        else:
            # import os
            # from multiprocessing import Process, freeze_support, set_start_method
            # def ns():
            #     # run these commands to set up the pyro solver server
            #     os.system('pyomo_ns')
            #     time.sleep(2)
            #     os.system('dispatch_srvr')
            #     time.sleep(2)
            #     # run this command as many times as the number of cores to allow parallel solves
            #     os.system('pyro_mip_server')
            #     time.sleep(2)
            #     os.system('pyro_mip_server')
            #     time.sleep(2)
            #     os.system('pyro_mip_server')
            #     time.sleep(2)
            #     os.system('pyro_mip_server')
            # freeze_support()
            # set_start_method('spawn')
            # p = Process(target=ns)
            # p.start()
            if l == 1:
                import           openTEPES_ParallelStageSolve
            else:
                importlib.reload(openTEPES_ParallelStageSolve)

        pTotalOCost      = pyo.value(mTEPES.pTotalOCost)
        pTotalOCost_L[l] = pTotalOCost

        # updating lower and upper bound
        Z_Lower =              pTotalTCost
        Z_Upper = min(Z_Upper, pTotalTCost - pyo.value(mMaster.vTheta) + pTotalOCost)
        print('Iteration ', l, 'Z_Lower ... ', Z_Lower)
        print('Iteration ', l, 'Z_Upper ... ', Z_Upper)

        for gc in mTEPES.gc:
            pGenerationInvestMarginalG_L[l,gc      ] =    mTEPES.pGenerationInvestMarginalG[gc      ]()
        for ec in mTEPES.ec:
            pGenerationInvestMarginalE_L[l,ec      ] =    mTEPES.pGenerationInvestMarginalE[ec      ]()
        for ni,nf,cc in mTEPES.lc:
            pNetworkInvestMarginalCap_L [l,ni,nf,cc] =    mTEPES.pNetworkInvestMarginalCap [ni,nf,cc]()
        for ni,nf,cc in mTEPES.lca:
            if pyo.value(mTEPES.vNetworkInvest[ni,nf,cc]) > 0.:
                pNetworkInvestMarginal2nd_L[l,ni,nf,cc] = mTEPES.pNetworkInvestMarginal2nd [ni,nf,cc]()

        # delete Benders cuts because they are regenerated in each iteration
        if l > 1:
            mMaster.del_component(mMaster.eBd_Cuts)
        # add one cut
        mMaster.ll.add(l)
        ll = mMaster.ll
        mMaster.eBd_Cuts = Constraint(mMaster.ll, rule=eBd_Cuts, doc='Benders cuts')
        # Benders cuts are considered as lazy
        # for ll in mMaster.ll:
        #     Solver.set_linear_constraint_attr(mMaster.eBd_Cuts[ll], 'Lazy', 1)

SolvingTime = time.time() - StartTime
StartTime   = time.time()
print('Benders decomposition                 ... ', round(SolvingTime), 's')

print('Total system cost [MEUR]                  ', Z_Upper)
