# %% openTEPES model with Linear AC load flow constraints
# %% libraries
import pandas        as pd
import time          # count clock time
import psutil        # access the number of CPUs
import pyomo.environ as pyo
from   pyomo.environ import Set, Var, Param, Binary, NonNegativeReals, RealSet, UnitInterval, Constraint, ConcreteModel, Objective, minimize, Suffix, DataPortal, TerminationCondition
from   pyomo.opt     import SolverFactory
from   collections   import defaultdict
import openTEPES_InputData as oT_ID
import builtins

# %% Model Formulation
# mTEPES.L         = RangeSet(20)

# #PARAM BAR
# mTEPES.gshb      = Param(model.BAR)
# mTEPES.bshb      = Param(model.BAR)


# #PARAM BRANCH
# mTEPES.a         = Param(model.RAM)
# mTEPES.r         = Param(model.RAM)
# mTEPES.x         = Param(model.RAM)
# mTEPES.z2        = Param(model.RAM, initialize = 0)
# mTEPES.fi        = Param(model.RAM)
# mTEPES.smax      = Param(model.RAM)
# mTEPES.bsh       = Param(model.RAM)

# mTEPES.bsh_half  = Param(model.RAM, initialize = model.bsh)
# mTEPES.fi_rad    = Param(model.RAM, initialize = model.fi)
# mTEPES.smax_pu   = Param(model.RAM, initialize = model.smax)
# mTEPES.a_init    = Param(model.RAM, initialize = model.a)

# #create PARAM LINEARIZATION
# mTEPES.delta_S   = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,                           doc='Delta Apparent Power Flow                                                   [GVA]')
# mTEPES.m         = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L,                 doc='Delta Active and Reactive Power Flow                                        [GVA]')

L                              = range(20)
pLineR                         = oT_ID.pLineX * 0 + 0.01
pLineBsh                       = oT_ID.pLineX * 0
pLineTAP                       = oT_ID.pLineX * 0 + 1
pLineFi                        = oT_ID.pLineX * 0
pLineSmax                      = oT_ID.pLineNTC * 1.5

pLineZ2                        = pLineR**2 + oT_ID.pLineX**2
pLineBsh                       = pLineBsh/2
pLineTAP                       = 1/pLineTAP
pLineFi                        = (pLineFi*3.14159265359)/180

pLineDelta_S                   = pLineSmax/len(L)

for ni,nf,cc in mTEPES.la:
    for k in L:
        pLineM[ni, nf, cc, k] = (2*k-1)*(pLineDelta_S[ni, nf, cc])

#create VAR LINEARIZATION
mTEPES.vDelta_P                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L,                   doc='Delta Active Power Flow                                        [GVA]')
mTEPES.vDelta_Q                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L,                   doc='Delta Reactive Power Flow                                      [GVA]')
mTEPES.vP_max                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,                             doc='Maximum bound of Active Power Flow                             [GVA]')
mTEPES.vP_min                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,                             doc='Minimum bound of Active Power Flow                             [GVA]')
mTEPES.vQ_max                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,                             doc='Maximum bound of Reactive Power Flow                           [GVA]')
mTEPES.vQ_min                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,                             doc='Minimum bound of Reactive Power Flow                           [GVA]')

#VARIABLES
mTEPES.vVoltage_sqr            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,nd:(Vmin**2,Vmax**2),          doc='voltage magnitude                                   [p.u.]')
mTEPES.vCurrentFlow_sqr        = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,*la:(0,pLineNTC[la]),          doc='Current flow                                          [A]')
mTEPES.vP                      = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,             doc='Active Power Flow                                               [GW]')
mTEPES.vQ                      = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,             doc='Reactive Power Flow                                             [GW]')
mTEPES.vReactiveTotalOutput    = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g , within=RealSet,             doc='Reactive total output of the unit                             [GVAr]')

def eTotalFCost(mTEPES):
    return mTEPES.vTotalFCost == sum(oT_ID.pNetFixedCost[lc] * mTEPES.vNetworkInvest   [lc] for lc in mTEPES.lc) + sum(oT_ID.pGenFixedCost[gc] * mTEPES.vGenerationInvest[gc] for gc in mTEPES.gc)
mTEPES.eTotalFCost = Constraint(rule=eTotalFCost, doc='total system fixed    cost [MEUR]')

def eTotalVCost(mTEPES):
    return mTEPES.vTotalVCost == (sum(oT_ID.pScenProb[sc] * oT_ID.pDuration[n] * oT_ID.pENSCost             * mTEPES.vENS        [sc,p,n,nd] for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd) +
                                  sum(oT_ID.pScenProb[sc] * oT_ID.pDuration[n] * oT_ID.pLinearVarCost  [nr] * mTEPES.vTotalOutput[sc,p,n,nr]                                                         +
                                      oT_ID.pScenProb[sc] * oT_ID.pDuration[n] * oT_ID.pConstantVarCost[nr] * mTEPES.vCommitment [sc,p,n,nr]                                                         +
                                      oT_ID.pScenProb[sc]                * oT_ID.pStartUpCost    [nr] * mTEPES.vStartUp    [sc,p,n,nr]                                                         +
                                      oT_ID.pScenProb[sc]                * oT_ID.pShutDownCost   [nr] * mTEPES.vShutDown   [sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr) )
mTEPES.eTotalVCost = Constraint(rule=eTotalVCost, doc='total system variable cost [MEUR]')

def eTotalECost(mTEPES):
    return mTEPES.vTotalECost == sum(oT_ID.pScenProb[sc] * oT_ID.pCO2Cost * oT_ID.pCO2EmissionRate[nr] * mTEPES.vTotalOutput[sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)
mTEPES.eTotalECost = Constraint(rule=eTotalECost, doc='total system emission cost [MEUR]')

def eTotalTCost(mTEPES):
    return mTEPES.vTotalFCost + mTEPES.vTotalVCost + mTEPES.vTotalECost
mTEPES.eTotalTCost = Objective(rule=eTotalTCost, sense=minimize, doc='total system cost [MEUR]')

GeneratingOFTime = time.time() - StartTime
StartTime        = time.time()
print('Generating objective function         ... ', round(GeneratingOFTime), 's')

#%% constraints
def eInstalGenCap(mTEPES,sc,p,n,gc):
    return mTEPES.vCommitment[sc,p,n,gc] <= mTEPES.vGenerationInvest[gc]
mTEPES.eInstalGenCap = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gc, rule=eInstalGenCap, doc='commitment  if installed unit [p.u.]')

print('eInstalGenCap         ... ', len(mTEPES.eInstalGenCap), ' rows')

def eInstalGenESS(mTEPES,sc,p,n,ec):
    return mTEPES.vTotalOutput[sc,p,n,ec] / oT_ID.pMaxPower[sc,p,n,ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalGenESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalGenESS, doc='output      if installed ESS unit [p.u.]')

print('eInstalGenESS         ... ', len(mTEPES.eInstalGenESS), ' rows')

def eInstalConESS(mTEPES,sc,p,n,ec):
    return mTEPES.vTotalOutput[sc,p,n,ec] / oT_ID.pMaxCharge      [ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalConESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalConESS, doc='consumption if installed ESS unit [p.u.]')

print('eInstalConESS         ... ', len(mTEPES.eInstalConESS), ' rows')

#%%
def eOperReserveUp(mTEPES,sc,p,n):
    if oT_ID.pOperReserveUp[sc,p,n]:
        return sum(mTEPES.vReserveUp  [sc,p,n,nr] for nr in mTEPES.nr) >= oT_ID.pOperReserveUp[sc,p,n]
mTEPES.eOperReserveUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveUp, doc='up   operating reserve [GW]')

print('eOperReserveUp        ... ', len(mTEPES.eOperReserveUp), ' rows')

def eOperReserveDw(mTEPES,sc,p,n):
    if oT_ID.pOperReserveDw[sc,p,n]:
        return sum(mTEPES.vReserveDown[sc,p,n,nr] for nr in mTEPES.nr) >= oT_ID.pOperReserveDw[sc,p,n]
mTEPES.eOperReserveDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveDw, doc='down operating reserve [GW]')

print('eOperReserveDw        ... ', len(mTEPES.eOperReserveDw), ' rows')

def eBalance(mTEPES,sc,p,n,nd):
    return (sum(mTEPES.vTotalOutput[sc,p,n,g] for g in oT_ID.pNode2Gen[nd]) - sum(mTEPES.vESSCharge[sc,p,n,es] for es in oT_ID.pNode2ESS[nd]) + mTEPES.vENS[sc,p,n,nd] == oT_ID.pDemand[nd][sc,p,n] +
        sum(mTEPES.vLineLosses[sc,p,n,nd,lout ] for lout  in oT_ID.loutl[nd]) + sum(mTEPES.vFlow[sc,p,n,nd,lout ] for lout  in oT_ID.lout[nd]) +
        sum(mTEPES.vLineLosses[sc,p,n,ni,nd,cc] for ni,cc in oT_ID.linl [nd]) - sum(mTEPES.vFlow[sc,p,n,ni,nd,cc] for ni,cc in oT_ID.lin [nd]))
mTEPES.eBalance = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, rule=eBalance, doc='load generation balance [GW]')

print('eBalance              ... ', len(mTEPES.eBalance), ' rows')

def eESSInventory(mTEPES,sc,p,n,es):
    if   mTEPES.n.ord(n) == oT_ID.pESSTimeStep[es]:
        return oT_ID.pESSInitialInventory[es]                                        + sum(oT_ID.pDuration[n2]*1e-3*(oT_ID.pESSEnergyInflows[es][sc,p,n2] - mTEPES.vTotalOutput[sc,p,n2,es] + oT_ID.pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-oT_ID.pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    elif mTEPES.n.ord(n) >  oT_ID.pESSTimeStep[es] and mTEPES.n.ord(n) % oT_ID.pESSTimeStep[es] == 0:
        return mTEPES.vESSInventory[sc,p,mTEPES.n.prev(n,oT_ID.pESSTimeStep[es]),es] + sum(oT_ID.pDuration[n2]*1e-3*(oT_ID.pESSEnergyInflows[es][sc,p,n2] - mTEPES.vTotalOutput[sc,p,n2,es] + oT_ID.pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-oT_ID.pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    else:
        return Constraint.Skip
mTEPES.eESSInventory = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, rule=eESSInventory, doc='ESS inventory balance [TWh]')

print('eESSInventory         ... ', len(mTEPES.eESSInventory), ' rows')

GeneratingRBITime = time.time() - StartTime
StartTime         = time.time()
print('Generating reserves/balance/inventory ... ', round(GeneratingRBITime), 's')

#%%
def eMaxOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   oT_ID.pOperReserveUp[sc,p,n] and oT_ID.pMaxPower2ndBlock[nr][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] + mTEPES.vReserveUp  [sc,p,n,nr]) / oT_ID.pMaxPower2ndBlock[nr][sc,p,n] <= mTEPES.vCommitment[sc,p,n,nr]
    else:
        return Constraint.Skip
mTEPES.eMaxOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMaxOutput2ndBlock, doc='max output of the second block of a committed unit [p.u.]')

print('eMaxOutput2ndBlock    ... ', len(mTEPES.eMaxOutput2ndBlock), ' rows')

def eMinOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   oT_ID.pOperReserveDw[sc,p,n] and oT_ID.pMaxPower2ndBlock[nr][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] - mTEPES.vReserveDown[sc,p,n,nr]) / oT_ID.pMaxPower2ndBlock[nr][sc,p,n] >= 0
    else:
        return Constraint.Skip
mTEPES.eMinOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMinOutput2ndBlock, doc='min output of the second block of a committed unit [p.u.]')

print('eMinOutput2ndBlock    ... ', len(mTEPES.eMinOutput2ndBlock), ' rows')

def eTotalOutput(mTEPES,sc,p,n,nr):
    return mTEPES.vTotalOutput[sc,p,n,nr] == oT_ID.pMinPower[nr][sc,p,n] * mTEPES.vCommitment[sc,p,n,nr] + mTEPES.vOutput2ndBlock[sc,p,n,nr]
mTEPES.eTotalOutput = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput, doc='total output of a unit [GW]')

print('eTotalOutput          ... ', len(mTEPES.eTotalOutput), ' rows')

def eUCStrShut(mTEPES,sc,p,n,nr):
    if n == mTEPES.n.first():
        return mTEPES.vCommitment[sc,p,n,nr] - oT_ID.pInitialUC[nr]                               == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
    else:
        return mTEPES.vCommitment[sc,p,n,nr] - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),nr] == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
mTEPES.eUCStrShut = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eUCStrShut, doc='relation among commitment startup and shutdown')

print('eUCStrShut            ... ', len(mTEPES.eUCStrShut), ' rows')

GeneratingGenConsTime = time.time() - StartTime
StartTime             = time.time()
print('Generating generation constraints     ... ', round(GeneratingGenConsTime), 's')

#%%
def eRampUp(mTEPES,sc,p,n,t):
    if   oT_ID.pRampUp[t] and oT_ID.pRampUp[t] < oT_ID.pMaxPower2ndBlock[t][sc,p,n] and n == mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(oT_ID.pInitialOutput[t]-oT_ID.pMinPower[t][sc,p,n],0)   + mTEPES.vReserveUp  [sc,p,n,t]) / oT_ID.pDuration[n] / oT_ID.pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    elif oT_ID.pRampUp[t] and oT_ID.pRampUp[t] < oT_ID.pMaxPower2ndBlock[t][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t] + mTEPES.vReserveUp  [sc,p,n,t]) / oT_ID.pDuration[n] / oT_ID.pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampUp, doc='maximum ramp up   [p.u.]')

print('eRampUp               ... ', len(mTEPES.eRampUp), ' rows')

def eRampDw(mTEPES,sc,p,n,t):
    if   oT_ID.pRampDw[t] and oT_ID.pRampDw[t] < oT_ID.pMaxPower2ndBlock[t][sc,p,n] and n == mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(oT_ID.pInitialOutput[t]-oT_ID.pMinPower[t][sc,p,n],0)   - mTEPES.vReserveDown[sc,p,n,t]) / oT_ID.pDuration[n] / oT_ID.pRampDw[t] >= - oT_ID.pInitialUC[t]                                                                                           + mTEPES.vShutDown[sc,p,n,t]
    elif oT_ID.pRampDw[t] and oT_ID.pRampDw[t] < oT_ID.pMaxPower2ndBlock[t][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t] - mTEPES.vReserveDown[sc,p,n,t]) / oT_ID.pDuration[n] / oT_ID.pRampDw[t] >= - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),t] + mTEPES.vShutDown[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampDw, doc='maximum ramp down [p.u.]')

print('eRampDw               ... ', len(mTEPES.eRampDw), ' rows')

GeneratingRampsTime = time.time() - StartTime
StartTime           = time.time()
print('Generating ramps   up/down            ... ', round(GeneratingRampsTime), 's')

#%%
def eMinUpTime(mTEPES,sc,p,n,t):
    if oT_ID.pUpTime[t] > 1 and mTEPES.n.ord(n) >= oT_ID.pUpTime[t]:
        return sum(mTEPES.vStartUp [sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-oT_ID.pUpTime[t]:mTEPES.n.ord(n)]) <=     mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinUpTime   = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinUpTime  , doc='minimum up   time [h]')

print('eMinUpTime            ... ', len(mTEPES.eMinUpTime), ' rows')

def eMinDownTime(mTEPES,sc,p,n,t):
    if oT_ID.pDwTime[t] > 1 and mTEPES.n.ord(n) >= oT_ID.pDwTime[t]:
        return sum(mTEPES.vShutDown[sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-oT_ID.pDwTime[t]:mTEPES.n.ord(n)]) <= 1 - mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinDownTime = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinDownTime, doc='minimum down time [h]')

print('eMinDownTime          ... ', len(mTEPES.eMinDownTime), ' rows')

GeneratingMinUDTime = time.time() - StartTime
StartTime           = time.time()
print('Generating minimum up/down time       ... ', round(GeneratingMinUDTime), 's')

#%%
def eInstalNetCap1(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / oT_ID.pLineNTC[ni,nf,cc] >= - mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eInstalNetCap1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap1, doc='maximum flow by installed network capacity [p.u.]')

print('eInstalNetCap1        ... ', len(mTEPES.eInstalNetCap1), ' rows')

def eInstalNetCap2(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / oT_ID.pLineNTC[ni,nf,cc] <=   mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eInstalNetCap2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap2, doc='maximum flow by installed network capacity [p.u.]')

print('eInstalNetCap2        ... ', len(mTEPES.eInstalNetCap2), ' rows')

def eKirchhoff2ndLawExst(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] == (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / oT_ID.pLineX[ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] * oT_ID.pSBase
mTEPES.eKirchhoff2ndLawExst = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lea, rule=eKirchhoff2ndLawExst, doc='flow for each AC existing  line [rad]')

print('eKirchhoff2ndLawExst  ... ', len(mTEPES.eKirchhoff2ndLawExst), ' rows')

def eKirchhoff2ndLawCnd1(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / oT_ID.pLineX[ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] * oT_ID.pSBase >= - 1 + mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eKirchhoff2ndLawCnd1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd1, doc='flow for each AC candidate line [rad]')

print('eKirchhoff2ndLawCnd1  ... ', len(mTEPES.eKirchhoff2ndLawCnd1), ' rows')

def eKirchhoff2ndLawCnd2(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / oT_ID.pLineX[ni,nf,cc] / oT_ID.pMaxFlow[ni,nf,cc] * oT_ID.pSBase <=   1 - mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eKirchhoff2ndLawCnd2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd2, doc='flow for each AC candidate line [rad]')

print('eKirchhoff2ndLawCnd2  ... ', len(mTEPES.eKirchhoff2ndLawCnd2), ' rows')

def eLineLosses1(mTEPES,sc,p,n,ni,nf,cc):
    if oT_ID.pIndNetLosses:
        return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >= - 0.5 * oT_ID.pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
mTEPES.eLineLosses1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses1, doc='ohmic losses for all the lines [GW]')

print('eLineLosses1          ... ', len(mTEPES.eLineLosses1), ' rows')

def eLineLosses2(mTEPES,sc,p,n,ni,nf,cc):
    if oT_ID.pIndNetLosses:
        return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >=   0.5 * oT_ID.pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
mTEPES.eLineLosses2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses2, doc='ohmic losses for all the lines [GW]')

print('eLineLosses2          ... ', len(mTEPES.eLineLosses2), ' rows')

GeneratingNetConsTime          = time.time() - StartTime
print('Generating network    constraints     ... ', round(GeneratingNetConsTime), 's')
