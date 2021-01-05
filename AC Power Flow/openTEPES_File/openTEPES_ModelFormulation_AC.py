# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.32 - May 26, 2020

import time
from   collections   import defaultdict
from   pyomo.environ import Constraint, Objective, minimize, tan, acos

print('Model formulation     ****')

StartTime = time.time()

def eTotalFCost(mTEPES):
   return mTEPES.vTotalFCost == sum(mTEPES.pNetFixedCost[lc] * mTEPES.vNetworkInvest   [lc] for lc in mTEPES.lc) + sum(mTEPES.pGenFixedCost[gc] * mTEPES.vGenerationInvest[gc] for gc in mTEPES.gc)
mTEPES.eTotalFCost = Constraint(rule=eTotalFCost, doc='total system fixed    cost [MEUR]')

def eTotalVCost(mTEPES):
    return mTEPES.vTotalVCost == (sum(mTEPES.pScenProb[sc] * mTEPES.pDuration[n] * mTEPES.pENSCost             * mTEPES.vENS        [sc,p,n,nd] * mTEPES.pDemand[sc,p,n,nd] * (1+tan(acos(mTEPES.pCapacitivePF))) for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd) +
                                  sum(mTEPES.pScenProb[sc] * mTEPES.pDuration[n] * mTEPES.pLinearVarCost  [nr] * mTEPES.vTotalOutput[sc,p,n,nr]                                                         +
                                      mTEPES.pScenProb[sc] * mTEPES.pDuration[n] * mTEPES.pConstantVarCost[nr] * mTEPES.vCommitment [sc,p,n,nr]                                                         +
                                      mTEPES.pScenProb[sc]                       * mTEPES.pStartUpCost    [nr] * mTEPES.vStartUp    [sc,p,n,nr]                                                         +
                                      mTEPES.pScenProb[sc]                       * mTEPES.pShutDownCost   [nr] * mTEPES.vShutDown   [sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr) )
mTEPES.eTotalVCost = Constraint(rule=eTotalVCost, doc='total system variable cost [MEUR]')

def eTotalECost(mTEPES):
    return mTEPES.vTotalECost == sum(mTEPES.pScenProb[sc] * mTEPES.pCO2Cost * mTEPES.pCO2EmissionRate[nr] * mTEPES.vTotalOutput[sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)
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
    return mTEPES.vTotalOutput[sc,p,n,ec] / mTEPES.pMaxPower[sc,p,n,ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalGenESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalGenESS, doc='output      if installed ESS unit [p.u.]')

print('eInstalGenESS         ... ', len(mTEPES.eInstalGenESS), ' rows')

def eInstalConESS(mTEPES,sc,p,n,ec):
    return mTEPES.vTotalOutput[sc,p,n,ec] / mTEPES.pMaxCharge      [ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalConESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalConESS, doc='consumption if installed ESS unit [p.u.]')

print('eInstalConESS         ... ', len(mTEPES.eInstalConESS), ' rows')

#%%
def eOperReserveUp(mTEPES,sc,p,n):
    if mTEPES.pOperReserveUp[sc,p,n]:
        return sum(mTEPES.vReserveUp  [sc,p,n,nr] for nr in mTEPES.nr) >= mTEPES.pOperReserveUp[sc,p,n]
mTEPES.eOperReserveUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveUp, doc='up   operating reserve [GW]')

print('eOperReserveUp        ... ', len(mTEPES.eOperReserveUp), ' rows')

def eOperReserveDw(mTEPES,sc,p,n):
    if mTEPES.pOperReserveDw[sc,p,n]:
        return sum(mTEPES.vReserveDown[sc,p,n,nr] for nr in mTEPES.nr) >= mTEPES.pOperReserveDw[sc,p,n]
mTEPES.eOperReserveDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveDw, doc='down operating reserve [GW]')

print('eOperReserveDw        ... ', len(mTEPES.eOperReserveDw), ' rows')

# incoming and outgoing lines (lin) (lout)
lin   = defaultdict(list)
linl  = defaultdict(list)
lout  = defaultdict(list)
loutl = defaultdict(list)
for ni,nf,cc in mTEPES.la:
    lin  [nf].append((ni,cc))
for ni,nf,cc in mTEPES.la:
    lout [ni].append((nf,cc))
for ni,nf,cc in mTEPES.ll:
    linl [nf].append((ni,cc))
for ni,nf,cc in mTEPES.ll:
    loutl[ni].append((nf,cc))

def eBalance_P(mTEPES,sc,p,n,nd):
    return (sum(mTEPES.vTotalOutput[sc,p,n,g] for g in mTEPES.g if (nd,g) in mTEPES.n2g) - sum(mTEPES.vESSCharge[sc,p,n,es] for es in mTEPES.es if (nd,es) in mTEPES.n2g) == mTEPES.pDemand[sc,p,n,nd] * (1 - mTEPES.vENS[sc,p,n,nd]) +
        sum(mTEPES.vP[sc,p,n,nd,lout ] + mTEPES.pLineR[nd,lout ] * mTEPES.vCurrentFlow_sqr[sc,p,n,nd,lout ] for lout  in lout[nd]) -
        sum(mTEPES.vP[sc,p,n,ni,nd,cc] for ni,cc in lin [nd]) + mTEPES.vVoltageMag_sqr[sc,p,n,nd] * mTEPES.pBusGshb[nd])
mTEPES.eBalance_P = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, rule=eBalance_P, doc='Active power generation balance [GW]')
        # sum(mTEPES.vFlow      [sc,p,n,nd,nf,cc] for nf,cc in mTEPES.nf*mTEPES.cc if (nd,nf,cc) in mTEPES.la ) -
        # sum(mTEPES.vFlow      [sc,p,n,ni,nd,cc] for ni,cc in mTEPES.ni*mTEPES.cc if (nd,ni,cc) in mTEPES.lin) +
        # sum(mTEPES.vLineLosses[sc,p,n,nd,nf,cc] for nf,cc in mTEPES.nf*mTEPES.cc if (nd,nf,cc) in mTEPES.ll ) +
        # sum(mTEPES.vLineLosses[sc,p,n,ni,nd,cc] for ni,cc in mTEPES.ni*mTEPES.cc if (nd,ni,cc) in mTEPES.lil) )
        # sum(mTEPES.vFlow      [sc,p,n,nd,nf,cc] for nf,cc in mTEPES.nf*mTEPES.cc if (nd,nf,cc) in mTEPES.la) -
        # sum(mTEPES.vFlow      [sc,p,n,ni,nd,cc] for ni,cc in mTEPES.ni*mTEPES.cc if (ni,nd,cc) in mTEPES.la) +
        # sum(mTEPES.vLineLosses[sc,p,n,nd,nf,cc] for nf,cc in mTEPES.nf*mTEPES.cc if (nd,nf,cc) in mTEPES.ll) +
        # sum(mTEPES.vLineLosses[sc,p,n,ni,nd,cc] for ni,cc in mTEPES.ni*mTEPES.cc if (ni,nd,cc) in mTEPES.ll) )
print('eBalance_P            ... ', len(mTEPES.eBalance_P), ' rows')

def eBalance_Q(mTEPES,sc,p,n,nd):
    return (sum(mTEPES.vReactiveTotalOutput[sc,p,n,g] for g in mTEPES.g if (nd,g) in mTEPES.n2g) == mTEPES.pDemand[sc,p,n,nd]* (tan(acos(mTEPES.pCapacitivePF))) * (1 - mTEPES.vENS[sc,p,n,nd]) +
        sum(mTEPES.vQ[sc,p,n,nd,lout ] - mTEPES.pLineBsh[nd,lout ] * mTEPES.vVoltageMag_sqr[sc,p,n,nd] + mTEPES.pLineX[nd,lout] * mTEPES.vCurrentFlow_sqr[sc,p,n,nd,lout ] for lout  in lout[nd]) -
        sum(mTEPES.vQ[sc,p,n,ni,nd,cc] + mTEPES.pLineBsh[ni,nd,cc] * mTEPES.vVoltageMag_sqr[sc,p,n,nd] for ni,cc in lin [nd]) - mTEPES.vVoltageMag_sqr[sc,p,n,nd] * mTEPES.pBusBshb[nd])
mTEPES.eBalance_Q = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, rule=eBalance_Q, doc='Reactive power generation balance [GVaR]')
print('eBalance_Q            ... ', len(mTEPES.eBalance_Q), ' rows')

def eESSInventory(mTEPES,sc,p,n,es):
    if   mTEPES.n.ord(n) == mTEPES.pESSTimeStep[es]:
        return mTEPES.pESSInitialInventory[es]                                        + sum(mTEPES.pDuration[n2]*1e-3*(mTEPES.pESSEnergyInflows[sc,p,n2,es] - mTEPES.vTotalOutput[sc,p,n2,es] + mTEPES.pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-mTEPES.pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    elif mTEPES.n.ord(n) >  mTEPES.pESSTimeStep[es] and mTEPES.n.ord(n) % mTEPES.pESSTimeStep[es] == 0:
        return mTEPES.vESSInventory[sc,p,mTEPES.n.prev(n,mTEPES.pESSTimeStep[es]),es] + sum(mTEPES.pDuration[n2]*1e-3*(mTEPES.pESSEnergyInflows[sc,p,n2,es] - mTEPES.vTotalOutput[sc,p,n2,es] + mTEPES.pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-mTEPES.pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    else:
        return Constraint.Skip
mTEPES.eESSInventory = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, rule=eESSInventory, doc='ESS inventory balance [TWh]')

print('eESSInventory         ... ', len(mTEPES.eESSInventory), ' rows')

GeneratingRBITime = time.time() - StartTime
StartTime         = time.time()
print('Generating reserves/balance/inventory ... ', round(GeneratingRBITime), 's')

#%%
def eMaxOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   mTEPES.pOperReserveUp[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] + mTEPES.vReserveUp  [sc,p,n,nr]) / mTEPES.pMaxPower2ndBlock[sc,p,n,nr] <= mTEPES.vCommitment[sc,p,n,nr]
    else:
        return Constraint.Skip
mTEPES.eMaxOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMaxOutput2ndBlock, doc='max output of the second block of a committed unit [p.u.]')

print('eMaxOutput2ndBlock    ... ', len(mTEPES.eMaxOutput2ndBlock), ' rows')

def eMinOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   mTEPES.pOperReserveDw[sc,p,n] and mTEPES.pMaxPower2ndBlock[sc,p,n,nr]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] - mTEPES.vReserveDown[sc,p,n,nr]) / mTEPES.pMaxPower2ndBlock[sc,p,n,nr] >= 0
    else:
        return Constraint.Skip
mTEPES.eMinOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMinOutput2ndBlock, doc='min output of the second block of a committed unit [p.u.]')

print('eMinOutput2ndBlock    ... ', len(mTEPES.eMinOutput2ndBlock), ' rows')

def eTotalOutput_P(mTEPES,sc,p,n,nr):
    return mTEPES.vTotalOutput[sc,p,n,nr] == mTEPES.pMinPower[sc,p,n,nr] * mTEPES.vCommitment[sc,p,n,nr] + mTEPES.vOutput2ndBlock[sc,p,n,nr]
mTEPES.eTotalOutput_P = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput_P, doc='total output of active power by a unit [GW]')

print('eTotalOutput_P        ... ', len(mTEPES.eTotalOutput_P), ' rows')

def eTotalOutput_Q_Capacitive(mTEPES,sc,p,n,nr):
    return mTEPES.vReactiveTotalOutput[sc,p,n,nr] <=  mTEPES.vTotalOutput[sc,p,n,nr] * (tan(acos(mTEPES.pCapacitivePF)))
mTEPES.eTotalOutput_Q1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput_Q_Capacitive, doc='total output of reactive power by a unit in the capacitive side [GW]')

print('eTotalOutput_Q1       ... ', len(mTEPES.eTotalOutput_Q1), ' rows')

def eTotalOutput_Q_Inductive(mTEPES,sc,p,n,nr):
    return mTEPES.vReactiveTotalOutput[sc,p,n,nr] >= -mTEPES.vTotalOutput[sc,p,n,nr] * (tan(acos(mTEPES.pInductivePF)))
mTEPES.eTotalOutput_Q2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput_Q_Inductive, doc='total output of reactive power by a unit in the inductive side [GVaR]')

print('eTotalOutput_Q2       ... ', len(mTEPES.eTotalOutput_Q2), ' rows')

def eUCStrShut(mTEPES,sc,p,n,nr):
    if n == mTEPES.n.first():
        return mTEPES.vCommitment[sc,p,n,nr] - mTEPES.pInitialUC[nr]                        == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
    else:
        return mTEPES.vCommitment[sc,p,n,nr] - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),nr] == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
mTEPES.eUCStrShut = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eUCStrShut, doc='relation among commitment startup and shutdown')

print('eUCStrShut            ... ', len(mTEPES.eUCStrShut), ' rows')

GeneratingGenConsTime = time.time() - StartTime
StartTime             = time.time()
print('Generating generation constraints     ... ', round(GeneratingGenConsTime), 's')

#%%
def eRampUp(mTEPES,sc,p,n,t):
    if   mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and n == mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(mTEPES.pInitialOutput[t]-mTEPES.pMinPower[sc,p,n,t],0) + mTEPES.vReserveUp  [sc,p,n,t]) / mTEPES.pDuration[n] / mTEPES.pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    elif mTEPES.pRampUp[t] and mTEPES.pRampUp[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t]             + mTEPES.vReserveUp  [sc,p,n,t]) / mTEPES.pDuration[n] / mTEPES.pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampUp, doc='maximum ramp up   [p.u.]')

print('eRampUp               ... ', len(mTEPES.eRampUp), ' rows')

def eRampDw(mTEPES,sc,p,n,t):
    if   mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t] and n == mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(mTEPES.pInitialOutput[t]-mTEPES.pMinPower[sc,p,n,t],0) - mTEPES.vReserveDown[sc,p,n,t]) / mTEPES.pDuration[n] / mTEPES.pRampDw[t] >= - mTEPES.pInitialUC[t]                                                                                           + mTEPES.vShutDown[sc,p,n,t]
    elif mTEPES.pRampDw[t] and mTEPES.pRampDw[t] < mTEPES.pMaxPower2ndBlock[sc,p,n,t]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t]             - mTEPES.vReserveDown[sc,p,n,t]) / mTEPES.pDuration[n] / mTEPES.pRampDw[t] >= - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),t] + mTEPES.vShutDown[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampDw, doc='maximum ramp down [p.u.]')

print('eRampDw               ... ', len(mTEPES.eRampDw), ' rows')

GeneratingRampsTime = time.time() - StartTime
StartTime           = time.time()
print('Generating ramps   up/down            ... ', round(GeneratingRampsTime), 's')

#%%
def eMinUpTime(mTEPES,sc,p,n,t):
    if mTEPES.pUpTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pUpTime[t]:
        return sum(mTEPES.vStartUp [sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-mTEPES.pUpTime[t]:mTEPES.n.ord(n)]) <=     mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinUpTime   = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinUpTime  , doc='minimum up   time [h]')

print('eMinUpTime            ... ', len(mTEPES.eMinUpTime), ' rows')

def eMinDownTime(mTEPES,sc,p,n,t):
    if mTEPES.pDwTime[t] > 1 and mTEPES.n.ord(n) >= mTEPES.pDwTime[t]:
        return sum(mTEPES.vShutDown[sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-mTEPES.pDwTime[t]:mTEPES.n.ord(n)]) <= 1 - mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinDownTime = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinDownTime, doc='minimum down time [h]')

print('eMinDownTime          ... ', len(mTEPES.eMinDownTime), ' rows')

GeneratingMinUDTime = time.time() - StartTime
StartTime           = time.time()
print('Generating minimum up/down time       ... ', round(GeneratingMinUDTime), 's')

#%%
def eInstalNetCap(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc] / (mTEPES.pLineNTC[ni,nf,cc] ** 2 / mTEPES.pVmax) <= mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eInstalNetCap = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap, doc='maximum flow by installed network capacity [p.u.]')

print('eInstalNetCap         ... ', len(mTEPES.eInstalNetCap), ' rows')

def eVoltageDiffExst(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vVoltageMag_sqr[sc,p,n,ni] * (mTEPES.pLineTAP[ni,nf,cc]**2) - mTEPES.vVoltageMag_sqr[sc,p,n,nf] == 
            mTEPES.pLineZ2[ni,nf,cc] * mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc] + 
            2 * (mTEPES.pLineR[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] + mTEPES.pLineX[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc]))
mTEPES.eVoltageDiffExst = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lea, rule=eVoltageDiffExst, doc='Voltage difference between nodes considering existent lines [p.u.]')

print('eVoltageDiffExst      ... ', len(mTEPES.eVoltageDiffExst), ' rows')

def eVoltageDiffCand1(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vVoltageMag_sqr[sc,p,n,ni] * (mTEPES.pLineTAP[ni,nf,cc]**2) - mTEPES.vVoltageMag_sqr[sc,p,n,nf] -
            mTEPES.pLineZ2[ni,nf,cc] * mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc] - 
            2 * (mTEPES.pLineR[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] + mTEPES.pLineX[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc])) <=  (mTEPES.pVmax**2 - mTEPES.pVmin**2) * mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eVoltageDiffCand1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eVoltageDiffCand1, doc='Voltage difference between nodes considering existent lines [p.u.]')

print('eVoltageDiffCand1     ... ', len(mTEPES.eVoltageDiffCand1), ' rows')

def eVoltageDiffCand2(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vVoltageMag_sqr[sc,p,n,ni] * (mTEPES.pLineTAP[ni,nf,cc]**2) - mTEPES.vVoltageMag_sqr[sc,p,n,nf] -
            mTEPES.pLineZ2[ni,nf,cc] * mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc] - 
            2 * (mTEPES.pLineR[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] + mTEPES.pLineX[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc])) >= -(mTEPES.pVmax**2 - mTEPES.pVmin**2) * mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eVoltageDiffCand2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eVoltageDiffCand2, doc='Voltage magnitude difference between nodes considering existent lines [p.u.]')

print('eVoltageDiffCand2     ... ', len(mTEPES.eVoltageDiffCand2), ' rows')

def eAngleDiffExst(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.pVnom**2 * (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) == 
            (mTEPES.pLineX[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] - mTEPES.pLineR[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc]))
mTEPES.eAngleDiffExst = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lea, rule=eAngleDiffExst, doc='Voltage angle difference between nodes considering existent lines [p.u.]')

print('eAngleDiffExst        ... ', len(mTEPES.eAngleDiffExst), ' rows')

def eAngleDiffCand1(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.pVnom**2 * (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) - 
            (mTEPES.pLineX[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] - mTEPES.pLineR[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc])) <=  2 * mTEPES.pMaxTheta[sc,p,n,ni] * mTEPES.pVmax**2 * mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eAngleDiffCand1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eAngleDiffCand1, doc='Voltage angle difference between nodes considering candidate lines [p.u.]')

print('eAngleDiffCand1       ... ', len(mTEPES.eAngleDiffCand1), ' rows')

def eAngleDiffCand2(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.pVnom**2 * (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) - 
            (mTEPES.pLineX[ni,nf,cc] * mTEPES.vP[sc,p,n,ni,nf,cc] - mTEPES.pLineR[ni,nf,cc] * mTEPES.vQ[sc,p,n,ni,nf,cc])) >= -2 * mTEPES.pMaxTheta[sc,p,n,ni] * mTEPES.pVmax**2 * mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eAngleDiffCand2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eAngleDiffCand2, doc='Voltage angle difference between nodes considering candidate lines [p.u.]')

print('eAngleDiffCand2       ... ', len(mTEPES.eAngleDiffCand2), ' rows')

def eCurrentFlow(mTEPES,sc,p,n,ni,nf,cc):
    return ((mTEPES.pVnom**2) * mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc] == sum(mTEPES.pLineM[ni,nf,cc,l] * mTEPES.vDelta_P[sc,p,n,ni,nf,cc,l] for l in mTEPES.L) + sum(mTEPES.pLineM[ni,nf,cc,l] * mTEPES.vDelta_Q[sc,p,n,ni,nf,cc,l] for l in mTEPES.L))
mTEPES.eCurrentFlow = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, rule=eCurrentFlow, doc='Linear constraint for current flow  [p.u.]')

print('eCurrentFlow          ... ', len(mTEPES.eCurrentFlow), ' rows')

def eCurrentFlow_LinearP1(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vP_max[sc,p,n,ni,nf,cc] - mTEPES.vP_min[sc,p,n,ni,nf,cc] == mTEPES.vP[sc,p,n,ni,nf,cc])
mTEPES.eCurrentFlow_LinearP1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, rule=eCurrentFlow_LinearP1, doc='Linear constraint for current flow relate to P [p.u.]')

print('eCurrentFlow_LinearP1 ... ', len(mTEPES.eCurrentFlow_LinearP1), ' rows')

def eCurrentFlow_LinearP2(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vP_max[sc,p,n,ni,nf,cc] + mTEPES.vP_min[sc,p,n,ni,nf,cc] == sum(mTEPES.vDelta_P[sc,p,n,ni,nf,cc,l] for l in mTEPES.L))
mTEPES.eCurrentFlow_LinearP2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, rule=eCurrentFlow_LinearP2, doc='Linear constraint for current flow relate to P [p.u.]')

print('eCurrentFlow_LinearP2 ... ', len(mTEPES.eCurrentFlow_LinearP2), ' rows')

def eCurrentFlow_LinearQ1(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vQ_max[sc,p,n,ni,nf,cc] - mTEPES.vQ_min[sc,p,n,ni,nf,cc] == mTEPES.vQ[sc,p,n,ni,nf,cc])
mTEPES.eCurrentFlow_LinearQ1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, rule=eCurrentFlow_LinearQ1, doc='Linear constraint for current flow relate to Q [p.u.]')

print('eCurrentFlow_LinearQ1 ... ', len(mTEPES.eCurrentFlow_LinearQ1), ' rows')

def eCurrentFlow_LinearQ2(mTEPES,sc,p,n,ni,nf,cc):
    return (mTEPES.vQ_max[sc,p,n,ni,nf,cc] + mTEPES.vQ_min[sc,p,n,ni,nf,cc] == sum(mTEPES.vDelta_Q[sc,p,n,ni,nf,cc,l] for l in mTEPES.L))
mTEPES.eCurrentFlow_LinearQ2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, rule=eCurrentFlow_LinearQ2, doc='Linear constraint for current flow relate to Q [p.u.]')

print('eCurrentFlow_LinearQ2 ... ', len(mTEPES.eCurrentFlow_LinearQ2), ' rows')

# def eExisteNetCap(mTEPES,sc,p,n,ni,nf,cc):
#     return mTEPES.vFlow[sc,p,n,ni,nf,cc] / mTEPES.pLineNTC[ni,nf,cc] <=   mTEPES.vNetworkInvest[ni,nf,cc]
# mTEPES.eInstalNetCap2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap2, doc='maximum flow by installed network capacity [p.u.]')

# print('eInstalNetCap2        ... ', len(mTEPES.eInstalNetCap2), ' rows')

# def eKirchhoff2ndLawExst(mTEPES,sc,p,n,ni,nf,cc):
#     return mTEPES.vFlow[sc,p,n,ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] == (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / mTEPES.pLineX[ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] * mTEPES.pSBase
# mTEPES.eKirchhoff2ndLawExst = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lea, rule=eKirchhoff2ndLawExst, doc='flow for each AC existing  line [rad]')

# print('eKirchhoff2ndLawExst  ... ', len(mTEPES.eKirchhoff2ndLawExst), ' rows')

# def eKirchhoff2ndLawCnd1(mTEPES,sc,p,n,ni,nf,cc):
#     return mTEPES.vFlow[sc,p,n,ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / mTEPES.pLineX[ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] * mTEPES.pSBase >= - 1 + mTEPES.vNetworkInvest[ni,nf,cc]
# mTEPES.eKirchhoff2ndLawCnd1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd1, doc='flow for each AC candidate line [rad]')

# print('eKirchhoff2ndLawCnd1  ... ', len(mTEPES.eKirchhoff2ndLawCnd1), ' rows')

# def eKirchhoff2ndLawCnd2(mTEPES,sc,p,n,ni,nf,cc):
#     return mTEPES.vFlow[sc,p,n,ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / mTEPES.pLineX[ni,nf,cc] / mTEPES.pMaxFlow[ni,nf,cc] * mTEPES.pSBase <=   1 - mTEPES.vNetworkInvest[ni,nf,cc]
# mTEPES.eKirchhoff2ndLawCnd2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd2, doc='flow for each AC candidate line [rad]')

# print('eKirchhoff2ndLawCnd2  ... ', len(mTEPES.eKirchhoff2ndLawCnd2), ' rows')

# def eLineLosses1(mTEPES,sc,p,n,ni,nf,cc):
#     if mTEPES.pIndNetLosses:
#         return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >= - 0.5 * mTEPES.pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
# mTEPES.eLineLosses1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses1, doc='ohmic losses for all the lines [GW]')

# print('eLineLosses1          ... ', len(mTEPES.eLineLosses1), ' rows')

# def eLineLosses2(mTEPES,sc,p,n,ni,nf,cc):
#     if mTEPES.pIndNetLosses:
#         return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >=   0.5 * mTEPES.pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
# mTEPES.eLineLosses2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses2, doc='ohmic losses for all the lines [GW]')

# print('eLineLosses2          ... ', len(mTEPES.eLineLosses2), ' rows')

GeneratingNetConsTime = time.time() - StartTime
StartTime             = time.time()
print('Generating network    constraints     ... ', round(GeneratingNetConsTime), 's')

# mTEPES.write('openTEPES_'+CaseName+'.lp', io_options={'symbolic_solver_labels': True})  # create lp-format file

WritingLPFileTime = time.time() - StartTime
StartTime         = time.time()
print('Writing LP file                       ... ', round(WritingLPFileTime), 's')
