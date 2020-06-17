# Open Generation and Transmission Operation and Expansion Planning Model with RES and  (openTEPES) - Version 1.6.33 - June 17, 2020

from   pyomo.environ import Constraint, tan, acos

def eTotalOutput_Q_Capacitive(mTEPES,sc,p,n,nr):
    return mTEPES.vReactiveTotalOutput[sc,p,n,nr] <=   mTEPES.vTotalOutput[sc,p,n,nr] * (tan(acos(mTEPES.pCapacitivePF)))
mTEPES.eTotalOutput_Q1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput_Q_Capacitive, doc='total output of reactive power by a unit in the capacitive side [GW]')

print('eTotalOutput_Q1       ... ', len(mTEPES.eTotalOutput_Q1), ' rows')

def eTotalOutput_Q_Inductive(mTEPES,sc,p,n,nr):
    return mTEPES.vReactiveTotalOutput[sc,p,n,nr] >= - mTEPES.vTotalOutput[sc,p,n,nr] * (tan(acos(mTEPES.pInductivePF)))
mTEPES.eTotalOutput_Q2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput_Q_Inductive, doc='total output of reactive power by a unit in the inductive side [GVaR]')

print('eTotalOutput_Q2       ... ', len(mTEPES.eTotalOutput_Q2), ' rows')

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
