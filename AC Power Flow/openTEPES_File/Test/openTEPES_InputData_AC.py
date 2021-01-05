# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import math
import pandas as pd
from   pyomo.environ import Param, RangeSet, NonNegativeReals, Var, RealSet

StartTime = time.time()

#%% AC Power Flow: Additional Sets
mTEPES.L                       = RangeSet(20)

#%% AC Power Flow: Additional Parameters
pVmin                = 0.95
pVnom                = 1.00
pVmax                = 1.05
pCapacitivePF        = 0.95                #Capacitive Power Factor
pInductivePF         = 0.99                #Inductive Power Factor

# pBusGshb           = dfNodeLocation.drop(['Latitude','Longitude'], axis=1).assign(gshb=0)
# pBusBshb           = dfNodeLocation.drop(['Latitude','Longitude'], axis=1).assign(bshb=0)

pBusGshb             = pd.Series([0.2]*len(mTEPES.nd), index=pd.Index(mTEPES.nd))
pBusBshb             = pd.Series([0.2]*len(mTEPES.nd), index=pd.Index(mTEPES.nd))

#%% Parameters: Branches
mTEPES.pBusGshb      = Param(mTEPES.nd, initialize=pBusGshb.to_dict(), within=NonNegativeReals, doc='Conductance'             )
mTEPES.pBusBshb      = Param(mTEPES.nd, initialize=pBusBshb.to_dict(), within=NonNegativeReals, doc='Susceptance'                              )
mTEPES.pLineFi       = Param(mTEPES.la, initialize=(0.*math.pi)/180,                            doc='xxxx'                                     )
mTEPES.pLineSmax     = Param(mTEPES.la, initialize=0.,                                          doc='xxxx',                        mutable=True)
mTEPES.pLineZ        = Param(mTEPES.la, initialize=0.,                                          doc='xxxx',                        mutable=True)
mTEPES.pLineZ2       = Param(mTEPES.la, initialize=0.,                                          doc='xxxx',                        mutable=True)
mTEPES.pLineDelta_S  = Param(mTEPES.la, initialize=0.,                 within=NonNegativeReals, doc='Delta of Smax splitted by L', mutable=True)
mTEPES.pLineM        = Param(mTEPES.la,mTEPES.L, initialize = 0.,      within=NonNegativeReals, doc='M partitions of Delta Smax' , mutable=True)

#%% Parameters: Power System
mTEPES.pVmin         = Param(           initialize=pVmin               , within=NonNegativeReals, doc='Minimum voltage magnitude')
mTEPES.pVnom         = Param(           initialize=pVnom               , within=NonNegativeReals, doc='Nominal voltage magnitude')
mTEPES.pVmax         = Param(           initialize=pVmax               , within=NonNegativeReals, doc='Maximum voltage magnitude')
mTEPES.pCapacitivePF = Param(           initialize=pCapacitivePF       , within=NonNegativeReals, doc='Capacitive power factor'  )
mTEPES.pInductivePF  = Param(           initialize=pInductivePF        , within=NonNegativeReals, doc='Inductive power factor'   )

for la in mTEPES.la:
    # he comentado estas dos líneas porque no están definidos como mutables estos parámetros, hay que modificarlos en los datos
    # mTEPES.pLineBsh    [la] = mTEPES.pLineBsh[la] / 2
    # mTEPES.pLineTAP    [la] = 1. / mTEPES.pLineTAP[la]
    mTEPES.pLineSmax   [la] = mTEPES.pLineNTC[la]*1.5
    mTEPES.pLineZ2     [la] = mTEPES.pLineR  [la]**2 + mTEPES.pLineX[la]**2
    mTEPES.pLineDelta_S[la] = mTEPES.pLineSmax[la] / len(mTEPES.L)
    for l in mTEPES.L:
        mTEPES.pLineM[la,l] = (2*l-1)*mTEPES.pLineDelta_S[la]

#%% AC Power Flow: Additional Variables
# def _bounds_Delta_S_rule(mTEPES, sc, p, n, la, L):
#     return (0, mTEPES.pLineDelta_S[la])
mTEPES.vDelta_P             = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,ni,nf,cc,l:(0, mTEPES.pLineDelta_S[ni,nf,cc]), doc='Delta Active Power Flow                                       [  GW]')
mTEPES.vDelta_Q             = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,ni,nf,cc,l:(0, mTEPES.pLineDelta_S[ni,nf,cc]), doc='Delta Reactive Power Flow                                     [GVaR]')
mTEPES.vP_max               = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=NonNegativeReals,                                                                            doc='Maximum bound of Active Power Flow                            [  GW]')
mTEPES.vP_min               = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=NonNegativeReals,                                                                            doc='Minimum bound of Active Power Flow                            [  GW]')
mTEPES.vQ_max               = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=NonNegativeReals,                                                                            doc='Maximum bound of Reactive Power Flow                          [GVaR]')
mTEPES.vQ_min               = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=NonNegativeReals,                                                                            doc='Minimum bound of Reactive Power Flow                          [GVaR]')

# #VARIABLES
mTEPES.vVoltageMag_sqr      = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd,           within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nd:(pVmin**2,pVmax**2),                                        doc='Voltage magnitude                                   [p.u.]')
mTEPES.vCurrentFlow_sqr     = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,*la:(0,(mTEPES.pLineNTC[la]**2/pVmax)),                               doc='Current flow                                           [A]')
mTEPES.vP                   = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=RealSet,                                                                                                     doc='Active Power Flow through lines                                 [GW]')
mTEPES.vQ                   = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la,           within=RealSet,                                                                                                     doc='Reactive Power Flow through lines                               [GW]')
mTEPES.vReactiveTotalOutput = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g ,           within=RealSet,          bounds=lambda mTEPES,sc,p,n,g :(-mTEPES.pMaxPower[sc,p,n,g ],mTEPES.pMaxPower[sc,p,n,g ]), doc='Total output of reactive power generators                     [GVAr]')

SettingUpDataTime = time.time() - StartTime
StartTime         = time.time()
print('Setting up input data                 ... ', round(SettingUpDataTime), 's')

# import numpy as np
# pLineR = mTEPES.pLineR
# pLineX = mTEPES.pLineX
# pLineZ = pLineR + pLineX*1j
# pLineY = 1/pLineZ
# Yb     = np.zeros((len(mTEPES.nd), len(mTEPES.nd)), dtype=complex)
#
# dfPosition = pd.DataFrame(columns=['Position'], index = mTEPES.nd)
# Count = 0
# for i in mTEPES.nd:
#     dfPosition['Position'][i] = Count
#     Count += 1
#
# mTEPES.dfPosition = dfPosition
# # print(Yb)
# for ni,nf,cc in mTEPES.la:
#     Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]] = Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]] + mTEPES.pLineY[ni,nf,cc]*mTEPES.pLineTAP[ni,nf,cc]
#     Yb[dfPosition['Position'][nf], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]]
#
# for k in mTEPES.nd:
#     for ni,nf,cc in mTEPES.la:
#         if ni == k:
#             Yb[dfPosition['Position'][ni], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni],dfPosition['Position'][ni]] + mTEPES.pLineY[ni,nf,cc]*mTEPES.pLineTAP[ni,nf,cc]**2 + mTEPES.pLineBsh[ni,nf,cc]*1j
#         elif nf == k:
#             Yb[dfPosition['Position'][ni], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni],dfPosition['Position'][ni]] + mTEPES.pLineY[ni,nf,cc]                              + mTEPES.pLineBsh[ni,nf,cc]*1j
#
# for k in mTEPES.nd:
#     Yb[dfPosition['Position'][k], dfPosition['Position'][k]] = Yb[dfPosition['Position'][k], dfPosition['Position'][k]] + mTEPES.pBusBshb[k] * 1j
#
# mTEPES.Yb = Yb
# # print(Yb)
#
# import numpy as np
# Ybarra                      = np.abs(mTEPES.Yb)/np.abs(mTEPES.Yb)
# Ybarra[np.isnan(Ybarra)]    = 0
#
# from cvxopt import spmatrix, amd
# from scipy.sparse import csr_matrix, find
# Ybarra                      = csr_matrix(Ybarra.real)
# #find(Ybarra)
# coo = Ybarra.tocoo()
# SP = spmatrix(coo.data, coo.row.tolist(), coo.col.tolist())
# isinstance(SP,spmatrix)
#
# from cvxopt import spmatrix, amd
# reordening = amd.order(SP)
# # print(reordening)
# Yorden                      = np.zeros((len(mTEPES.nd), len(mTEPES.nd)), dtype=complex)
# for ni in mTEPES.nd:
#     for nf in mTEPES.nd:
#         Yorden[mTEPES.dfPosition['Position'][ni], mTEPES.dfPosition['Position'][nf]] = mTEPES.Yb[reordening[mTEPES.dfPosition['Position'][ni]],reordening[mTEPES.dfPosition['Position'][nf]]]
#
# # print(Yorden)
#
# mTEPES.Yorden = Yorden
#
# # from cvxopt import spmatrix, amd, normal
# # from chompack import symbolic, cspmatrix, cholesky
# # pattern=cholesky(abs(Yorden.imag))