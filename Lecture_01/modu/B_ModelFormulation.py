# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.7.26 - December 28, 2020

import time
from   collections   import defaultdict
from   pyomo.environ import Constraint, Objective, minimize, maximize

def ModelFormulation(model):
    print('Model formulation           ****')

    StartTime = time.time()

    def eTotalGCost(model):
        return model.vTotalGCost == (sum(model.pBidPrice[d]*model.vTotalDem[d] for d in model.d)
                                     - sum(model.pOfferPrice[g]*model.vTotalGen[g] for g in model.g))
    model.eTotalGCost = Constraint(rule=eTotalGCost, doc='social welfare [M$]')

    def eTotalTCost(model):
        return model.vTotalGCost
    model.eTotalTCost = Objective(rule=eTotalTCost, sense=maximize, doc='total system cost [M$]')

    GeneratingOFTime = time.time() - StartTime
    StartTime        = time.time()
    print('Generating objective function         ... ', round(GeneratingOFTime), 's')

    #%% constraints

    # incoming and outgoing lines (lin) (lout)
    lin   = defaultdict(list)
    lout  = defaultdict(list)
    for ni,nf,cc in model.la:
        lin  [nf].append((ni,cc))
        lout [ni].append((nf,cc))

    #%%
    def eBalance(model, nd):
        return (sum(model.vTotalGen[g] for g in model.g if (nd,g) in model.n2g) == sum(model.vTotalDem[d] for d in model.d if (nd,d) in model.n2d) +
                sum(model.vFlow[nd,lout ] for lout  in lout[nd]) -
                sum(model.vFlow[ni,nd,cc] for ni,cc in lin [nd]))
    model.eBalance = Constraint(model.nd, rule=eBalance, doc='load generation balance [GW]')

    print('eBalance              ... ', len(model.eBalance), ' rows')

    def eKirchhoff2ndLawExst(model,ni,nf,cc):
        return model.vFlow[ni,nf,cc] / max(model.pBigMFlowBck[ni,nf,cc](),model.pBigMFlowFrw[ni,nf,cc]()) - (model.vTheta[ni] - model.vTheta[nf]) * model.pLineB[ni,nf,cc] / max(model.pBigMFlowBck[ni,nf,cc](),model.pBigMFlowFrw[ni,nf,cc]()) * model.pSBase == 0
    model.eKirchhoff2ndLawExst = Constraint(model.la, rule=eKirchhoff2ndLawExst, doc='flow for each AC existing  line [rad]')

    print('eKirchhoff2ndLawExst  ... ', len(model.eKirchhoff2ndLawExst), ' rows')


    GeneratingNetConsTime = time.time() - StartTime
    StartTime             = time.time()
    print('Generating network    constraints     ... ', round(GeneratingNetConsTime), 's')
