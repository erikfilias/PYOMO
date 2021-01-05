# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.31 - May 17, 2020

import pandas as pd

# these bounds on vTheta and vFlow do not reduce the problem size after the presolve. They may impact in the crossover but it depends on the case they improve or worsen

# maximum angular difference between any two nodes
import networkx as nx
NetworkGraph = nx.MultiGraph()
NetworkGraph.add_nodes_from(list(mTEPES.nd))
NetworkGraph.add_weighted_edges_from(((ni,nf,pLineX[ni,nf,cc]*pLineNTC[ni,nf,cc]/pSBase) for ni,nf,cc in mTEPES.le))
pMaxThetaDiff = pd.DataFrame.from_dict(dict(nx.all_pairs_dijkstra_path_length(NetworkGraph, cutoff=None, weight='weight')))
pMaxTheta     = pd.DataFrame(pd.concat([pMaxThetaDiff.loc[mTEPES.rf]]*len(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd))).reset_index(drop=True).set_index(pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd)))

# maximum line flow based on the electrical distance
for ni,nf,cc in mTEPES.lca:
    pMaxFlow[ni,nf,cc] = min(pLineNTC[ni,nf,cc], pMaxThetaDiff.loc[ni,nf]*pSBase/pLineX[ni,nf,cc])

# these constraints are too many and they increase the size of the problem and, accordingly, the solution time

#%% maximum angular difference between any two nodes
def eMaxThetaDiff(mTEPES,sc,p,n,ni,nf):
    if ni > nf and nf != mTEPES.pReferenceNode:
        return (-pMaxThetaDiff.loc[ni,nf], vTheta[sc,p,n,ni] - vTheta[sc,p,n,nf], pMaxThetaDiff.loc[ni,nf])
    else:
        return Constraint.Skip
mTEPES.eMaxThetaDiff = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ni, mTEPES.nf, rule=eMaxThetaDiff, doc='maximum angle difference [rad]')

print('eMaxThetaDiff         ... ', len(mTEPES.eMaxThetaDiff), ' rows')
