# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.7.26 - December 28, 2020

import time
import math
import pandas        as pd
from   pyomo.environ import DataPortal, Set, Param, Var, Binary, NonNegativeReals, Reals, UnitInterval, Boolean, Any
import os


def InputData(CaseName,model):
    print('Input data                  ****')

    # get current directory
    model.path = os.getcwd()
    print("Current Directory", model.path)

    StartTime = time.time()
    #%% reading data from CSV
    dfParameter          = pd.read_csv(model.path + '/par/m_data_parameter_'            +CaseName+'.csv', index_col=[0    ])
    dfDemand             = pd.read_csv(model.path + '/par/m_data_demand_'               +CaseName+'.csv', index_col=[0    ])
    dfGeneration         = pd.read_csv(model.path + '/par/m_data_generation_'           +CaseName+'.csv', index_col=[0    ])
    dfNetwork            = pd.read_csv(model.path + '/par/m_data_network_'              +CaseName+'.csv', index_col=[0,1,2])

    # substitute NaN by 0
    dfDemand.fillna            (0.0, inplace=True)
    dfGeneration.fillna        (0.0, inplace=True)
    dfNetwork.fillna           (0.0, inplace=True)

    # show some statistics of the data
    print('Demand                       \n', dfDemand.describe()            )
    print('Generation                   \n', dfGeneration.describe()        )
    print('Network                      \n', dfNetwork.describe()           )

    #%% reading the sets
    dictSets = DataPortal()
    dictSets.load(filename=model.path+'/set/m_dict_generation_'  +CaseName+'.csv', set='g'   , format='set')
    dictSets.load(filename=model.path+'/set/m_dict_node_'        +CaseName+'.csv', set='nd'  , format='set')
    dictSets.load(filename=model.path+'/set/m_dict_circuit_'     +CaseName+'.csv', set='cc'  , format='set')
    dictSets.load(filename=model.path+'/set/m_dict_demand_'      +CaseName+'.csv', set='d'   , format='set')

    model.dd   = Set(initialize=dictSets['d' ],   ordered=False, doc='demands'       )
    model.gg   = Set(initialize=dictSets['g' ],   ordered=False, doc='units'         )
    model.nd   = Set(initialize=dictSets['nd'],   ordered=False, doc='nodes'         )
    model.ni   = Set(initialize=dictSets['nd'],   ordered=False, doc='nodes'         )
    model.nf   = Set(initialize=dictSets['nd'],   ordered=False, doc='nodes'         )
    model.cc   = Set(initialize=dictSets['cc'],   ordered=False, doc='circuits'      )

    # %% parameters
    pSBase               = dfParameter['SBase'        ][0] * 1e-3                                                                           # base power                          [GW]
    pReferenceNode       = dfParameter['ReferenceNode'][0]                                                                                  # reference node

    # %% demand parameters
    pDemToNode           = dfDemand     ['Node'                ]                                                                            # demand
    pLoad                = dfDemand     ['Load'                ] * 1e-3                                                                     # demand                              [GW]
    pBidPrice            = dfDemand     ['Price'               ] * 1e-3                                                                     # demand                              [$/GWh]


    #%% generation parameters
    pGenToNode          = dfGeneration  ['Node'                ]                                                                            # generator location in node
    pMaxPower           = dfGeneration  ['Capacity'            ] * 1e-3                                                                     # Maximum power                               [GW]
    pOfferPrice         = dfGeneration  ['Price'               ] * 1e-3                                                                     # Offer Price                                 [$/GWh]

    # %% network parameters
    pLineB              = dfNetwork     ['Susceptance'         ].sort_index()                                                               # susceptance                                 [p.u.]
    pLineNTCFrw         = dfNetwork     ['TTC'                 ] * 1e-3 * dfNetwork['SecurityFactor' ]                                      # net transfer capacity in forward  direction [GW]
    pLineNTCBck         = dfNetwork     ['TTCBck'              ] * 1e-3 * dfNetwork['SecurityFactor' ]                                      # net transfer capacity in backward direction [GW]

    # replace pLineNTCBck = 0.0 by pLineNTCFrw
    pLineNTCBck = pLineNTCBck.where(pLineNTCBck > 0.0, other=pLineNTCFrw)
    # replace pLineNTCFrw = 0.0 by pLineNTCBck
    pLineNTCFrw = pLineNTCFrw.where(pLineNTCFrw > 0.0, other=pLineNTCBck)

    print(pLineB)
    ReadingDataTime = time.time() - StartTime
    StartTime       = time.time()
    print('Reading    input data                 ... ', round(ReadingDataTime), 's')

    # %% defining subsets:
    model.d = Set(initialize=model.dd, ordered=False, doc='demanding    units',
                  filter=lambda model, dd: dd in model.dd and pLoad[dd] > 0.0)
    model.g = Set(initialize=model.gg, ordered=False, doc='generating    units',
                   filter=lambda model, gg: gg in model.gg and pMaxPower[gg] > 0.0)
    model.la = Set(initialize=model.ni * model.nf * model.cc, ordered=False, doc='all           lines',
                    filter=lambda model, ni, nf, cc: (ni, nf, cc) in pLineB)
    model.rf = Set(initialize=model.nd, ordered=True, doc='reference node',
                    filter=lambda model, nd: nd in pReferenceNode)

    # %% inverse index node to demand
    pNodeToDem = pDemToNode.reset_index().set_index('Node').set_axis(['Demand'], axis=1, inplace=False)[
        ['Demand']]
    pNodeToDem = pNodeToDem.loc[pNodeToDem['Demand'].isin(model.d)]

    pNode2Dem = pNodeToDem.reset_index()
    pNode2Dem['Y/N'] = 1
    pNode2Dem = pNode2Dem.set_index(['Node', 'Demand'])['Y/N']

    model.n2d = Set(initialize=model.nd * model.d, doc='node   to demand',
                     filter=lambda model, nd, d: (nd, d) in pNode2Dem)

    # %% inverse index node to generator
    pNodeToGen = pGenToNode.reset_index().set_index('Node').set_axis(['Generator'], axis=1, inplace=False)[
        ['Generator']]
    pNodeToGen = pNodeToGen.loc[pNodeToGen['Generator'].isin(model.g)]

    pNode2Gen = pNodeToGen.reset_index()
    pNode2Gen['Y/N'] = 1
    pNode2Gen = pNode2Gen.set_index(['Node', 'Generator'])['Y/N']

    model.n2g = Set(initialize=model.nd * model.g, doc='node   to generator',
                    filter=lambda model, nd, g: (nd, g) in pNode2Gen)

    # BigM maximum flow to be used in the Kirchhoff's 2nd law disjunctive constraint
    pBigMFlowBck = pLineNTCBck * 0
    pBigMFlowFrw = pLineNTCFrw * 0
    for la in model.la:
        pBigMFlowBck.loc[la] = pLineNTCBck[la]
        pBigMFlowFrw.loc[la] = pLineNTCFrw[la]

    model.pSBase                = Param(          initialize=pSBase                     , within=NonNegativeReals)

    model.pMaxPower             = Param(model.gg, initialize=pMaxPower.to_dict()        , within=NonNegativeReals, doc='Maximum power'                )
    model.pOfferPrice           = Param(model.gg, initialize=pOfferPrice.to_dict()      , within=NonNegativeReals, doc='Offer price'                  )

    model.pLoad                 = Param(model.dd, initialize=pLoad.to_dict()            , within=NonNegativeReals, doc='Maximum power'                )
    model.pBidPrice             = Param(model.dd, initialize=pBidPrice.to_dict()        , within=NonNegativeReals, doc='Offer price'                  )
   

    model.pLineB                = Param(                               model.la, initialize=pLineB.to_dict()                   , within=NonNegativeReals, doc='Susceptance'                    )
    model.pLineNTCFrw           = Param(                               model.la, initialize=pLineNTCFrw.to_dict()              , within=NonNegativeReals, doc='NTC forward'                    )
    model.pLineNTCBck           = Param(                               model.la, initialize=pLineNTCBck.to_dict()              , within=NonNegativeReals, doc='NTC backward'                   )
    model.pBigMFlowBck          = Param(                               model.la, initialize=pBigMFlowBck.to_dict()             , within=NonNegativeReals, doc='Maximum backward capacity', mutable=True)
    model.pBigMFlowFrw          = Param(                               model.la, initialize=pBigMFlowFrw.to_dict()             , within=NonNegativeReals, doc='Maximum forward  capacity', mutable=True)

    #%% variables
    model.vTotalGCost           = Var(                                      within=NonNegativeReals,                                                                                              doc='total variable generation  operation cost        [$ ]')
    model.vTotalGen             = Var(model.g                             , within=NonNegativeReals, bounds=lambda model,g : (0,model.pMaxPower        [g ]),                                     doc='total output of the gen unit                     [GW]')
    model.vTotalDem             = Var(model.d                             , within=NonNegativeReals, bounds=lambda model,d : (0,model.pLoad            [d ]),                                     doc='total output of the dem unit                     [GW]')

    model.vFlow                 = Var(model.la                            , within=Reals,            bounds=lambda model,*la: (-model.pLineNTCBck[la], model.pLineNTCFrw[la]),                    doc='flow                                             [GW]')
    model.vTheta                = Var(model.nd                            , within=Reals,            bounds=lambda model, nd: (-math.pi/2, math.pi/2),                                            doc='voltage angle                                   [rad]')

    # fixing the voltage angle of the reference node for each scenario, period, and load level
    model.vTheta[model.rf.first()].fix(0.0)

    SettingUpDataTime = time.time() - StartTime
    StartTime         = time.time()
    print('Setting up input data                 ... ', round(SettingUpDataTime), 's')
