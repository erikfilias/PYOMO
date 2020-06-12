# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.32 - May 26, 2020

import time
import pandas        as pd
from   pyomo.environ import DataPortal, Set, Param, Var, Binary, NonNegativeReals, RealSet, UnitInterval, Boolean, RangeSet

print('Input data            ****')

StartTime = time.time()

#%% reading data from CSV
dfOption             = pd.read_csv('oT_Data_Option_'            +CaseName+'.csv', index_col=[0    ])
dfParameter          = pd.read_csv('oT_Data_Parameter_'         +CaseName+'.csv', index_col=[0    ])
dfScenario           = pd.read_csv('oT_Data_Scenario_'          +CaseName+'.csv', index_col=[0    ])
dfDuration           = pd.read_csv('oT_Data_Duration_'          +CaseName+'.csv', index_col=[0    ])
dfDemand             = pd.read_csv('oT_Data_Demand_'            +CaseName+'.csv', index_col=[0,1,2])
dfOperatingReserve   = pd.read_csv('oT_Data_OperatingReserve_'  +CaseName+'.csv', index_col=[0,1,2])
dfGeneration         = pd.read_csv('oT_Data_Generation_'        +CaseName+'.csv', index_col=[0    ])
dfVariableMaxPower   = pd.read_csv('oT_Data_VariableGeneration_'+CaseName+'.csv', index_col=[0,1,2])
dfVariableMinStorage = pd.read_csv('oT_Data_VariableMinStorage_'+CaseName+'.csv', index_col=[0,1,2])
dfVariableMaxStorage = pd.read_csv('oT_Data_VariableMaxStorage_'+CaseName+'.csv', index_col=[0,1,2])
dfESSEnergyInflows   = pd.read_csv('oT_Data_ESSEnergyInflows_'  +CaseName+'.csv', index_col=[0,1,2])
dfNodeLocation       = pd.read_csv('oT_Data_NodeLocation_'      +CaseName+'.csv', index_col=[0    ])
dfNetwork            = pd.read_csv('oT_Data_Network_'           +CaseName+'.csv', index_col=[0,1,2])

# substitute NaN by 0
dfOption.fillna            (0, inplace=True)
dfParameter.fillna         (0, inplace=True)
dfScenario.fillna          (0, inplace=True)
dfDuration.fillna          (0, inplace=True)
dfDemand.fillna            (0, inplace=True)
dfOperatingReserve.fillna  (0, inplace=True)
dfGeneration.fillna        (0, inplace=True)
dfVariableMaxPower.fillna  (0, inplace=True)
dfVariableMinStorage.fillna(0, inplace=True)
dfVariableMaxStorage.fillna(0, inplace=True)
dfESSEnergyInflows.fillna  (0, inplace=True)
dfNodeLocation.fillna      (0, inplace=True)
dfNetwork.fillna           (0, inplace=True)

#%% reading the sets
dictSets = DataPortal()
dictSets.load(filename='oT_Dict_Scenario_'    +CaseName+'.csv', set='sc'  , format='set')
dictSets.load(filename='oT_Dict_Period_'      +CaseName+'.csv', set='p'   , format='set')
dictSets.load(filename='oT_Dict_LoadLevel_'   +CaseName+'.csv', set='n'   , format='set')
dictSets.load(filename='oT_Dict_Generation_'  +CaseName+'.csv', set='g'   , format='set')
dictSets.load(filename='oT_Dict_Technology_'  +CaseName+'.csv', set='gt'  , format='set')
dictSets.load(filename='oT_Dict_Storage_'     +CaseName+'.csv', set='st'  , format='set')
dictSets.load(filename='oT_Dict_Node_'        +CaseName+'.csv', set='nd'  , format='set')
dictSets.load(filename='oT_Dict_Zone_'        +CaseName+'.csv', set='zn'  , format='set')
dictSets.load(filename='oT_Dict_Area_'        +CaseName+'.csv', set='ar'  , format='set')
dictSets.load(filename='oT_Dict_Region_'      +CaseName+'.csv', set='rg'  , format='set')
dictSets.load(filename='oT_Dict_Circuit_'     +CaseName+'.csv', set='cc'  , format='set')
dictSets.load(filename='oT_Dict_Line_'        +CaseName+'.csv', set='lt'  , format='set')

dictSets.load(filename='oT_Dict_NodeToZone_'  +CaseName+'.csv', set='ndzn', format='set')
dictSets.load(filename='oT_Dict_ZoneToArea_'  +CaseName+'.csv', set='znar', format='set')
dictSets.load(filename='oT_Dict_AreaToRegion_'+CaseName+'.csv', set='arrg', format='set')

mTEPES.sc = Set(initialize=dictSets['sc'], ordered=True,  doc='scenarios'   )
mTEPES.p  = Set(initialize=dictSets['p' ], ordered=True,  doc='periods'     )
mTEPES.nn = Set(initialize=dictSets['n' ], ordered=True,  doc='load levels' )
mTEPES.gg = Set(initialize=dictSets['g' ], ordered=False, doc='units'       )
mTEPES.gt = Set(initialize=dictSets['gt'], ordered=False, doc='technologies')
mTEPES.st = Set(initialize=dictSets['st'], ordered=False, doc='ESS types'   )
mTEPES.nd = Set(initialize=dictSets['nd'], ordered=False, doc='nodes'       )
mTEPES.ni = Set(initialize=dictSets['nd'], ordered=False, doc='nodes'       )
mTEPES.nf = Set(initialize=dictSets['nd'], ordered=False, doc='nodes'       )
mTEPES.zn = Set(initialize=dictSets['zn'], ordered=False, doc='zones'       )
mTEPES.ar = Set(initialize=dictSets['ar'], ordered=False, doc='areas'       )
mTEPES.rg = Set(initialize=dictSets['rg'], ordered=False, doc='regions'     )
mTEPES.cc = Set(initialize=dictSets['cc'], ordered=False, doc='circuits'    )
mTEPES.lt = Set(initialize=dictSets['lt'], ordered=False, doc='line types'  )

#%% parameters
pIndBinGenInvest    = dfOption   ['IndBinGenInvest'][0].astype('int')                                                                # Indicator of binary generation expansion decisions, 0 continuous - 1 binary
pIndBinNetInvest    = dfOption   ['IndBinNetInvest'][0].astype('int')                                                                # Indicator of binary network    expansion decisions, 0 continuous - 1 binary
pIndBinGenOperat    = dfOption   ['IndBinGenOperat'][0].astype('int')                                                                # Indicator of binary generation operation decisions, 0 continuous - 1 binary
pIndNetLosses       = dfOption   ['IndNetLosses'   ][0].astype('int')                                                                # Indicator of network losses,                        0 lossless   - 1 ohmic losses
pENSCost            = dfParameter['ENSCost'        ][0] * 1e-3                                                                       # cost of energy not served           [MEUR/GWh]
pCO2Cost            = dfParameter['CO2Cost'        ][0]                                                                              # cost of CO2 emission                [EUR/CO2 ton]
pSBase              = dfParameter['SBase'          ][0] * 1e-3                                                                       # base power                          [GW]
pReferenceNode      = dfParameter['ReferenceNode'  ][0]                                                                              # reference node
pTimeStep           = dfParameter['TimeStep'       ][0].astype('int')                                                                # duration of the unit time step      [h]

pScenProb           = dfScenario          ['Probability'  ]                                                                          # probabilities of scenarios          [p.u.]
pDuration           = dfDuration          ['Duration'     ] * pTimeStep                                                              # duration of load levels             [h]
pDemand             = dfDemand            [list(mTEPES.nd)] * 1e-3                                                                   # demand                              [GW]
pOperReserveUp      = dfOperatingReserve  ['Up'           ] * 1e-3                                                                   # upward   operating reserve          [GW]
pOperReserveDw      = dfOperatingReserve  ['Down'         ] * 1e-3                                                                   # downward operating reserve          [GW]
pVariableMaxPower   = dfVariableMaxPower  [list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable maximum power      [GW]
pVariableMinStorage = dfVariableMinStorage[list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable minimum storage    [TWh]
pVariableMaxStorage = dfVariableMaxStorage[list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable maximum storage    [TWh]
pESSEnergyInflows   = dfESSEnergyInflows  [list(mTEPES.gg)] * 1e-3                                                                   # dynamic energy inflows              [GW]

# compute the demand as the mean over the time step load levels and assign it to active load levels. Idem for operating reserve, variable max power, variable min and max storage capacity and inflows
pDemand             = pDemand.rolling            (pTimeStep).mean()
pOperReserveUp      = pOperReserveUp.rolling     (pTimeStep).mean()
pOperReserveDw      = pOperReserveDw.rolling     (pTimeStep).mean()
pVariableMaxPower   = pVariableMaxPower.rolling  (pTimeStep).mean()
pVariableMinStorage = pVariableMinStorage.rolling(pTimeStep).mean()
pVariableMaxStorage = pVariableMaxStorage.rolling(pTimeStep).mean()
pESSEnergyInflows   = pESSEnergyInflows.rolling  (pTimeStep).mean()

pDemand.fillna            (0, inplace=True)
pOperReserveUp.fillna     (0, inplace=True)
pOperReserveDw.fillna     (0, inplace=True)
pVariableMaxPower.fillna  (0, inplace=True)
pVariableMinStorage.fillna(0, inplace=True)
pVariableMaxStorage.fillna(0, inplace=True)
pESSEnergyInflows.fillna  (0, inplace=True)

if pTimeStep > 1:
    # assign duration 0 to load levels not being considered, active load levels are at the end of every pTimeStep
    for i in range(pTimeStep-2,-1,-1):
        pDuration[range(i,len(mTEPES.nn),pTimeStep)] = 0
    # drop rows with duration 0
    pDemand             = pDemand.loc            [pDemand.index            [range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pOperReserveUp      = pOperReserveUp.loc     [pOperReserveUp.index     [range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pOperReserveDw      = pOperReserveDw.loc     [pOperReserveDw.index     [range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pVariableMaxPower   = pVariableMaxPower.loc  [pVariableMaxPower.index  [range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pVariableMinStorage = pVariableMinStorage.loc[pVariableMinStorage.index[range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pVariableMaxStorage = pVariableMaxStorage.loc[pVariableMaxStorage.index[range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]
    pESSEnergyInflows   = pESSEnergyInflows.loc  [pESSEnergyInflows.index  [range(pTimeStep-1,len(mTEPES.sc*mTEPES.nn),pTimeStep)]]

#%% generation parameters
pGenToNode             = dfGeneration['Node'              ]                                                                          # generator location in node
pGenToTechnology       = dfGeneration['Technology'        ]                                                                          # generator association to technology
pMustRun               = dfGeneration['MustRun'           ]                                                                          # must-run unit                       [Yes]
pRatedMaxPower         = dfGeneration['MaxPower'          ] * 1e-3                                                                   # rated maximum power                 [GW]
pRatedMinPower         = dfGeneration['MinPower'          ] * 1e-3                                                                   # rated minimum power                 [GW]
pLinearVarCost         = dfGeneration['LinearVarCost'     ] * 1e-3 * dfGeneration['FuelCost'] + dfGeneration['OMVarCost'] * 1e-3     # linear   term variable cost         [MEUR/GWh]
pConstantVarCost       = dfGeneration['ConstantVarCost'   ] * 1e-6 * dfGeneration['FuelCost']                                        # constant term variable cost         [MEUR/h]
pStartUpCost           = dfGeneration['StartUpCost'       ] * 1e-6 * dfGeneration['FuelCost']                                        # startup  cost                       [MEUR]
pShutDownCost          = dfGeneration['ShutDownCost'      ] * 1e-6 * dfGeneration['FuelCost']                                        # shutdown cost                       [MEUR]
pRampUp                = dfGeneration['RampUp'            ] * 1e-3                                                                   # ramp up   rate                      [GW/h]
pRampDw                = dfGeneration['RampDown'          ] * 1e-3                                                                   # ramp down rate                      [GW/h]
pCO2EmissionRate       = dfGeneration['CO2EmissionRate'   ] * 1e-3                                                                   # emission  rate                      [t CO2/MWh]
pUpTime                = dfGeneration['UpTime'            ]                                                                          # minimum up   time                   [h]
pDwTime                = dfGeneration['DownTime'          ]                                                                          # minimum down time                   [h]
pGenFixedCost          = dfGeneration['FixedCost'         ] *        dfGeneration['FixedChargeRate']                                 # generation fixed cost               [MEUR]
pIndBinUnitInvest      = dfGeneration['BinaryInvestment'  ]                                                                          # binary unit investment decision     [Yes]
pMaxCharge             = dfGeneration['MaxCharge'         ] * 1e-3                                                                   # maximum ESS charge                  [GW]
pESSInitialInventory   = dfGeneration['InitialStorage'    ] * 1e-3                                                                   # initial ESS storage                 [TWh]
pESSMaxStorageCapacity = dfGeneration['MaxStorageCapacity'] * 1e-3                                                                   # maximum ESS storage capacity        [TWh]
pESSMinStorageCapacity = dfGeneration['MinStorageCapacity'] * 1e-3                                                                   # minimum ESS storage capacity        [TWh]
pEfficiency            = dfGeneration['Efficiency'        ]                                                                          #         ESS efficiency              [p.u.]
pStorageType           = dfGeneration['StorageType'       ]                                                                          #         ESS type

pNodeLat              = dfNodeLocation['Latitude'         ]                                                                          # node latitude                       [º]
pNodeLon              = dfNodeLocation['Longitude'        ]                                                                          # node longitude                      [º]

pLineType              = dfNetwork   ['LineType'          ]                                                                          # line type
pLineVoltage           = dfNetwork   ['Voltage'           ]                                                                          # line voltage                        [kV]
pLineLossFactor        = dfNetwork   ['LossFactor'        ]                                                                          # loss factor                         [p.u.]
pLineR                 = dfNetwork   ['Resistance'        ]                                                                          # resistance                           [p.u.]
pLineX                 = dfNetwork   ['Reactance'         ]                                                                          # reactance                           [p.u.]
pLineBsh               = dfNetwork   ['Susceptance'       ]                                                                          # Susceptance                           [p.u.]
pLineTAP               = dfNetwork   ['TAP'               ]                                                                          # TAP changer                           [p.u.]
pLineNTC               = dfNetwork   ['TTC'               ] * 1e-3 * dfNetwork['SecurityFactor' ]                                    # net transfer capacity               [GW]
pNetFixedCost          = dfNetwork   ['FixedCost'         ] *        dfNetwork['FixedChargeRate']                                    # network    fixed cost               [MEUR]
pIndBinLineInvest      = dfNetwork   ['BinaryInvestment'  ]                                                                          # binary line    investment decision  [Yes]

ReadingDataTime = time.time() - StartTime
StartTime       = time.time()
print('Reading    input data                 ... ', round(ReadingDataTime), 's')

#%% defining subsets: active load levels (n), thermal units (t), ESS units (es), candidate gen units (gc), candidate ESS units (ec), all the lines (la), candidate lines (lc) and lines with losses (ll)
mTEPES.n  = Set(initialize=mTEPES.nn,                     ordered=True , doc='load levels'        , filter=lambda mTEPES,nn      : nn        in mTEPES.nn and pDuration      [nn] >  0)
mTEPES.n2 = Set(initialize=mTEPES.nn,                     ordered=True , doc='load levels'        , filter=lambda mTEPES,nn      : nn        in mTEPES.nn and pDuration      [nn] >  0)
mTEPES.g  = Set(initialize=mTEPES.gg,                     ordered=False, doc='generating    units', filter=lambda mTEPES,gg      : gg        in mTEPES.gg and pRatedMaxPower [gg] >  0)
mTEPES.t  = Set(initialize=mTEPES.g ,                     ordered=False, doc='thermal       units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pLinearVarCost  [g] >  0)
mTEPES.r  = Set(initialize=mTEPES.g ,                     ordered=False, doc='RES           units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pLinearVarCost  [g] == 0 and pESSMaxStorageCapacity[g] == 0)
mTEPES.es = Set(initialize=mTEPES.g ,                     ordered=False, doc='ESS           units', filter=lambda mTEPES,g       : g         in mTEPES.g  and                              pESSMaxStorageCapacity[g] >  0)
mTEPES.gc = Set(initialize=mTEPES.g ,                     ordered=False, doc='candidate     units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pGenFixedCost   [g] >  0)
mTEPES.ec = Set(initialize=mTEPES.es,                     ordered=False, doc='candidate ESS units', filter=lambda mTEPES,es      : es        in mTEPES.es and pGenFixedCost  [es] >  0)
mTEPES.la = Set(initialize=mTEPES.ni*mTEPES.nf*mTEPES.cc, ordered=False, doc='all           lines', filter=lambda mTEPES,ni,nf,cc:(ni,nf,cc) in               pLineX                  )
mTEPES.lc = Set(initialize=mTEPES.la,                     ordered=False, doc='candidate     lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost  [la] >  0)
mTEPES.cd = Set(initialize=mTEPES.la,                     ordered=False, doc='          DC  lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost  [la] >  0 and pLineType[la] == 'DC')
mTEPES.ed = Set(initialize=mTEPES.la,                     ordered=False, doc='          DC  lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost  [la] == 0 and pLineType[la] == 'DC')
mTEPES.ll = Set(initialize=mTEPES.la,                     ordered=False, doc='loss          lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pLineLossFactor[la] >  0 and pIndNetLosses >   0  )
mTEPES.rf = Set(initialize=mTEPES.nd,                     ordered=True , doc='reference node'     , filter=lambda mTEPES,nd      : nd        in               pReferenceNode          )

# non-RES units, they can contribute to the operating reserve and can be committed
mTEPES.nr = mTEPES.g - mTEPES.r

# existing lines (le)
mTEPES.le = mTEPES.la - mTEPES.lc

# candidate and existing lines in AC (lca and lea)
mTEPES.lca = mTEPES.lc - mTEPES.cd
mTEPES.lea = mTEPES.le - mTEPES.ed

# # input lines
# mTEPES.lin = Set(initialize=mTEPES.nf*mTEPES.ni*mTEPES.cc, doc='input line', filter=lambda mTEPES,nf,ni,cc: (ni,nf,cc) in mTEPES.la)
# mTEPES.lil = Set(initialize=mTEPES.nf*mTEPES.ni*mTEPES.cc, doc='input line', filter=lambda mTEPES,nf,ni,cc: (ni,nf,cc) in mTEPES.ll)

# line type
pLineType = pLineType.reset_index()
pLineType['Y/N'] = 1
pLineType = pLineType.set_index(['level_0','level_1','level_2','LineType'])['Y/N']

mTEPES.pLineType = Set(initialize=mTEPES.ni*mTEPES.nf*mTEPES.cc*mTEPES.lt, doc='line type', filter=lambda mTEPES,ni,nf,cc,lt: (ni,nf,cc,lt) in pLineType)

#%% inverse index node to generator
pNodeToGen = pGenToNode.reset_index().set_index('Node').set_axis(['Generator'], axis=1, inplace=False)[['Generator']]
pNodeToGen = pNodeToGen.loc[pNodeToGen['Generator'].isin(mTEPES.g)]

pNode2Gen = pNodeToGen.reset_index()
pNode2Gen['Y/N'] = 1
pNode2Gen = pNode2Gen.set_index(['Node','Generator'])['Y/N']

mTEPES.n2g = Set(initialize=mTEPES.nd*mTEPES.g, doc='node to generator', filter=lambda mTEPES,nd,g: (nd,g) in pNode2Gen)

#%% inverse index generator to technology
mTEPES.t2g = pGenToTechnology.reset_index().set_index('Technology').set_axis(['Generator'], axis=1, inplace=False)['Generator']

# minimum and maximum variable power
pVariableMaxPower = pVariableMaxPower.replace(0, float('nan'))
pMinPower         = pd.DataFrame([pRatedMinPower]*len(pVariableMaxPower.index), index=pd.MultiIndex.from_tuples(pVariableMaxPower.index), columns=pRatedMinPower.index)
pMaxPower         = pd.DataFrame([pRatedMaxPower]*len(pVariableMaxPower.index), index=pd.MultiIndex.from_tuples(pVariableMaxPower.index), columns=pRatedMaxPower.index)
pMaxPower         = pVariableMaxPower.where(pVariableMaxPower < pMaxPower, other=pMaxPower)
pMaxPower2ndBlock = pMaxPower - pMinPower

# minimum and maximum variable storage capacity
pVariableMinStorage = pVariableMinStorage.replace(0, float('nan'))
pVariableMaxStorage = pVariableMaxStorage.replace(0, float('nan'))
pESSMinStorage      = pd.DataFrame([pESSMinStorageCapacity]*len(pVariableMinStorage.index), index=pd.MultiIndex.from_tuples(pVariableMinStorage.index), columns=pESSMinStorageCapacity.index)
pESSMaxStorage      = pd.DataFrame([pESSMaxStorageCapacity]*len(pVariableMaxStorage.index), index=pd.MultiIndex.from_tuples(pVariableMaxStorage.index), columns=pESSMaxStorageCapacity.index)
pESSMinStorage      = pVariableMinStorage.where(pVariableMinStorage > pESSMinStorage, other=pESSMinStorage)
pESSMaxStorage      = pVariableMaxStorage.where(pVariableMaxStorage < pESSMaxStorage, other=pESSMaxStorage)

# minimum up and down time converted to an integer number of time steps
pUpTime = round(pUpTime/pTimeStep).astype('int')
pDwTime = round(pDwTime/pTimeStep).astype('int')

#%% definition of the time-steps leap to observe the stored energy at ESS
pESSTimeStep = (pUpTime*0).astype('int')
for es in mTEPES.es:
    if  pStorageType[es] == 'Daily'  :
        pESSTimeStep[es] = int(  24/pTimeStep)
    if  pStorageType[es] == 'Weekly' :
        pESSTimeStep[es] = int( 168/pTimeStep)
    if  pStorageType[es] == 'Monthly':
        pESSTimeStep[es] = int( 672/pTimeStep)
    if  pStorageType[es] == 'Yearly' :
        pESSTimeStep[es] = int(8736/pTimeStep)

mTEPES.pStorageType = Set(mTEPES.es*mTEPES.st, initialize=pStorageType)

# drop levels with duration 0
pDuration = pDuration.loc[pDuration.index[range(pTimeStep-1,len(mTEPES.nn),pTimeStep)]]

# values < 1e-5 times the maximum system demand are converted to 0
pEpsilon = pDemand.max().sum()*1e-5
pDemand          [pDemand           < pEpsilon] = 0
pOperReserveUp   [pOperReserveUp    < pEpsilon] = 0
pOperReserveDw   [pOperReserveDw    < pEpsilon] = 0
pMinPower        [pMinPower         < pEpsilon] = 0
pMaxPower        [pMaxPower         < pEpsilon] = 0
pMaxPower2ndBlock[pMaxPower2ndBlock < pEpsilon] = 0
pMaxCharge       [pMaxCharge        < pEpsilon] = 0
pESSEnergyInflows[pESSEnergyInflows < pEpsilon] = 0
pESSMinStorage   [pESSMinStorage    < pEpsilon] = 0
pESSMaxStorage   [pESSMaxStorage    < pEpsilon] = 0
pLineNTC         [pLineNTC          < pEpsilon] = 0

# maximum flow is the net transfer capacity but to be used in the Kirchhoff's 2nd law constraint
pMaxFlow = pLineNTC

# maximum voltage angle
# pMaxTheta = pDemand*0 + float('inf')
# pMaxTheta = pDemand*0 + float(1.570795)

# improve boundS in maximum flow and theta
# pIndImproveNetworkModel = 1
# if pIndImproveNetworkModel == 1:
#     import openTEPES_NetworkModel

# this option avoids a warning in the following assignments
pd.options.mode.chained_assignment = None

mTEPES.pIndBinGenInvest     = Param(initialize=pIndBinGenInvest, within=Boolean, mutable=True)
mTEPES.pIndBinNetInvest     = Param(initialize=pIndBinNetInvest, within=Boolean, mutable=True)
mTEPES.pIndBinGenOperat     = Param(initialize=pIndBinGenOperat, within=Boolean, mutable=True)
mTEPES.pIndNetLosses        = Param(initialize=pIndNetLosses   , within=Boolean)

mTEPES.pENSCost             = Param(initialize=pENSCost        , within=NonNegativeReals)
mTEPES.pCO2Cost             = Param(initialize=pCO2Cost        , within=NonNegativeReals)
mTEPES.pSBase               = Param(initialize=pSBase          , within=NonNegativeReals)

# print(pDemand)
# print(pDuration)
# mTEPES.n.pprint()
# pDemand = pDemand.iloc[[pDemand.index            [range(1,len(mTEPES.n))]]] 
pDemand                     = pDemand.iloc[  range(0,len(mTEPES.n))]
pDuration                   = pDuration.iloc[range(0,len(mTEPES.n))]
pOperReserveUp              = pOperReserveUp.iloc[range(0,len(mTEPES.n))]
pOperReserveDw              = pOperReserveDw.iloc[range(0,len(mTEPES.n))]
pMinPower                   = pMinPower.iloc[range(0,len(mTEPES.n))]
pMaxPower                   = pMaxPower.iloc[range(0,len(mTEPES.n))]
pMaxPower2ndBlock           = pMaxPower2ndBlock.iloc[range(0,len(mTEPES.n))]
pESSEnergyInflows           = pESSEnergyInflows.iloc[range(0,len(mTEPES.n))]
pESSMinStorage              = pESSMinStorage.iloc[range(0,len(mTEPES.n))]
pESSMaxStorage              = pESSMaxStorage.iloc[range(0,len(mTEPES.n))]
# ppMaxTheta                  = pMaxTheta.iloc[range(1,len(mTEPES.n))]

pMaxTheta = pDemand*0 + float(1.570795)

# print(pDemand)
# print(pDuration)
# mTEPES.n.pprint()

mTEPES.pDemand              = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, initialize=pDemand.stack().to_dict()          , within=NonNegativeReals, doc='Demand'                       )
mTEPES.pScenProb            = Param(mTEPES.sc,                                initialize=pScenProb.to_dict()                , within=NonNegativeReals, doc='Probability'                  )
mTEPES.pDuration            = Param(                     mTEPES.n,            initialize=pDuration.to_dict()                , within=NonNegativeReals, doc='Duration'                     )
mTEPES.pNodeLon             = Param(                               mTEPES.nd, initialize=pNodeLon.to_dict()                                          , doc='Longitude'                    )
mTEPES.pNodeLat             = Param(                               mTEPES.nd, initialize=pNodeLat.to_dict()                                          , doc='Latitude'                     )
mTEPES.pOperReserveUp       = Param(mTEPES.sc, mTEPES.p, mTEPES.n,            initialize=pOperReserveUp.to_dict()           , within=NonNegativeReals, doc='Upward   operating reserve'   )
mTEPES.pOperReserveDw       = Param(mTEPES.sc, mTEPES.p, mTEPES.n,            initialize=pOperReserveDw.to_dict()           , within=NonNegativeReals, doc='Downward operating reserve'   )
mTEPES.pMinPower            = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pMinPower.stack().to_dict()        , within=NonNegativeReals, doc='Minimum power'                )
mTEPES.pMaxPower            = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pMaxPower.stack().to_dict()        , within=NonNegativeReals, doc='Maximum power'                )
mTEPES.pMaxPower2ndBlock    = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pMaxPower2ndBlock.stack().to_dict(), within=NonNegativeReals, doc='Second block'                 )
mTEPES.pESSEnergyInflows    = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pESSEnergyInflows.stack().to_dict(), within=NonNegativeReals, doc='ESS Energy inflows'           )
mTEPES.pESSMinStorage       = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pESSMinStorage.stack().to_dict()   , within=NonNegativeReals, doc='ESS Minimum stoarage capacity')
mTEPES.pESSMaxStorage       = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.gg, initialize=pESSMaxStorage.stack().to_dict()   , within=NonNegativeReals, doc='ESS Maximum stoarage capacity')

mTEPES.pLinearVarCost       = Param(                               mTEPES.gg, initialize=pLinearVarCost.to_dict()           , within=NonNegativeReals, doc='Linear   variable cost'       )
mTEPES.pConstantVarCost     = Param(                               mTEPES.gg, initialize=pConstantVarCost.to_dict()         , within=NonNegativeReals, doc='Constant variable cost'       )
mTEPES.pStartUpCost         = Param(                               mTEPES.gg, initialize=pStartUpCost.to_dict()             , within=NonNegativeReals, doc='Startup  cost'                )
mTEPES.pShutDownCost        = Param(                               mTEPES.gg, initialize=pShutDownCost.to_dict()            , within=NonNegativeReals, doc='Shutdown cost'                )
mTEPES.pRampUp              = Param(                               mTEPES.gg, initialize=pRampUp.to_dict()                  , within=NonNegativeReals, doc='Ramp up   rate'               )
mTEPES.pRampDw              = Param(                               mTEPES.gg, initialize=pRampDw.to_dict()                  , within=NonNegativeReals, doc='Ramp down rate'               )
mTEPES.pUpTime              = Param(                               mTEPES.gg, initialize=pUpTime.to_dict()                  , within=NonNegativeReals, doc='Up   time'                    )
mTEPES.pDwTime              = Param(                               mTEPES.gg, initialize=pDwTime.to_dict()                  , within=NonNegativeReals, doc='Down time'                    )
mTEPES.pMaxCharge           = Param(                               mTEPES.gg, initialize=pMaxCharge.to_dict()               , within=NonNegativeReals, doc='Maximum charge power'         )
mTEPES.pGenFixedCost        = Param(                               mTEPES.gg, initialize=pGenFixedCost.to_dict()            , within=NonNegativeReals, doc='Generation fixed cost'        )
mTEPES.pIndBinUnitInvest    = Param(                               mTEPES.gg, initialize=pIndBinUnitInvest.to_dict()        , within=NonNegativeReals, doc='Binary investment decision'   )
mTEPES.pEfficiency          = Param(                               mTEPES.gg, initialize=pEfficiency.to_dict()              , within=NonNegativeReals, doc='Round-trip efficiency'        )
mTEPES.pCO2EmissionRate     = Param(                               mTEPES.gg, initialize=pCO2EmissionRate.to_dict()         , within=NonNegativeReals, doc='CO2 Emission rate'            )
mTEPES.pESSTimeStep         = Param(                               mTEPES.gg, initialize=pESSTimeStep.to_dict()             , within=NonNegativeReals, doc='ESS Storage cycle'            )
mTEPES.pESSInitialInventory = Param(                               mTEPES.gg, initialize=pESSInitialInventory.to_dict()     , within=NonNegativeReals, doc='ESS Initial storage', mutable=True)

mTEPES.pLineLossFactor      = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineLossFactor.to_dict()          , within=NonNegativeReals, doc='Loss factor'                  )
mTEPES.pLineX               = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineX.to_dict()                   , within=NonNegativeReals, doc='Reactance'                    )
mTEPES.pLineVoltage         = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineVoltage.to_dict()             , within=NonNegativeReals, doc='Voltage'                      )
mTEPES.pLineNTC             = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineNTC.to_dict()                 , within=NonNegativeReals, doc='NTC'                          )
mTEPES.pNetFixedCost        = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pNetFixedCost.to_dict()            , within=NonNegativeReals, doc='Network fixed cost'           )
mTEPES.pIndBinLineInvest    = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pIndBinLineInvest.to_dict()        , within=NonNegativeReals, doc='Binary investment decision'   )
mTEPES.pMaxFlow             = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pMaxFlow.to_dict()                 , within=NonNegativeReals, doc='Maximum capacity'             )
mTEPES.pMaxTheta            = Param(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ni, initialize=pMaxTheta.stack().to_dict()        , within=NonNegativeReals, doc='Maximum voltage angle'        )

#%% variables
mTEPES.vTotalFCost           = Var(                                          within=NonNegativeReals,                                                                                                     doc='total system fixed    cost                     [MEUR]')
mTEPES.vTotalVCost           = Var(                                          within=NonNegativeReals,                                                                                                     doc='total system variable cost                     [MEUR]')
mTEPES.vTotalECost           = Var(                                          within=NonNegativeReals,                                                                                                     doc='total system emission cost                     [MEUR]')
mTEPES.vTotalOutput          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g , within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,g :(0,mTEPES.pMaxPower        [sc,p,n,g ]),                             doc='total output of the unit                         [GW]')
mTEPES.vOutput2ndBlock       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,mTEPES.pMaxPower2ndBlock[sc,p,n,nr]),                             doc='second block of the unit                         [GW]')
mTEPES.vReserveUp            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,mTEPES.pMaxPower2ndBlock[sc,p,n,nr]),                             doc='upward   operating reserve                       [GW]')
mTEPES.vReserveDown          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,mTEPES.pMaxPower2ndBlock[sc,p,n,nr]),                             doc='downward operating reserve                       [GW]')
mTEPES.vESSInventory         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,es:(mTEPES.pESSMinStorage[sc,p,n,es],mTEPES.pESSMaxStorage[sc,p,n,es]), doc='ESS inventory                                   [TWh]')
mTEPES.vESSSpillage          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals,                                                                                                     doc='ESS spillage                                    [TWh]')
mTEPES.vESSCharge            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,es:(0,mTEPES.pMaxCharge       [es]       ),                             doc='ESS charge power                                 [GW]')
mTEPES.vENS                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nd:(0,1),                             doc='energy not served in node                        [GW]')

if mTEPES.pIndBinGenOperat == 0:
    mTEPES.vCommitment       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                               doc='commitment of the unit                          [0,1]')
    mTEPES.vStartUp          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                               doc='startup    of the unit                          [0,1]')
    mTEPES.vShutDown         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                               doc='shutdown   of the unit                          [0,1]')
else:
    mTEPES.vCommitment       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                               doc='commitment of the unit                          {0,1}')
    mTEPES.vStartUp          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                               doc='startup    of the unit                          {0,1}')
    mTEPES.vShutDown         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                               doc='shutdown   of the unit                          {0,1}')

if mTEPES.pIndBinGenInvest == 0:
    mTEPES.vGenerationInvest = Var(                               mTEPES.gc, within=NonNegativeReals, bounds=lambda mTEPES,       g :(0,1),                                                               doc='generation investment decision exists in a year [0,1]')
else:
    mTEPES.vGenerationInvest = Var(                               mTEPES.gc, within=Binary,                                                                                                               doc='generation investment decision exists in a year {0,1}')

if mTEPES.pIndBinNetInvest == 0:
    mTEPES.vNetworkInvest    = Var(                               mTEPES.lc, within=NonNegativeReals, bounds=lambda mTEPES,      *lc:(0,1),                                                               doc='network    investment decision exists in a year [0,1]')
else:
    mTEPES.vNetworkInvest    = Var(                               mTEPES.lc, within=Binary,                                                                                                               doc='network    investment decision exists in a year {0,1}')

# relax binary condition in generation and network investment decisions
for gc in mTEPES.gc:
    if mTEPES.pIndBinGenInvest != 0 and not (mTEPES.pIndBinUnitInvest[gc] == 'Yes' or mTEPES.pIndBinUnitInvest[gc] == 'YES' or mTEPES.pIndBinUnitInvest[gc] == 'yes' or mTEPES.pIndBinUnitInvest[gc] == 'Y' or mTEPES.pIndBinUnitInvest[gc] == 'y'):
        mTEPES.vGenerationInvest[gc].domain = UnitInterval
for lc in mTEPES.lc:
    if mTEPES.pIndBinNetInvest != 0 and not (mTEPES.pIndBinLineInvest[lc] == 'Yes' or mTEPES.pIndBinLineInvest[lc] == 'YES' or mTEPES.pIndBinLineInvest[lc] == 'yes' or mTEPES.pIndBinLineInvest[lc] == 'Y' or mTEPES.pIndBinLineInvest[lc] == 'y'):
        mTEPES.vNetworkInvest   [lc].domain = UnitInterval

# mTEPES.vLineLosses           = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,*ll:(0,0.5*mTEPES.pLineLossFactor[ll]*mTEPES.pLineNTC[ll]),     doc='half line losses                                 [GW]')
# mTEPES.vFlow                 = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,          bounds=lambda mTEPES,sc,p,n,*la:(-mTEPES.pLineNTC[la],mTEPES.pLineNTC[la]),                 doc='flow                                             [GW]')
mTEPES.vTheta                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=RealSet,          bounds=lambda mTEPES,sc,p,n, nd:(-mTEPES.pMaxTheta[sc,p,n,nd],mTEPES.pMaxTheta[sc,p,n,nd]), doc='voltage angle                                   [rad]')

for nr in mTEPES.nr:
    if pMustRun[nr] == 'Yes' or pMustRun[nr] == 'YES' or pMustRun[nr] == 'yes' or pMustRun[nr] == 'Y' or pMustRun[nr] == 'y':
        pMustRun[nr] = 1
    else:
        pMustRun[nr] = 0

for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
    if pMustRun[nr] == 1:
        mTEPES.vCommitment[sc,p,n,nr].fix(1)

mTEPES.pMustRun  = Param(mTEPES.gg, initialize=pMustRun.to_dict())

# fix the must-run units and their output
mTEPES.pInitialOutput = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
mTEPES.pInitialUC     = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
pSystemOutput  = 0
for nr in mTEPES.nr:
    n1 = next(iter(mTEPES.sc*mTEPES.p*mTEPES.n))
    if pSystemOutput < sum(pDemand[nd][n1] for nd in mTEPES.nd) and mTEPES.pMustRun[nr] == 1:
        mTEPES.pInitialOutput[nr] = mTEPES.pMaxPower[n1,nr]
        mTEPES.pInitialUC    [nr] = 1
        pSystemOutput      = pSystemOutput + mTEPES.pInitialOutput[nr]

# thermal and variable units ordered by increasing variable cost
mTEPES.go = pLinearVarCost.sort_values().index

# determine the initial committed units and their output
for go in mTEPES.go:
    n1 = next(iter(mTEPES.sc*mTEPES.p*mTEPES.n))
    if pSystemOutput < sum(pDemand[nd][n1] for nd in mTEPES.nd) and mTEPES.pMustRun[go] != 1:
        mTEPES.pInitialOutput[go] = mTEPES.pMaxPower[n1,go]
        mTEPES.pInitialUC    [go] = 1
        pSystemOutput      = pSystemOutput + mTEPES.pInitialOutput[go]

# fixing the voltage angle of the reference node for each scenario, period, and load level
for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
    mTEPES.vTheta[sc,p,n,mTEPES.rf.first()].fix(0)

# fixing the ESS inventory at the last load level
for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
    mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(pESSInitialInventory[es])

# fixing the ESS inventory at the end of the following pESSTimeStep (weekly, yearly), i.e., for daily ESS is fixed at the end of the week, for weekly/monthly ESS is fixed at the end of the year
for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
     if pStorageType[es] == 'Daily'   and mTEPES.n.ord(n) % ( 168/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Weekly'  and mTEPES.n.ord(n) % (8736/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Monthly' and mTEPES.n.ord(n) % (8736/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])



#%% AC Power Flow: Additional Sets
mTEPES.L                       = RangeSet(20)

#%% AC Power Flow: Additional Parameters
pVmin                          = 0.95
pVnom                          = 1.00
pVmax                          = 1.05
pCapacitivePF                  = 0.95                #Capacitive Power Factor
pInductivePF                   = 0.99                #Inductive Power Factor

pBusGshb                       = dfNodeLocation.drop(['Latitude','Longitude'], axis=1).assign(gshb=0)
pBusBshb                       = dfNodeLocation.drop(['Latitude','Longitude'], axis=1).assign(bshb=0)

pLineBsh                       = pLineX * 0
pLineTAP                       = pLineX * 0 + 1
pLineFi                        = pLineX * 0



# pLineR                         = pLineX * 0 + 0.01
pLineSmax                      = pLineNTC * 1.5

pLineZ2                        = pLineR**2 + pLineX**2
pLineBsh                       = pLineBsh/2
pLineTAP                       = 1/pLineTAP
pLineFi                        = (pLineFi*3.14159265359)/180

#%% Parameters: Branches
mTEPES.pBusGshb                = Param(         mTEPES.nd, initialize=pBusGshb['gshb'].to_dict()                 , within=NonNegativeReals, doc='Conductance'                   )
mTEPES.pBusBshb                = Param(         mTEPES.nd, initialize=pBusBshb['bshb'].to_dict()                 , within=NonNegativeReals, doc='Susceptance'                   )

#%% Parameters: Branches

mTEPES.pLineR                  = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineR.to_dict()                 , within=NonNegativeReals, doc='Resistance'                   )
mTEPES.pLineSmax               = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineSmax.to_dict()              , within=NonNegativeReals, doc='Apparent power capacity'      )
mTEPES.pLineZ2                 = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineZ2.to_dict()                , within=NonNegativeReals, doc='Squared impedance'            )
mTEPES.pLineBsh                = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineBsh.to_dict()               , within=NonNegativeReals, doc='Susceptance'                  )
mTEPES.pLineTAP                = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineTAP.to_dict()               , within=NonNegativeReals, doc='Tap changer'                  )
mTEPES.pLineFi                 = Param(         mTEPES.ni, mTEPES.nf, mTEPES.cc, initialize=pLineFi.to_dict()                , within=NonNegativeReals, doc='Phase shifter'                )

#%% Parameters: Power System
mTEPES.pVmin                   = Param(         initialize=pVmin                , within=NonNegativeReals, doc='Minimum voltage magnitude'                )
mTEPES.pVnom                   = Param(         initialize=pVnom                , within=NonNegativeReals, doc='Nominal voltage magnitude'                )
mTEPES.pVmax                   = Param(         initialize=pVmax                , within=NonNegativeReals, doc='Maximum voltage magnitude'                )
mTEPES.pCapacitivePF           = Param(         initialize=pCapacitivePF        , within=NonNegativeReals, doc='Capacitive power factor'                  )
mTEPES.pInductivePF            = Param(         initialize=pInductivePF         , within=NonNegativeReals, doc='Inductive power factor'                   )

#%% Parameters: Linearization
def delta_S_init(mTEPES,i,j,cc):
    a = (pLineSmax[i,j,cc])/len(mTEPES.L)
    return a
mTEPES.pLineDelta_S            = Param(         mTEPES.la, initialize = delta_S_init                   , within=NonNegativeReals, doc='Delta of Smax splitted by L'  )

def m_init(mTEPES,i,j,cc,l):
    a = (2*l-1)*(mTEPES.pLineDelta_S[i,j,cc])
    return a
mTEPES.pLineM                  = Param( mTEPES.la,mTEPES.L, initialize = m_init                   , within=NonNegativeReals, doc='M partitions of Delta Smax'   )


#%% AC Power Flow: Additional Variables
# def _bounds_Delta_S_rule(mTEPES, sc, p, n, la, L):
#     return (0, mTEPES.pLineDelta_S[la])
mTEPES.vDelta_P                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,ni,nf,cc,l:(0, mTEPES.pLineDelta_S[ni,nf,cc]),          doc='Delta Active Power Flow                                       [  GW]')
mTEPES.vDelta_Q                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, mTEPES.L, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,ni,nf,cc,l:(0, mTEPES.pLineDelta_S[ni,nf,cc]),          doc='Delta Reactive Power Flow                                     [GVaR]')
mTEPES.vP_max                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    doc='Maximum bound of Active Power Flow                            [  GW]')
mTEPES.vP_min                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    doc='Minimum bound of Active Power Flow                            [  GW]')
mTEPES.vQ_max                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    doc='Maximum bound of Reactive Power Flow                          [GVaR]')
mTEPES.vQ_min                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    doc='Minimum bound of Reactive Power Flow                          [GVaR]')

# #VARIABLES
mTEPES.vVoltageMag_sqr         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,nd:(pVmin**2,pVmax**2),          doc='Voltage magnitude                                   [p.u.]')
mTEPES.vCurrentFlow_sqr        = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=NonNegativeReals,    bounds=lambda mTEPES,sc,p,n,*la:(0,(pLineNTC[la]**2/pVmax)),         doc='Current flow                                           [A]')
mTEPES.vP                      = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,             doc='Active Power Flow through lines                                 [GW]')
mTEPES.vQ                      = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,             doc='Reactive Power Flow through lines                               [GW]')
mTEPES.vReactiveTotalOutput    = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g , within=RealSet,             bounds=lambda mTEPES,sc,p,n,g :(-mTEPES.pMaxPower        [sc,p,n,g ],mTEPES.pMaxPower        [sc,p,n,g ]), doc='Total output of reactive power generators                     [GVAr]')

SettingUpDataTime = time.time() - StartTime
StartTime         = time.time()
print('Setting up input data                 ... ', round(SettingUpDataTime), 's')

import numpy as np
pLineZ = pLineR + pLineX*1j
pLineY = 1/pLineZ
Yb     = np.zeros((len(mTEPES.nd), len(mTEPES.nd)), dtype=complex)

dfPosition = pd.DataFrame(columns=['Position'], index = mTEPES.nd)
Count = 0
for i in mTEPES.nd:
    dfPosition['Position'][i] = Count
    Count += 1

mTEPES.dfPosition = dfPosition
print(Yb)
for ni,nf,cc in mTEPES.la:
    Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]] = Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]] + pLineY[ni,nf,cc]*pLineTAP[ni,nf,cc]
    Yb[dfPosition['Position'][nf], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni], dfPosition['Position'][nf]]

for k in mTEPES.nd:
    for ni,nf,cc in mTEPES.la:
        if ni == k:
            Yb[dfPosition['Position'][ni], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni],dfPosition['Position'][ni]] + pLineY[ni,nf,cc]*pLineTAP[ni,nf,cc]**2 + pLineBsh[ni,nf,cc]*1j
        elif nf == k:
            Yb[dfPosition['Position'][ni], dfPosition['Position'][ni]] = Yb[dfPosition['Position'][ni],dfPosition['Position'][ni]] + pLineY[ni,nf,cc]                       + pLineBsh[ni, nf, cc] * 1j

for k in mTEPES.nd:
    Yb[dfPosition['Position'][k], dfPosition['Position'][k]] = Yb[dfPosition['Position'][k], dfPosition['Position'][k]] + mTEPES.pBusBshb[k] * 1j

mTEPES.Yb = Yb
print(Yb)

import numpy as np
Ybarra                      = np.abs(mTEPES.Yb)/np.abs(mTEPES.Yb)
Ybarra[np.isnan(Ybarra)]    = 0

from cvxopt import spmatrix, amd
from scipy.sparse import csr_matrix, find
Ybarra                      = csr_matrix(Ybarra.real)
#find(Ybarra)
coo = Ybarra.tocoo()
SP = spmatrix(coo.data, coo.row.tolist(), coo.col.tolist())
isinstance(SP,spmatrix)

from cvxopt import spmatrix, amd
reordening = amd.order(SP)
# print(reordening)
Yorden                      = np.zeros((len(mTEPES.nd), len(mTEPES.nd)), dtype=complex)
for ni in mTEPES.nd:
    for nf in mTEPES.nd:
        Yorden[mTEPES.dfPosition['Position'][ni], mTEPES.dfPosition['Position'][nf]] = mTEPES.Yb[reordening[mTEPES.dfPosition['Position'][ni]],reordening[mTEPES.dfPosition['Position'][nf]]]

print(Yorden)

# from cvxopt import spmatrix, amd, normal
# from chompack import symbolic, cspmatrix, cholesky
# pattern=cholesky(abs(Yorden.imag))