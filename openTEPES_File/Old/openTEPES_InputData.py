#%% libraries
import pandas        as pd
import time          # count clock time
import psutil        # access the number of CPUs
import pyomo.environ as pyo
from   pyomo.environ import Set, Var, Binary, NonNegativeReals, RealSet, UnitInterval, Constraint, ConcreteModel, Objective, minimize, Suffix, DataPortal, TerminationCondition
from   pyomo.opt     import SolverFactory
from   collections   import defaultdict

StartTime = time.time()

#%% reading the sets
dictSets = DataPortal()
dictSets.load(filename='oT_Dict_Scenario_'    +CaseName+'.csv', set='sc'  , format='set')
dictSets.load(filename='oT_Dict_Period_'      +CaseName+'.csv', set='p'   , format='set')
dictSets.load(filename='oT_Dict_LoadLevel_'   +CaseName+'.csv', set='n'   , format='set')
dictSets.load(filename='oT_Dict_Generation_'  +CaseName+'.csv', set='g'   , format='set')
dictSets.load(filename='oT_Dict_Storage_'     +CaseName+'.csv', set='st'  , format='set')
dictSets.load(filename='oT_Dict_Technology_'  +CaseName+'.csv', set='gt'  , format='set')
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

#%% reading data from CSV
dfOption             = pd.read_csv('oT_Data_Option_'            +CaseName+'.csv', index_col=[0    ])
dfParameter          = pd.read_csv('oT_Data_Parameter_'         +CaseName+'.csv', index_col=[0    ])
dfDuration           = pd.read_csv('oT_Data_Duration_'          +CaseName+'.csv', index_col=[0]    )
dfScenario           = pd.read_csv('oT_Data_Scenario_'          +CaseName+'.csv', index_col=[0    ])
dfNodeLocation       = pd.read_csv('oT_Data_NodeLocation_'      +CaseName+'.csv', index_col=[0    ])
dfDuration           = pd.read_csv('oT_Data_Duration_'          +CaseName+'.csv', index_col=[0    ])
dfDemand             = pd.read_csv('oT_Data_Demand_'            +CaseName+'.csv', index_col=[0,1,2])
dfOperatingReserve   = pd.read_csv('oT_Data_OperatingReserve_'  +CaseName+'.csv', index_col=[0,1,2])
dfGeneration         = pd.read_csv('oT_Data_Generation_'        +CaseName+'.csv', index_col=[0    ])
dfVariableMaxPower   = pd.read_csv('oT_Data_VariableGeneration_'+CaseName+'.csv', index_col=[0,1,2])
dfVariableMinStorage = pd.read_csv('oT_Data_VariableMinStorage_'+CaseName+'.csv', index_col=[0,1,2])
dfVariableMaxStorage = pd.read_csv('oT_Data_VariableMaxStorage_'+CaseName+'.csv', index_col=[0,1,2])
dfESSEnergyInflows   = pd.read_csv('oT_Data_ESSEnergyInflows_'  +CaseName+'.csv', index_col=[0,1,2])
dfNetwork            = pd.read_csv('oT_Data_Network_'           +CaseName+'.csv', index_col=[0,1,2])

# substitute NaN by 0
dfOption.fillna            (0, inplace=True)
dfParameter.fillna         (0, inplace=True)
dfDuration.fillna          (0, inplace=True)
dfScenario.fillna          (0, inplace=True)
dfNodeLocation.fillna      (0, inplace=True)
dfDuration.fillna          (0, inplace=True)
dfDemand.fillna            (0, inplace=True)
dfOperatingReserve.fillna  (0, inplace=True)
dfGeneration.fillna        (0, inplace=True)
dfVariableMaxPower.fillna  (0, inplace=True)
dfVariableMinStorage.fillna(0, inplace=True)
dfVariableMaxStorage.fillna(0, inplace=True)
dfESSEnergyInflows.fillna  (0, inplace=True)
dfNetwork.fillna           (0, inplace=True)

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

pDuration           = dfDuration          ['Duration'     ] * pTimeStep                                                              # duration of load levels             [h]
pScenProb           = dfScenario          ['Probability'  ]                                                                          # probabilities of scenarios          [p.u.]
pDemand             = dfDemand            [list(mTEPES.nd)] * 1e-3                                                                   # demand                              [GW]
pOperReserveUp      = dfOperatingReserve  ['Up'           ] * 1e-3                                                                   # operating reserve up                [GW]
pOperReserveDw      = dfOperatingReserve  ['Down'         ] * 1e-3                                                                   # operating reserve down              [GW]
pVariableMaxPower   = dfVariableMaxPower  [list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable maximum power      [GW]
pVariableMinStorage = dfVariableMinStorage[list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable minimum storage    [TWh]
pVariableMaxStorage = dfVariableMaxStorage[list(mTEPES.gg)] * 1e-3                                                                   # dynamic variable maximum storage    [TWh]
pESSEnergyInflows   = dfESSEnergyInflows  [list(mTEPES.gg)] * 1e-3                                                                   # dynamic energy inflows              [GW]
pNodeLat            = dfNodeLocation      ['Latitude'     ]                                                                          # node latitude                       [º]
pNodeLon            = dfNodeLocation      ['Longitude'    ]                                                                          # node longitude                      [º]

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

pLineType              = dfNetwork   ['LineType'          ]                                                                          # line type
pLineVoltage           = dfNetwork   ['Voltage'           ]                                                                          # line voltage                        [kV]
pLineLossFactor        = dfNetwork   ['LossFactor'        ]                                                                          # loss factor                         [p.u.]
pLineX                 = dfNetwork   ['Reactance'         ]                                                                          # reactance                           [p.u.]
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

# non-RES units
mTEPES.nr = mTEPES.g - mTEPES.r

# existing lines (le)
mTEPES.le = mTEPES.la - mTEPES.lc

# candidate and existing lines in AC (lca) (lea)
mTEPES.lca = mTEPES.lc - mTEPES.cd
mTEPES.lea = mTEPES.le - mTEPES.ed

#%% inverse index node to generator
pNodeToGen = pGenToNode.reset_index().set_index('Node').set_axis(['Generator'], axis=1, inplace=False)[['Generator']]
pNodeToGen = pNodeToGen.loc[pNodeToGen['Generator'].isin(mTEPES.g )]
pNodeToESS = pNodeToGen.loc[pNodeToGen['Generator'].isin(mTEPES.es)]

pNode2Gen = defaultdict(list)
pNode2ESS = defaultdict(list)
for nd in mTEPES.nd:
    if nd in pNodeToGen['Generator']:
        pNode2Gen[nd] = pNodeToGen.loc[[nd],'Generator'].tolist()
    if nd in pNodeToESS['Generator']:
        pNode2ESS[nd] = pNodeToESS.loc[[nd],'Generator'].tolist()

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

# values < 10 kW/kWh are converted to 0
pDemand          [pDemand           < 10e-6] = 0
pOperReserveUp   [pOperReserveUp    < 10e-6] = 0
pOperReserveDw   [pOperReserveDw    < 10e-6] = 0
pMinPower        [pMinPower         < 10e-6] = 0
pMaxPower        [pMaxPower         < 10e-6] = 0
pMaxPower2ndBlock[pMaxPower2ndBlock < 10e-6] = 0
pMaxCharge       [pMaxCharge        < 10e-6] = 0
pESSEnergyInflows[pESSEnergyInflows < 10e-6] = 0
pESSMinStorage   [pESSMinStorage    < 10e-9] = 0
pESSMaxStorage   [pESSMaxStorage    < 10e-9] = 0
pLineNTC         [pLineNTC          < 10e-6] = 0

# maximum flow is the net transfer capacity but to be used in the Kirchhoff's 2nd law constraint
pMaxFlow = pLineNTC

# maximum voltage angle
pMaxTheta = pDemand*0 + float('inf')

# this option avoids a warning in the following assignments
pd.options.mode.chained_assignment = None

#%% variables
mTEPES.vTotalFCost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system fixed    cost                     [MEUR]')
mTEPES.vTotalVCost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system variable cost                     [MEUR]')
mTEPES.vTotalECost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system emission cost                     [MEUR]')
mTEPES.vTotalOutput          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g , within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,g :(0,pMaxPower          [g ][sc,p,n]),                     doc='total output of the unit                         [GW]')
mTEPES.vOutput2ndBlock       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='second block of the unit                         [GW]')
mTEPES.vReserveUp            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='operating reserve up                             [GW]')
mTEPES.vReserveDown          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='operating reserve down                           [GW]')
mTEPES.vESSInventory         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,es:(pESSMinStorage[es][sc,p,n],pESSMaxStorage[es][sc,p,n]), doc='ESS inventory                                   [TWh]')
mTEPES.vESSSpillage          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals,                                                                                         doc='ESS spillage                                    [TWh]')
mTEPES.vESSCharge            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,es:(0,pMaxCharge         [es]        ),                     doc='ESS charge power                                 [GW]')
mTEPES.vENS                  = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nd:(0,pDemand            [nd][sc,p,n]),                     doc='energy not served in node                        [GW]')

if pIndBinGenOperat == 0:
    mTEPES.vCommitment       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                   doc='commitment of the unit                          [0,1]')
    mTEPES.vStartUp          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                   doc='startup    of the unit                          [0,1]')
    mTEPES.vShutDown         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,1),                                                   doc='shutdown   of the unit                          [0,1]')
else:
    mTEPES.vCommitment       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                   doc='commitment of the unit                          {0,1}')
    mTEPES.vStartUp          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                   doc='startup    of the unit                          {0,1}')
    mTEPES.vShutDown         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=Binary,                                                                                                   doc='shutdown   of the unit                          {0,1}')

if pIndBinGenInvest == 0:
    mTEPES.vGenerationInvest = Var(                               mTEPES.gc, within=NonNegativeReals, bounds=lambda mTEPES,       g :(0,1),                                                   doc='generation investment decision exists in a year [0,1]')
else:
    mTEPES.vGenerationInvest = Var(                               mTEPES.gc, within=Binary,                                                                                                   doc='generation investment decision exists in a year {0,1}')

if pIndBinNetInvest == 0:
    mTEPES.vNetworkInvest    = Var(                               mTEPES.lc, within=NonNegativeReals, bounds=lambda mTEPES,      *lc:(0,1),                                                   doc='network    investment decision exists in a year [0,1]')
else:
    mTEPES.vNetworkInvest    = Var(                               mTEPES.lc, within=Binary,                                                                                                   doc='network    investment decision exists in a year {0,1}')

# relax binary condition in generation and network investment decisions
for gc in mTEPES.gc:
    if pIndBinGenInvest != 0 and not (pIndBinUnitInvest[gc] == 'Yes' or pIndBinUnitInvest[gc] == 'YES' or pIndBinUnitInvest[gc] == 'yes' or pIndBinUnitInvest[gc] == 'Y' or pIndBinUnitInvest[gc] == 'y'):
        mTEPES.vGenerationInvest[gc].domain = UnitInterval
for lc in mTEPES.lc:
    if pIndBinNetInvest != 0 and not (pIndBinLineInvest[lc] == 'Yes' or pIndBinLineInvest[lc] == 'YES' or pIndBinLineInvest[lc] == 'yes' or pIndBinLineInvest[lc] == 'Y' or pIndBinLineInvest[lc] == 'y'):
        mTEPES.vNetworkInvest   [lc].domain = UnitInterval

mTEPES.vLineLosses           = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,*ll:(0,0.5*pLineLossFactor[ll]*pLineNTC[ll]),               doc='half line losses                                 [GW]')
mTEPES.vFlow                 = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,          bounds=lambda mTEPES,sc,p,n,*la:(-pLineNTC[la],pLineNTC[la]),                           doc='flow                                             [GW]')
mTEPES.vTheta                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=RealSet,          bounds=lambda mTEPES,sc,p,n,nd:(-pMaxTheta[nd][sc,p,n],pMaxTheta[nd][sc,p,n]),          doc='voltage angle                                   [rad]')

# minimum up and down time converted to an integer number of time steps
pUpTime = round(pUpTime/pTimeStep).astype('int')
pDwTime = round(pDwTime/pTimeStep).astype('int')

for nr in mTEPES.nr:
    if pMustRun[nr] == 'Yes' or pMustRun[nr] == 'YES' or pMustRun[nr] == 'yes' or pMustRun[nr] == 'Y' or pMustRun[nr] == 'y':
        pMustRun[nr] = 1
    else:
        pMustRun[nr] = 0

# fixing the UC of must-run units to 1
for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr:
    if pMustRun[nr] == 1:
        mTEPES.vCommitment[sc,p,n,nr].fix(1)

# fix the must-run units and their output
pInitialOutput = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
pInitialUC     = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
pSystemOutput  = 0
for nr in mTEPES.nr:
    n1 = next(iter(mTEPES.sc*mTEPES.p*mTEPES.n))
    if pSystemOutput < sum(pDemand[nd][n1] for nd in mTEPES.nd) and pMustRun[nr] == 1:
        pInitialOutput[nr] = pMaxPower[nr][n1]
        pInitialUC    [nr] = 1
        pSystemOutput      = pSystemOutput + pInitialOutput[nr]

# thermal and variable units ordered by increasing variable cost
mTEPES.go = pLinearVarCost.sort_values().index

# determine the initial committed units and their output
for go in mTEPES.go:
    n1 = next(iter(mTEPES.sc*mTEPES.p*mTEPES.n))
    if pSystemOutput < sum(pDemand[nd][n1] for nd in mTEPES.nd) and pMustRun[go] != 1:
        pInitialOutput[go] = pMaxPower[go][n1]
        pInitialUC    [go] = 1
        pSystemOutput      = pSystemOutput + pInitialOutput[go]

# fixing the voltage angle of the reference node for each scenario, period, and load level
for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n:
    mTEPES.vTheta[sc,p,n,mTEPES.rf.first()].fix(0)

# fixing the ESS inventory at the last load level
for sc,p,es in mTEPES.sc*mTEPES.p*mTEPES.es:
    mTEPES.vESSInventory[sc,p,mTEPES.n.last(),es].fix(pESSInitialInventory[es])

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

# fixing the ESS inventory at the end of the following pESSTimeStep (weekly, yearly), i.e., for daily ESS is fixed at the end of the week, for weekly/monthly ESS is fixed at the end of the year
for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es:
     if pStorageType[es] == 'Daily'   and mTEPES.n.ord(n) % ( 168/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Weekly'  and mTEPES.n.ord(n) % (8736/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Monthly' and mTEPES.n.ord(n) % (8736/pTimeStep) == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])

SettingUpDataTime = time.time() - StartTime
StartTime         = time.time()
print('Setting up input data                 ... ', round(SettingUpDataTime), 's')
