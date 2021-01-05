#                     GNU GENERAL PUBLIC LICENSE
#                        Version 3, 29 June 2007
#
#  Copyright (C) 2007 Free Software Foundation, Inc. <https://fsf.org/>
#  Everyone is permitted to copy and distribute verbatim copies
#  of this license document, but changing it is not allowed.
#
#                             Preamble
#
#   The GNU General Public License is a free, copyleft license for
# software and other kinds of works.
#
#   The licenses for most software and other practical works are designed
# to take away your freedom to share and change the works.  By contrast,
# the GNU General Public License is intended to guarantee your freedom to
# share and change all versions of a program--to make sure it remains free
# software for all its users.  We, the Free Software Foundation, use the
# GNU General Public License for most of our software; it applies also to
# any other work released this way by its authors.  You can apply it to
# your programs, too.
#
#   When we speak of free software, we are referring to freedom, not
# price.  Our General Public Licenses are designed to make sure that you
# have the freedom to distribute copies of free software (and charge for
# them if you wish), that you receive source code or can get it if you
# want it, that you can change the software or use pieces of it in new
# free programs, and that you know you can do these things.
#
#   To protect your rights, we need to prevent others from denying you
# these rights or asking you to surrender the rights.  Therefore, you have
# certain responsibilities if you distribute copies of the software, or if
# you modify it: responsibilities to respect the freedom of others.
#
#   For example, if you distribute copies of such a program, whether
# gratis or for a fee, you must pass on to the recipients the same
# freedoms that you received.  You must make sure that they, too, receive
# or can get the source code.  And you must show them these terms so they
# know their rights.
#
#   Developers that use the GNU GPL protect your rights with two steps:
# (1) assert copyright on the software, and (2) offer you this License
# giving you legal permission to copy, distribute and/or modify it.
#
#   For the developers' and authors' protection, the GPL clearly explains
# that there is no warranty for this free software.  For both users' and
# authors' sake, the GPL requires that modified versions be marked as
# changed, so that their problems will not be attributed erroneously to
# authors of previous versions.
#
#   Some devices are designed to deny users access to install or run
# modified versions of the software inside them, although the manufacturer
# can do so.  This is fundamentally incompatible with the aim of
# protecting users' freedom to change the software.  The systematic
# pattern of such abuse occurs in the area of products for individuals to
# use, which is precisely where it is most unacceptable.  Therefore, we
# have designed this version of the GPL to prohibit the practice for those
# products.  If such problems arise substantially in other domains, we
# stand ready to extend this provision to those domains in future versions
# of the GPL, as needed to protect the freedom of users.
#
#   Finally, every program is threatened constantly by software patents.
# States should not allow patents to restrict development and use of
# software on general-purpose computers, but in those that do, we wish to
# avoid the special danger that patents applied to a free program could
# make it effectively proprietary.  To prevent this, the GPL assures that
# patents cannot be used to render the program non-free.

# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.25 - April 8, 2020
# simplicity and transparency in power systems planning

# Developed by

#    Andres Ramos, Erik Alvarez
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu
#    Erik.Alvarez@comillas.edu
#    https://www.iit.comillas.edu/aramos/Ramos_CV.htm

#    with the very valuable collaboration from David Dominguez (david.dominguez@comillas.edu) and Alejandro Rodríguez (argallego@comillas.edu), our local Python gurus

#%% Libraries
import pandas        as pd
import time          # count clock time
import psutil        # access the number of CPUs
import pyomo.environ as pyo
from   pyomo.environ import Set, Var, Binary, NonNegativeReals, RealSet, Constraint, ConcreteModel, Objective, minimize, Suffix, DataPortal
from   pyomo.opt     import SolverFactory
from   collections   import defaultdict

# CANDIDATE DISCOVERY
import math
import os

StartTime = time.time()

CaseName = '9n'                              # To select the case
#%% model declaration
mTEPES = ConcreteModel('Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.25 - April 8, 2020')

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

# CANDIDATE DISCOVERY - use an ordered set for nodes
mTEPES.nd = Set(initialize=dictSets['nd'], ordered=True, doc='nodes'       )

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
dfDuration           = pd.read_csv('oT_Data_Duration_'          +CaseName+'.csv', index_col=[0    ])
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

# CANDIDATE DISCOVERY
pIndCandidateDiscovery       = dfOption   ['IndCandidateDiscovery'   ][0].astype('int')                                                                # Indicator of network losses,           0 Only user-defined candidates   - 1 discover candidates

if pIndCandidateDiscovery==1:
    pIndBinNetInvestSaved = pIndBinNetInvest
    pIndBinNetInvest = 0

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
pVariableMinStorage = dfVariableMinStorage[list(mTEPES.gg)]                                                                          # dynamic variable minimum storage    [GWh]
pVariableMaxStorage = dfVariableMaxStorage[list(mTEPES.gg)]                                                                          # dynamic variable maximum storage    [GWh]
pESSEnergyInflows   = dfESSEnergyInflows  [list(mTEPES.gg)] * 1e-3                                                                   # dynamic energy inflows              [GWh]
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
pGenToNode             = dfGeneration['Node'              ]                                                                             # generator location in node
pGenToTechnology       = dfGeneration['Technology'        ]                                                                             # generator association to technology
pRatedMaxPower         = dfGeneration['MaxPower'          ] * 1e-3                                                                      # rated maximum power                 [GW]
pRatedMinPower         = dfGeneration['MinPower'          ] * 1e-3                                                                      # rated minimum power                 [GW]
pLinearVarCost         = dfGeneration['LinearVarCost'     ] * 1e-3 * dfGeneration['FuelCost'] + dfGeneration['OMVarCost'] * 1e-3        # linear   term variable cost         [MEUR/GWh]
pConstantVarCost       = dfGeneration['ConstantVarCost'   ] * 1e-6 * dfGeneration['FuelCost']                                           # constant term variable cost         [MEUR/h]
pStartUpCost           = dfGeneration['StartUpCost'       ] * 1e-6 * dfGeneration['FuelCost']                                           # startup  cost                       [MEUR]
pShutDownCost          = dfGeneration['ShutDownCost'      ] * 1e-6 * dfGeneration['FuelCost']                                           # shutdown cost                       [MEUR]
pRampUp                = dfGeneration['RampUp'            ] * 1e-3                                                                      # ramp up   rate                      [GW/h]
pRampDw                = dfGeneration['RampDown'          ] * 1e-3                                                                      # ramp down rate                      [GW/h]
pCO2EmissionRate       = dfGeneration['CO2EmissionRate'   ] * 1e-3                                                                      # emission  rate                      [t CO2/MWh]
pUpTime                = dfGeneration['UpTime'            ]                                                                             # minimum up   time                   [h]
pDwTime                = dfGeneration['DownTime'          ]                                                                             # minimum down time                   [h]
pGenFixedCost          = dfGeneration['FixedCost'         ] *        dfGeneration['FixedChargeRate']                                    # generation fixed cost               [MEUR]
pMaxCharge             = dfGeneration['MaxCharge'         ] * 1e-3                                                                      # maximum ESS charge                  [GW]
pESSInitialInventory   = dfGeneration['InitialStorage'    ]                                                                             # initial ESS storage                 [GWh]
pESSMaxStorageCapacity = dfGeneration['MaxStorageCapacity']                                                                             # maximum ESS storage capacity        [GWh]
pESSMinStorageCapacity = dfGeneration['MinStorageCapacity']                                                                             # minimum ESS storage capacity        [GWh]
pEfficiency            = dfGeneration['Efficiency'        ]                                                                             #         ESS efficiency              [p.u.]
pStorageType           = dfGeneration['StorageType'       ]                                                                             #         ESS type

pLineType              = dfNetwork   ['LineType'          ]                                                                             # line type
pLineVoltage           = dfNetwork   ['Voltage'           ]                                                                             # line voltage                        [kV]
pLineLossFactor        = dfNetwork   ['LossFactor'        ]                                                                             # loss factor                         [p.u.]
pLineX                 = dfNetwork   ['Reactance'         ]                                                                             # reactance                           [p.u.]
pLineNTC               = dfNetwork   ['TTC'               ] * 1e-3 * dfNetwork['SecurityFactor' ]                                       # net transfer capacity               [GW]
pNetFixedCost          = dfNetwork   ['FixedCost'         ] *        dfNetwork['FixedChargeRate']                                       # network    fixed cost               [MEUR]


ReadingDataTime = time.time() - StartTime
StartTime       = time.time()
print('Reading    input data                 ... ', round(ReadingDataTime), 's')

#%% defining subsets: active load levels (n), thermal units (t), ESS units (es), candidate gen units (gc), candidate ESS units (ec), all the lines (la), candidate lines (lc) and lines with losses (ll)
mTEPES.n  = Set(initialize=mTEPES.nn,                     ordered=True , doc='load levels'        , filter=lambda mTEPES,nn      : nn        in mTEPES.nn and pDuration         [nn] >  0)
mTEPES.n2 = Set(initialize=mTEPES.nn,                     ordered=True , doc='load levels'        , filter=lambda mTEPES,nn      : nn        in mTEPES.nn and pDuration         [nn] >  0)
mTEPES.g  = Set(initialize=mTEPES.gg,                     ordered=False, doc='generating    units', filter=lambda mTEPES,gg      : gg        in mTEPES.gg and pRatedMaxPower    [gg] >  0)
mTEPES.t  = Set(initialize=mTEPES.g ,                     ordered=False, doc='thermal       units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pLinearVarCost     [g] >  0)
mTEPES.r  = Set(initialize=mTEPES.g ,                     ordered=False, doc='RES           units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pLinearVarCost     [g] == 0 and pESSMaxStorageCapacity[g] == 0)
mTEPES.es = Set(initialize=mTEPES.g ,                     ordered=False, doc='ESS           units', filter=lambda mTEPES,g       : g         in mTEPES.g  and                                 pESSMaxStorageCapacity[g] >  0)
mTEPES.gc = Set(initialize=mTEPES.g ,                     ordered=False, doc='candidate     units', filter=lambda mTEPES,g       : g         in mTEPES.g  and pGenFixedCost      [g] >  0)
mTEPES.ec = Set(initialize=mTEPES.es,                     ordered=False, doc='candidate ESS units', filter=lambda mTEPES,es      : es        in mTEPES.es and pGenFixedCost     [es] >  0)
mTEPES.la = Set(initialize=mTEPES.ni*mTEPES.nf*mTEPES.cc, ordered=False, doc='all           lines', filter=lambda mTEPES,ni,nf,cc:(ni,nf,cc) in               pLineX                     )
mTEPES.lc = Set(initialize=mTEPES.la,                     ordered=False, doc='candidate     lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost     [la] >  0)
mTEPES.cd = Set(initialize=mTEPES.la,                     ordered=False, doc='          DC  lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost     [la] >  0 and pLineType[la] == 'DC')
mTEPES.ed = Set(initialize=mTEPES.la,                     ordered=False, doc='          DC  lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pNetFixedCost     [la] == 0 and pLineType[la] == 'DC')
mTEPES.ll = Set(initialize=mTEPES.la,                     ordered=False, doc='loss          lines', filter=lambda mTEPES,*la     : la        in mTEPES.la and pLineLossFactor   [la] >  0 and pIndNetLosses >   0  )
mTEPES.rf = Set(initialize=mTEPES.nd,                     ordered=True , doc='reference node'     , filter=lambda mTEPES,nd      : nd        in               pReferenceNode             )

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

# values < 10 kW are converted to 0
pDemand          [pDemand           < 10e-6] = 0
pOperReserveUp   [pOperReserveUp    < 10e-6] = 0
pOperReserveDw   [pOperReserveDw    < 10e-6] = 0
pMinPower        [pMinPower         < 10e-6] = 0
pMaxPower        [pMaxPower         < 10e-6] = 0
pMaxPower2ndBlock[pMaxPower2ndBlock < 10e-6] = 0
pMaxCharge       [pMaxCharge        < 10e-6] = 0
pESSEnergyInflows[pESSEnergyInflows < 10e-6] = 0
pESSMinStorage   [pESSMinStorage    < 10e-6] = 0
pESSMaxStorage   [pESSMaxStorage    < 10e-6] = 0
pLineNTC         [pLineNTC          < 10e-6] = 0



# maximum flow is the net transfer capacity but to be used in the Kirchhoff's 2nd law constraint
pMaxFlow = pLineNTC

# this option avoids a warning in the following assignments
pd.options.mode.chained_assignment = None

# rounding the minimum up and down time steps to AN integer number of time steps
pUpTime = round(pUpTime/pTimeStep).astype('int')
pDwTime = round(pDwTime/pTimeStep).astype('int')

# thermal and variable units ordered by increasing variable cost
mTEPES.go = pLinearVarCost.sort_values().index

# determine the initial committed units and their output
pInitialOutput = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
pInitialUC     = pd.Series([0.]*len(mTEPES.gg), dfGeneration.index)
pSystemOutput  = 0
for go in mTEPES.go:
    n1 = next(iter(mTEPES.sc*mTEPES.p*mTEPES.n))
    if pSystemOutput < sum(pDemand[nd][n1] for nd in mTEPES.nd):
        pInitialOutput[go] = pMaxPower[go][n1]
        pInitialUC    [go] = 1
        pSystemOutput      = pSystemOutput + pInitialOutput[go]

#%% variables
mTEPES.vTotalFCost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system fixed    cost                     [MEUR]')
mTEPES.vTotalVCost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system variable cost                     [MEUR]')
mTEPES.vTotalECost           = Var(                                          within=NonNegativeReals,                                                                                         doc='total system emission cost                     [MEUR]')
mTEPES.vTotalOutput          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.g , within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,g :(0,pMaxPower          [g ][sc,p,n]),                     doc='total output of the unit                         [GW]')
mTEPES.vOutput2ndBlock       = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='second block of the unit                         [GW]')
mTEPES.vReserveUp            = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='operating reserve up                             [GW]')
mTEPES.vReserveDown          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,nr:(0,pMaxPower2ndBlock  [nr][sc,p,n]),                     doc='operating reserve down                           [GW]')
mTEPES.vESSInventory         = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,es:(pESSMinStorage[es][sc,p,n],pESSMaxStorage[es][sc,p,n]), doc='ESS inventory                                   [GWh]')
mTEPES.vESSSpillage          = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, within=NonNegativeReals,                                                                                         doc='ESS spillage                                    [GWh]')
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

mTEPES.vLineLosses           = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, within=NonNegativeReals, bounds=lambda mTEPES,sc,p,n,*ll:(0,0.5*pLineLossFactor[ll]*pLineNTC[ll]),               doc='half line losses                                 [GW]')
mTEPES.vFlow                 = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.la, within=RealSet,          bounds=lambda mTEPES,sc,p,n,*la:(-pLineNTC[la],pLineNTC[la]),                           doc='flow                                             [GW]')
mTEPES.vTheta                = Var(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, within=RealSet,                                                                                                  doc='voltage angle                                   [rad]')

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
     if pStorageType[es] == 'Daily'   and mTEPES.n.ord(n) %  168/pTimeStep == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Weekly'  and mTEPES.n.ord(n) % 8736/pTimeStep == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])
     if pStorageType[es] == 'Monthly' and mTEPES.n.ord(n) % 8736/pTimeStep == 0:
         mTEPES.vESSInventory[sc,p,n,es].fix(pESSInitialInventory[es])

SettingUpDataTime = time.time() - StartTime
StartTime         = time.time()
print('Setting up input data                 ... ', round(SettingUpDataTime), 's')

def eTotalFCost(mTEPES):
   return mTEPES.vTotalFCost == sum(pNetFixedCost[lc] * mTEPES.vNetworkInvest   [lc] for lc in mTEPES.lc) + sum(pGenFixedCost[gc] * mTEPES.vGenerationInvest[gc] for gc in mTEPES.gc)
mTEPES.eTotalFCost = Constraint(rule=eTotalFCost, doc='total system fixed    cost [MEUR]')

def eTotalVCost(mTEPES):
    return mTEPES.vTotalVCost == (sum(pScenProb[sc] * pDuration[n] * pENSCost             * mTEPES.vENS        [sc,p,n,nd] for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd) +
                                  sum(pScenProb[sc] * pDuration[n] * pLinearVarCost  [nr] * mTEPES.vTotalOutput[sc,p,n,nr]                                                         +
                                      pScenProb[sc] * pDuration[n] * pConstantVarCost[nr] * mTEPES.vCommitment [sc,p,n,nr]                                                         +
                                      pScenProb[sc]                * pStartUpCost    [nr] * mTEPES.vStartUp    [sc,p,n,nr]                                                         +
                                      pScenProb[sc]                * pShutDownCost   [nr] * mTEPES.vShutDown   [sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr) )
mTEPES.eTotalVCost = Constraint(rule=eTotalVCost, doc='total system variable cost [MEUR]')

def eTotalECost(mTEPES):
    return mTEPES.vTotalECost == sum(pScenProb[sc] * pCO2Cost * pCO2EmissionRate[nr] * mTEPES.vTotalOutput[sc,p,n,nr] for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)
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
    return mTEPES.vTotalOutput[sc,p,n,ec] / pMaxPower[sc,p,n,ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalGenESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalGenESS, doc='output      if installed ESS unit [p.u.]')

print('eInstalGenESS         ... ', len(mTEPES.eInstalGenESS), ' rows')

def eInstalConESS(mTEPES,sc,p,n,ec):
    return mTEPES.vTotalOutput[sc,p,n,ec] / pMaxCharge      [ec] <= mTEPES.vGenerationInvest[ec]
mTEPES.eInstalConESS = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ec, rule=eInstalConESS, doc='consumption if installed ESS unit [p.u.]')

print('eInstalConESS         ... ', len(mTEPES.eInstalConESS), ' rows')

#%%
def eOperReserveUp(mTEPES,sc,p,n):
    if pOperReserveUp[sc,p,n]:
        return sum(mTEPES.vReserveUp  [sc,p,n,nr] for nr in mTEPES.nr) >= pOperReserveUp[sc,p,n]
    else:
        return Constraint.Skip
mTEPES.eOperReserveUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveUp, doc='up   operating reserve [GW]')

print('eOperReserveUp        ... ', len(mTEPES.eOperReserveUp), ' rows')

def eOperReserveDw(mTEPES,sc,p,n):
    if pOperReserveDw[sc,p,n]:
        return sum(mTEPES.vReserveDown[sc,p,n,nr] for nr in mTEPES.nr) >= pOperReserveDw[sc,p,n]
    else:
        return Constraint.Skip
mTEPES.eOperReserveDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, rule=eOperReserveDw, doc='down operating reserve [GW]')

print('eOperReserveDw        ... ', len(mTEPES.eOperReserveDw), ' rows')

def eBalance(mTEPES,sc,p,n,nd):
    return (sum(mTEPES.vTotalOutput[sc,p,n,g] for g in pNode2Gen[nd]) - sum(mTEPES.vESSCharge[sc,p,n,es] for es in pNode2ESS[nd]) + mTEPES.vENS[sc,p,n,nd] == pDemand[nd][sc,p,n] +
        sum(mTEPES.vLineLosses[sc,p,n,nd,lout ] for lout  in loutl[nd]) + sum(mTEPES.vFlow[sc,p,n,nd,lout ] for lout  in lout[nd]) +
        sum(mTEPES.vLineLosses[sc,p,n,ni,nd,cc] for ni,cc in linl [nd]) - sum(mTEPES.vFlow[sc,p,n,ni,nd,cc] for ni,cc in lin [nd]))
mTEPES.eBalance = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nd, rule=eBalance, doc='load generation balance [GW]')

print('eBalance              ... ', len(mTEPES.eBalance), ' rows')

def eESSInventory(mTEPES,sc,p,n,es):
    if   mTEPES.n.ord(n) == pESSTimeStep[es]:
        return pESSInitialInventory[es]                                        + sum(pDuration[n2]*(pESSEnergyInflows[es][sc,p,n2] - mTEPES.vTotalOutput[sc,p,n2,es] + pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    elif mTEPES.n.ord(n) >  pESSTimeStep[es] and mTEPES.n.ord(n) % pESSTimeStep[es] == 0:
        return mTEPES.vESSInventory[sc,p,mTEPES.n.prev(n,pESSTimeStep[es]),es] + sum(pDuration[n2]*(pESSEnergyInflows[es][sc,p,n2] - mTEPES.vTotalOutput[sc,p,n2,es] + pEfficiency[es]*mTEPES.vESSCharge[sc,p,n2,es]) for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-pESSTimeStep[es]:mTEPES.n.ord(n)]) == mTEPES.vESSInventory[sc,p,n,es] + mTEPES.vESSSpillage[sc,p,n,es]
    else:
        return Constraint.Skip
mTEPES.eESSInventory = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.es, rule=eESSInventory, doc='ESS inventory balance [GWh]')

print('eESSInventory         ... ', len(mTEPES.eESSInventory), ' rows')

GeneratingRBITime = time.time() - StartTime
StartTime         = time.time()
print('Generating reserves/balance/inventory ... ', round(GeneratingRBITime), 's')

#%%
def eMaxOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n] and n != mTEPES.n.last():
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] + mTEPES.vReserveUp  [sc,p,n,nr]) / pMaxPower2ndBlock[nr][sc,p,n] <= mTEPES.vCommitment[sc,p,n,nr] - mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,mTEPES.n.next(n),nr]
    elif pOperReserveUp[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n] and n == mTEPES.n.last():
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] + mTEPES.vReserveUp  [sc,p,n,nr]) / pMaxPower2ndBlock[nr][sc,p,n] <= mTEPES.vCommitment[sc,p,n,nr] - mTEPES.vStartUp[sc,p,n,nr]
    else:
        return Constraint.Skip
mTEPES.eMaxOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMaxOutput2ndBlock, doc='max output of the second block of a committed unit [p.u.]')

print('eMaxOutput2ndBlock    ... ', len(mTEPES.eMaxOutput2ndBlock), ' rows')

def eMinOutput2ndBlock(mTEPES,sc,p,n,nr):
    if   pOperReserveDw[sc,p,n] and pMaxPower2ndBlock[nr][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,nr] - mTEPES.vReserveDown[sc,p,n,nr]) / pMaxPower2ndBlock[nr][sc,p,n] >= 0
    else:
        return Constraint.Skip
mTEPES.eMinOutput2ndBlock = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eMinOutput2ndBlock, doc='min output of the second block of a committed unit [p.u.]')

print('eMinOutput2ndBlock    ... ', len(mTEPES.eMinOutput2ndBlock), ' rows')

def eTotalOutput(mTEPES,sc,p,n,nr):
    return mTEPES.vTotalOutput[sc,p,n,nr] == pMinPower[nr][sc,p,n] * mTEPES.vCommitment[sc,p,n,nr] + mTEPES.vOutput2ndBlock[sc,p,n,nr]
mTEPES.eTotalOutput = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eTotalOutput, doc='total output of a unit [GW]')

print('eTotalOutput          ... ', len(mTEPES.eTotalOutput), ' rows')

def eUCStrShut(mTEPES,sc,p,n,nr):
    if mTEPES.n.first():
        return mTEPES.vCommitment[sc,p,n,nr] - pInitialUC[nr]                               == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
    else:
        return mTEPES.vCommitment[sc,p,n,nr] - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),nr] == mTEPES.vStartUp[sc,p,n,nr] - mTEPES.vShutDown[sc,p,n,nr]
mTEPES.eUCStrShut = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.nr, rule=eUCStrShut, doc='relation among commitment startup and shutdown')

print('eUCStrShut            ... ', len(mTEPES.eUCStrShut), ' rows')

GeneratingGenConsTime = time.time() - StartTime
StartTime             = time.time()
print('Generating generation constraints     ... ', round(GeneratingGenConsTime), 's')

#%%
def eRampUp(mTEPES,sc,p,n,t):
    if   pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n] and mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(pInitialOutput[t]-pMinPower[t][sc,p,n],0)   + mTEPES.vReserveUp  [sc,p,n,t]) / pDuration[n] / pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    elif pRampUp[t] and pRampUp[t] < pMaxPower2ndBlock[t][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t] + mTEPES.vReserveUp  [sc,p,n,t]) / pDuration[n] / pRampUp[t] <=   mTEPES.vCommitment[sc,p,n,t] - mTEPES.vStartUp[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampUp = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampUp, doc='maximum ramp up   [p.u.]')

print('eRampUp               ... ', len(mTEPES.eRampUp), ' rows')

def eRampDw(mTEPES,sc,p,n,t):
    if   pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n] and mTEPES.n.first():
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - max(pInitialOutput[t]-pMinPower[t][sc,p,n],0)   - mTEPES.vReserveDown[sc,p,n,t]) / pDuration[n] / pRampDw[t] >= - pInitialUC[t]                               + mTEPES.vShutDown[sc,p,n,t]
    elif pRampDw[t] and pRampDw[t] < pMaxPower2ndBlock[t][sc,p,n]:
        return (mTEPES.vOutput2ndBlock[sc,p,n,t] - mTEPES.vOutput2ndBlock[sc,p,mTEPES.n.prev(n),t] - mTEPES.vReserveDown[sc,p,n,t]) / pDuration[n] / pRampDw[t] >= - mTEPES.vCommitment[sc,p,mTEPES.n.prev(n),t] + mTEPES.vShutDown[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eRampDw = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eRampDw, doc='maximum ramp down [p.u.]')

print('eRampDw               ... ', len(mTEPES.eRampDw), ' rows')

GeneratingRampsTime = time.time() - StartTime
StartTime           = time.time()
print('Generating ramps   up/down            ... ', round(GeneratingRampsTime), 's')

#%%
def eMinUpTime(mTEPES,sc,p,n,t):
    if pUpTime[t] > 1 and mTEPES.n.ord(n) >= pUpTime[t]:
        return sum(mTEPES.vStartUp [sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-pUpTime[t]:mTEPES.n.ord(n)]) <=     mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinUpTime   = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinUpTime  , doc='minimum up   time [h]')

print('eMinUpTime            ... ', len(mTEPES.eMinUpTime), ' rows')

def eMinDownTime(mTEPES,sc,p,n,t):
    if pDwTime[t] > 1 and mTEPES.n.ord(n) >= pDwTime[t]:
        return sum(mTEPES.vShutDown[sc,p,n2,t] for n2 in list(mTEPES.n2)[mTEPES.n.ord(n)-pDwTime[t]:mTEPES.n.ord(n)]) <= 1 - mTEPES.vCommitment[sc,p,n,t]
    else:
        return Constraint.Skip
mTEPES.eMinDownTime = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.t, rule=eMinDownTime, doc='minimum down time [h]')

print('eMinDownTime          ... ', len(mTEPES.eMinDownTime), ' rows')

GeneratingMinUDTime = time.time() - StartTime
StartTime           = time.time()
print('Generating minimum up/down time       ... ', round(GeneratingMinUDTime), 's')

#%%
def eInstalNetCap1(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / pLineNTC[ni,nf,cc] >= - mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eInstalNetCap1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap1, doc='maximum flow by installed network capacity [p.u.]')

print('eInstalNetCap1        ... ', len(mTEPES.eInstalNetCap1), ' rows')

def eInstalNetCap2(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / pLineNTC[ni,nf,cc] <=   mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eInstalNetCap2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lc, rule=eInstalNetCap2, doc='maximum flow by installed network capacity [p.u.]')

print('eInstalNetCap2        ... ', len(mTEPES.eInstalNetCap2), ' rows')

def eKirchhoff2ndLawExst(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / pMaxFlow[ni,nf,cc] == (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / pLineX[ni,nf,cc] / pMaxFlow[ni,nf,cc] * pSBase
mTEPES.eKirchhoff2ndLawExst = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lea, rule=eKirchhoff2ndLawExst, doc='flow for each AC existing  line [rad]')

print('eKirchhoff2ndLawExst  ... ', len(mTEPES.eKirchhoff2ndLawExst), ' rows')

def eKirchhoff2ndLawCnd1(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / pLineX[ni,nf,cc] / pMaxFlow[ni,nf,cc] * pSBase >= - 1 + mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eKirchhoff2ndLawCnd1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd1, doc='flow for each AC candidate line [rad]')

print('eKirchhoff2ndLawCnd1  ... ', len(mTEPES.eKirchhoff2ndLawCnd1), ' rows')

def eKirchhoff2ndLawCnd2(mTEPES,sc,p,n,ni,nf,cc):
    return mTEPES.vFlow[sc,p,n,ni,nf,cc] / pMaxFlow[ni,nf,cc] -  (mTEPES.vTheta[sc,p,n,ni] - mTEPES.vTheta[sc,p,n,nf]) / pLineX[ni,nf,cc] / pMaxFlow[ni,nf,cc] * pSBase <=   1 - mTEPES.vNetworkInvest[ni,nf,cc]
mTEPES.eKirchhoff2ndLawCnd2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.lca, rule=eKirchhoff2ndLawCnd2, doc='flow for each AC candidate line [rad]')

print('eKirchhoff2ndLawCnd2  ... ', len(mTEPES.eKirchhoff2ndLawCnd2), ' rows')

def eLineLosses1(mTEPES,sc,p,n,ni,nf,cc):
    if pIndNetLosses:
        return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >= - 0.5 * pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
mTEPES.eLineLosses1 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses1, doc='ohmic losses for all the lines [GW]')

print('eLineLosses1          ... ', len(mTEPES.eLineLosses1), ' rows')

def eLineLosses2(mTEPES,sc,p,n,ni,nf,cc):
    if pIndNetLosses:
        return mTEPES.vLineLosses[sc,p,n,ni,nf,cc] >=   0.5 * pLineLossFactor[ni,nf,cc] * mTEPES.vFlow[sc,p,n,ni,nf,cc]
mTEPES.eLineLosses2 = Constraint(mTEPES.sc, mTEPES.p, mTEPES.n, mTEPES.ll, rule=eLineLosses2, doc='ohmic losses for all the lines [GW]')

print('eLineLosses2          ... ', len(mTEPES.eLineLosses2), ' rows')

GeneratingNetConsTime = time.time() - StartTime
StartTime             = time.time()
print('Generating network    constraints     ... ', round(GeneratingNetConsTime), 's')

#%% solving the problem
mTEPES.write('openTEPES_'+CaseName+'.lp', io_options={'symbolic_solver_labels':True})  # create lp-format file
Solver = SolverFactory('gurobi')                                                       # select solver
Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
#Solver.options['ResultFile'   ] = 'openTEPES_'+CaseName+'.ilp'                        # introduced to show results of IIS
#Solver.options['Method'       ] = 2                                                   # barrier method
Solver.options['MIPGap'        ] = 0.03
Solver.options['Threads'       ] = (psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2
#Solver.options['TimeLimit'     ] =    40000
#Solver.options['IterationLimit'] = 10000000
if pIndBinGenInvest + pIndBinNetInvest + pIndBinGenOperat == 0:
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual = Suffix(direction=Suffix.IMPORT)
SolverResults = Solver.solve(mTEPES, tee=True)                                         # tee=True displays the output of the solver
SolverResults.write()     
# summary of the solver results

# CANDIDATE DISCOVERY
# The way to use this is to run it and wait for the file UpgradedNetwork to be generated
# Then, upgraded network is fed to the model in a subsequent iteration
# Transformers, conversors could be added if we add node voltage
# We could introduce the more sophisticated version based on angles, but it seems that pMaxTheta is not defined yet
# Bound on candidates is established by hand, we could also introduce it in the options file
# In this version a linear relaxation of the problem is solved first. Normally, we would restrict investment at all, but I could not find a parameter for this. Is there a reason why it is not defined? 


pBoundOnCandidates = 50

if pIndCandidateDiscovery==1:    

# Initializing 
    
    dfCables            = pd.read_csv('oT_Data_Cables_'           +CaseName+'.csv', index_col=[0])
    dfCables.fillna(0, inplace=True)

    pLinearDistance = pd.DataFrame(0, index = range(len(nd)*len(nd)), columns = ['ni', 'nf', 'LinearDistance'])
    pPotentialBenefit = pd.DataFrame(0, index = range(len(nd)*len(nd)), columns = ['ni', 'nf', 'PotentialBenefit'])
    pPotentialBenefitRatio = pd.DataFrame(0, index = range(len(nd)*len(nd)), columns = ['ni', 'nf', 'cc', 'PotentialBenefitRatio'])   
    pPotentialBenefitRatioSorted = pd.DataFrame(0, index = range(len(nd)*len(nd)), columns = ['ni', 'nf', 'cc', 'PotentialBenefitRatio'])   

    AuxIndex  = 0
    AuxIndex2 = 0

# For every pair of nodes, calculating potential benefit and cost
    for ni,nf in mTEPES.nd*mTEPES.nd:
        if mTEPES.nd.ord(ni) < mTEPES.nd.ord(nf):
            ThisDistance = 2 * math.asin(math.sqrt(math.sin((pNodeLat[nf]-pNodeLat[ni])*math.pi/180/2)**2 + 
                                    math.cos(pNodeLat[ni]*math.pi/180)*math.cos(pNodeLat[nf]*math.pi/180)*(math.sin((pNodeLon[nf]-
                                    pNodeLon[ni])*math.pi/180/2)**2)))
            pLinearDistance.loc[AuxIndex] = [ni, nf, ThisDistance]
            
# Potential Benefit per unit capacity
            ThisBenefit = sum(abs(mTEPES.dual[mTEPES.eBalance[sc,p,n,nf]]-mTEPES.dual[mTEPES.eBalance[sc,p,n,ni]])*1e3/pScenProb[sc]/pDuration[n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n) 
            pPotentialBenefit.loc[AuxIndex] = [ni, nf, ThisBenefit]                      
            
            for cabletype in mTEPES.cc:
                
                if ThisDistance <= dfCables.loc[cabletype,'MaximumDistanceCT'].item() and ThisDistance >= dfCables.loc[cabletype,'MinimumDistanceCT'].item():
                    
# Cost of candidate
                    ThisCost = dfCables.loc[cabletype,'FixedCostperKmCT'].item()*ThisDistance
# Potential Benefit ratio using the capacity
                    ThisCapacity = dfCables.loc[cabletype,'TTCCT'].item()
                    if ThisCost > 0:
                        pPotentialBenefitRatio.loc[AuxIndex2] = [ni, nf, cabletype, ThisBenefit * ThisCapacity/(ThisCost+0.0000001)] 
                        AuxIndex2 = AuxIndex2 + 1
                    
            AuxIndex = AuxIndex +1
            
    pLinearDistance.set_index(['ni','nf'])
    pPotentialBenefit.set_index(['ni','nf'])
    pPotentialBenefitRatio.set_index(['ni','nf','cc'])   
       
    pNetworkNew = pd.DataFrame(0, range(pBoundOnCandidates), columns = [ 'ni', 'nf', 'cc', 'LineType','Voltage','LossFactor','Reactance','TTC','SecurityFactor','FixedCost','FixedChargeRate'])
    
# Sorting and selecting candidates
    pPBR  = pd.read_csv('PotentialBenefitRatio'+'.csv', index_col=[0, 1 ,2, 3])
    pPotentialBenefitRatioSorted = pPBR.sort_values(by = ['PotentialBenefitRatio'], ascending = False)
    pSorted = pPotentialBenefitRatioSorted[0:pBoundOnCandidates]

# hasta aquí ok
    listni = pSorted.index.get_level_values('ni')
    listnf = pSorted.index.get_level_values('nf')
    listnc = pSorted.index.get_level_values('cc')

# Adding the information to the lines - candidate lines
    pNumberLines = len(dfNetwork.index)
    for i in range(len(pSorted)):    
    
        ni = listni[i]
        nf = listnf[i]
        cc = listnc[i]
        
        ThisDistance = pLinearDistance.loc[(pLinearDistance.ni==ni)].loc[(pLinearDistance.nf==nf), 'LinearDistance'].item()
        ThisLineType = dfCables.loc[cc,'LineType']
        ThisVoltage = dfCables.loc[cc,'Voltage']
        ThisLossFactor = dfCables.loc[cc,'LossFactor'] 
        ThisReactance = dfCables.loc[cc,'XperKmCT'] * ThisDistance
        ThisTTC = dfCables.loc[cc,'TTCCT']
        ThisSecurityFactor = dfCables.loc[cc,'SecurityFactor']
        ThisFixedCost = dfCables.loc[cc,'FixedCostperKmCT'] * ThisDistance
        ThisFixedChargeRate= dfCables.loc[cc,'FixedChargeRate']
        pNetworkNew.loc[i-1] = [ni,nf,cc,ThisLineType,ThisVoltage,ThisLossFactor,ThisReactance,ThisTTC,ThisSecurityFactor,ThisFixedCost,ThisFixedChargeRate]
        dfNetwork.loc[(ni,nf,cc)] = [ThisLineType,ThisVoltage,ThisLossFactor,ThisReactance,ThisTTC,ThisSecurityFactor,ThisFixedCost,ThisFixedChargeRate]
    
    dfNetwork.to_csv('UpgradedNetwork_' + CaseName + '.csv')
#%% fix values of binary variables to get dual variables and solve it again
if pIndBinGenInvest*len(mTEPES.gc) + pIndBinNetInvest*len(mTEPES.lc) + pIndBinGenOperat:
    if pIndBinGenInvest:
        for gc in mTEPES.gc:
            mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
    if pIndBinNetInvest:
        for lc in mTEPES.lc:
            mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
    if pIndBinGenOperat:
        for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
            mTEPES.vCommitment[sc,p,n,t].fix(mTEPES.vCommitment[sc,p,n,t]())
            mTEPES.vStartUp   [sc,p,n,t].fix(mTEPES.vStartUp   [sc,p,n,t]())
            mTEPES.vShutDown  [sc,p,n,t].fix(mTEPES.vShutDown  [sc,p,n,t]())
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual   = Suffix(direction=Suffix.IMPORT)
    SolverResults = Solver.solve(mTEPES, tee=True)                                     # tee=True displays the output of the solver
    SolverResults.write()                                                              # summary of the solver results

SolvingTime = time.time() - StartTime
StartTime   = time.time()
print('Solving                               ... ', round(SolvingTime), 's')

print('Objective function value                  ', mTEPES.eTotalTCost.expr())

#%% inverse index generator to technology
pTechnologyToGen = pGenToTechnology.reset_index().set_index('Technology').set_axis(['Generator'], axis=1, inplace=False)['Generator']

import matplotlib.pyplot as plt

#%% outputting the investment decisions
if len(mTEPES.gc):
    OutputResults = pd.DataFrame.from_dict(mTEPES.vGenerationInvest.extract_values(), orient='index', columns=[str(mTEPES.vGenerationInvest)])
    OutputResults.index.names = ['Generator']
    OutputResults.to_csv('oT_Result_GenerationInvestment_'+CaseName+'.csv', index=True, header=True)
if len(mTEPES.lc):
    OutputResults = pd.DataFrame.from_dict(mTEPES.vNetworkInvest.extract_values(), orient='index', columns=[str(mTEPES.vNetworkInvest)])
    OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
    OutputResults.index.names = ['InitialNode','FinalNode','Circuit']
    pd.pivot_table(OutputResults, values=str(mTEPES.vNetworkInvest), index=['InitialNode','FinalNode'], columns=['Circuit'], fill_value=0).to_csv('oT_Result_NetworkInvestment_'+CaseName+'.csv', sep=',')

#%% outputting the generation operation
OutputResults = pd.Series(data=[mTEPES.vCommitment[sc,p,n,nr]() for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t)))
OutputResults.to_frame(name='p.u.').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='p.u.').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationCommitment_'+CaseName+'.csv', sep=',')
OutputResults = pd.Series(data=[mTEPES.vStartUp   [sc,p,n,nr]() for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t)))
OutputResults.to_frame(name='p.u.').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='p.u.').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationStartUp_'   +CaseName+'.csv', sep=',')
OutputResults = pd.Series(data=[mTEPES.vShutDown  [sc,p,n,nr]() for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t)))
OutputResults.to_frame(name='p.u.').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='p.u.').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationShutDown_'  +CaseName+'.csv', sep=',')

if sum(pOperReserveUp[sc,p,n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n):
    OutputResults = pd.Series(data=[mTEPES.vReserveUp  [sc,p,n,nr]()*1e3 for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationReserveUp_'  +CaseName+'.csv', sep=',')

if sum(pOperReserveDw[sc,p,n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n):
    OutputResults = pd.Series(data=[mTEPES.vReserveDown[sc,p,n,nr]()*1e3 for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationReserveDown_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*1e3                        for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults.to_frame(name='MW ').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW ').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationOutput_'+CaseName+'.csv', sep=',')

if len(mTEPES.r):
    OutputResults = pd.Series(data=[(pMaxPower[r][sc,p,n]-mTEPES.vTotalOutput[sc,p,n,r]())*1e3 for sc,p,n,r in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.r], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.r)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_RESCurtailment_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*pDuration[n]                   for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationEnergy_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,nr]()*pCO2EmissionRate[nr]*1e3      for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t)))
OutputResults.to_frame(name='tCO2').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='tCO2').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationEmission_'+CaseName+'.csv', sep=',')

#%% outputting the ESS operation
if len(mTEPES.es):
    OutputResults = pd.Series(data=[mTEPES.vESSCharge   [sc,p,n,es]()*1e3              for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSChargeOutput_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[mTEPES.vESSCharge   [sc,p,n,es]()*pDuration[n]     for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSChargeEnergy_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(pTechnologyToGen[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSTechnologyEnergy_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[mTEPES.vESSInventory[sc,p,n,es]()                  for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh', dropna=False).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSInventory_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[mTEPES.vESSSpillage [sc,p,n,es]()                  for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh', dropna=False).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSSpillage_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series({Key:OptimalSolution.value*1e3 for Key,OptimalSolution in mTEPES.vESSCharge.items()})
    OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(pTechnologyToGen[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyCharge_'+CaseName+'.csv', sep=',')

    TechnologyCharge = OutputResults.loc[:,:,:,:]

OutputResults = pd.Series({Key:OptimalSolution.value*1e3 for Key,OptimalSolution in mTEPES.vTotalOutput.items()})
OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(pTechnologyToGen[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyOutput_'+CaseName+'.csv', sep=',')

TechnologyOutput = OutputResults.loc[:,:,:,:]

for sc,p in mTEPES.sc*mTEPES.p:
    fig, fg = plt.subplots()
    fg.stackplot(range(len(mTEPES.n)),  TechnologyOutput.loc[sc,p,:,:].values.reshape(len(mTEPES.n),len(mTEPES.gt)).transpose().tolist(), labels=list(mTEPES.gt))
    #fg.stackplot(range(len(mTEPES.n)), -TechnologyCharge.loc[sc,p,:,:].values.reshape(len(mTEPES.n),len(mTEPES.gt)).transpose().tolist(), labels=list(mTEPES.gt))
    fg.plot(range(len(mTEPES.n)), pDemand.sum(axis=1)*1e3, label='Demand', linewidth=1, color='k')
    fg.set(xlabel='Hours', ylabel='MW')
    plt.title(sc)
    fg.tick_params(axis='x', rotation=90)
    fg.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('oT_TechnologyOutput_'+sc+'_'+p+'_'+CaseName+'.png', bbox_inches='tight')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*pDuration[n] for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(pTechnologyToGen[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyEnergy_'+CaseName+'.csv', sep=',')

#%% outputting the network operation
OutputResults = pd.DataFrame.from_dict(mTEPES.vFlow.extract_values(), orient='index', columns=[str(mTEPES.vFlow)])
OutputResults *= 1e3
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults, values=str(mTEPES.vFlow), index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_NetworkFlow_' + CaseName + '.csv', sep=',')

OutputResults = pd.Series(data=[abs(mTEPES.vFlow[sc,p,n,ni,nf,cc]()/pLineNTC[ni,nf,cc]) for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la], index= pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la)))
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults.to_frame(name='pu'), values='pu', index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_NetworkUtilization_'+CaseName+'.csv', sep=',')

OutputResults = pd.DataFrame.from_dict(mTEPES.vTheta.extract_values(), orient='index', columns=[str(mTEPES.vTheta)])
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','Node']
pd.pivot_table(OutputResults, values=str(mTEPES.vTheta), index=['Scenario','Period','LoadLevel'], columns=['Node'], fill_value=0).to_csv('oT_Result_NetworkAngle_'+CaseName+'.csv', sep=',')

OutputResults = pd.DataFrame.from_dict(mTEPES.vENS.extract_values(), orient='index', columns=[str(mTEPES.vENS)])
OutputResults *= 1e3
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','Node']
pd.pivot_table(OutputResults, values=str(mTEPES.vENS), index=['Scenario','Period','LoadLevel'], columns=['Node'], fill_value=0).to_csv('oT_Result_NetworkENS_'+CaseName+'.csv', sep=',')

#%% outputting the LSRMC
OutputResults = pd.Series(data=[mTEPES.dual[mTEPES.eBalance[sc,p,n,nd]]*1e3/pScenProb[sc]/pDuration[n] for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd)))
OutputResults.to_frame(name='LSRMC').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='LSRMC').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_LSRMC_'+CaseName+'.csv', sep=',')

LSRMC = OutputResults.loc[:,:]

fig, fg = plt.subplots()
for nd in mTEPES.nd:
    fg.plot(range(len(mTEPES.sc*mTEPES.p*mTEPES.n)), LSRMC[:,:,:,nd], label=nd)
fg.set(xlabel='Hours', ylabel='EUR/MWh')
fg.set_ybound(lower=0, upper=100)
plt.title('LSRMC')
fg.tick_params(axis='x', rotation=90)
fg.legend()
plt.tight_layout()
plt.show()
plt.savefig('oT_LSRMC_'+CaseName+'.png', bbox_inches='tight')

WritingResultsTime = time.time() - StartTime
StartTime          = time.time()
print('Writing output results                ... ', round(WritingResultsTime), 's')

# CANDIDATE DISCOVERY
# I get rid of geopandas as it is not working right now...

    
    
