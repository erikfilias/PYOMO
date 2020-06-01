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

# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.30 - May 8, 2020
# simplicity and transparency in power systems planning

# Developed by

#    Andres Ramos, Erik Alvarez, Sara Lumbreras
#    Instituto de Investigacion Tecnologica
#    Escuela Tecnica Superior de Ingenieria - ICAI
#    UNIVERSIDAD PONTIFICIA COMILLAS
#    Alberto Aguilera 23
#    28015 Madrid, Spain
#    Andres.Ramos@comillas.edu
#    Erik.Alvarez@comillas.edu
#    https://www.iit.comillas.edu/aramos/Ramos_CV.htm

#    with the very valuable collaboration from David Dominguez (david.dominguez@comillas.edu) and Alejandro Rodriguez (argallego@comillas.edu), our local Python gurus

#%% libraries
import pandas        as pd
import time          # count clock time
import psutil        # access the number of CPUs
import pyomo.environ as pyo
from   pyomo.environ import Set, Var, Binary, NonNegativeReals, RealSet, UnitInterval, Constraint, ConcreteModel, Objective, minimize, Suffix, DataPortal, TerminationCondition
from   pyomo.opt     import SolverFactory
from   collections   import defaultdict
# import shared_stuff
import builtins

builtins.StartTime = time.time()

builtins.CaseName = 'sSEP'                              # To select the case

#%% model declaration
builtins.mTEPES = ConcreteModel('Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.30 - May 8, 2020')

#%% Input Data
import openTEPES_InputData as oT_ID

#%% Model Formulation
builtins.StartTime = time.time()

import openTEPES_ModelFormulation as oT_MF

#%% Writting LP File
# mTEPES.write('openTEPES_'+CaseName+'.lp', io_options={'symbolic_solver_labels': True})  # create lp-format file

builtins.StartTime             = time.time()
builtins.WritingLPFileTime     = time.time() - StartTime
StartTime                      = time.time()
print('Writing LP file                       ... ', round(WritingLPFileTime), 's')

#%% solving the problem
Solver = SolverFactory('gurobi')                                                       # select solver
Solver.options['LogFile'       ] = 'openTEPES_'+CaseName+'.log'
#Solver.options['IISFile'      ] = 'openTEPES_'+CaseName+'.ilp'                        # should be uncommented to show results of IIS
Solver.options['Method'        ] = 2                                                   # barrier method
Solver.options['MIPGap'        ] = 0.03
Solver.options['Threads'       ] = int((psutil.cpu_count(logical=True) + psutil.cpu_count(logical=False))/2)
#Solver.options['TimeLimit'     ] =    40000
#Solver.options['IterationLimit'] = 10000000
if oT_ID.pIndBinGenInvest*len(mTEPES.gc) + oT_ID.pIndBinNetInvest*len(mTEPES.lc) + oT_ID.pIndBinGenOperat == 0:
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual = Suffix(direction=Suffix.IMPORT)
SolverResults = Solver.solve(mTEPES, tee=True)                                         # tee=True displays the log of the solver
SolverResults.write()                                                                  # summary of the solver results

#%% fix values of binary variables to get dual variables and solve it again
if oT_ID.pIndBinGenInvest*len(mTEPES.gc) + oT_ID.pIndBinNetInvest*len(mTEPES.lc) + oT_ID.pIndBinGenOperat:
    if oT_ID.pIndBinGenInvest*len(mTEPES.gc):
        for gc in mTEPES.gc:
            mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
    if oT_ID.pIndBinNetInvest*len(mTEPES.lc):
        for lc in mTEPES.lc:
            mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
    if oT_ID.pIndBinGenOperat:
        for sc,p,n,t in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.t:
            mTEPES.vCommitment[sc,p,n,t].fix(mTEPES.vCommitment[sc,p,n,t]())
            mTEPES.vStartUp   [sc,p,n,t].fix(mTEPES.vStartUp   [sc,p,n,t]())
            mTEPES.vShutDown  [sc,p,n,t].fix(mTEPES.vShutDown  [sc,p,n,t]())
    Solver.options['relax_integrality'] = 1                                            # introduced to show results of the dual variables
    mTEPES.dual   = Suffix(direction=Suffix.IMPORT)
    SolverResults = Solver.solve(mTEPES, tee=True)                                     # tee=True displays the log of the solver
    SolverResults.write()                                                              # summary of the solver results

builtins.SolvingTime = time.time() - StartTime
builtins.StartTime = time.time()
print('Solving                               ... ', round(SolvingTime), 's')

print('Total system cost [MEUR]                  ', mTEPES.eTotalTCost.expr())

#%% Output Data

import openTEPES_OutputData as oT_OD


print('Total time                            ... ', round(oT_ID.ReadingDataTime + oT_ID.SettingUpDataTime + oT_MF.GeneratingOFTime + oT_MF.GeneratingRBITime + oT_MF.GeneratingGenConsTime + oT_MF.GeneratingRampsTime + oT_MF.GeneratingMinUDTime + oT_MF.GeneratingNetConsTime + WritingLPFileTime + SolvingTime + oT_OD.WritingResultsTime + oT_OD.PlottingNetMapsTime), 's')