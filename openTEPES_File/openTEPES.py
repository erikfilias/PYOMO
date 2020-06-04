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

# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.32 - May 26, 2020
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
import time
import builtins
from   pyomo.environ import ConcreteModel

builtins.StartTime = time.time()

builtins.CaseName = '9n'                               # To select the case

#%% model declaration
builtins.mTEPES = ConcreteModel('Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.32 - May 26, 2020')

import openTEPES_InputData_AC

import openTEPES_ModelFormulation_AC
#
# # candidate discovery
# pIndCandidateDiscovery = 0
# if pIndCandidateDiscovery == 1:
#     pIndBinGenInvestSaved   = mTEPES.pIndBinGenInvest
#     pIndBinNetInvestSaved   = mTEPES.pIndBinNetInvest
#     pIndBinGenOperatSaved   = mTEPES.pIndBinGenOperat
#     mTEPES.pIndBinGenInvest = 0
#     mTEPES.pIndBinNetInvest = 0
#     mTEPES.pIndBinGenOperat = 0
#     # pIndMemorySolving = 0
#     # if pIndMemorySolving == 0:
#     #     import openTEPES_ProblemSolving
#     # else:
#     #     import openTEPES_MemorySolving
#     import openTEPES_ProblemSolving
#     import openTEPES_CandidateDiscovery
#     mTEPES.pIndBinGenInvest = pIndBinGenInvestSaved
#     mTEPES.pIndBinNetInvest = pIndBinNetInvestSaved
#     mTEPES.pIndBinGenOperat = pIndBinGenOperatSaved
#
# # solve lossless first and fix the generation and network investment decisions
# pIndLosslessSolving = 0
# if pIndLosslessSolving == 1 and mTEPES.pIndNetLosses == 1:
#     import openTEPES_LosslessSolve
#
# # solve the problem in memory (persistent) or writing the lp file
# pIndMemorySolving = 0
# if pIndMemorySolving == 0:
#     import openTEPES_ProblemSolving
# else:
#     import openTEPES_MemorySolving
#
# # stage solving with expansion decisions fixed
# pIndStageSolving        = 0
# pIndSequentialSolving   = 0
# builtins.pStageDuration = 168        # duration of the stage (weekly or monthly is what makes sense from an system operation point of view
# if pIndStageSolving == 1:
#     if pIndSequentialSolving == 1:
#         import openTEPES_SequentialStageSolve
#     else:
#         # run these commands in a python terminal, each one individually, to set up the pyro solver server
#         # os.fork('pyomo_ns')
#         # os.fork('dispatch_srvr')
#         # run this command as many times as the number of cores to allow parallel solves
#         # os.fork('pyro_mip_server')
#         import openTEPES_ParallelStageSolve
#         # pid = os.getpid()
#         # os.kill(pid)
#
# pIndOutputResults = 0
# if pIndOutputResults == 1:
#     import openTEPES_OutputResults
#
# pIndVisualization = 0
# if pIndVisualization == 1:
#     import openTEPES_Visualization
#
# TotalTime = time.time() - StartTime
# print('Total time                            ... ', round(TotalTime), 's')
