# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

print('Lossless solve')

StartTime = time.time()

# deactivate loss constraints
mTEPES.eLineLosses1.deactivate()
mTEPES.eLineLosses2.deactivate()

pIndMemoryolving = 1
if pIndMemoryolving == 0:
    import openTEPES_ProblemSolving
else:
    import openTEPES_MemorySolving

SolvingTime = time.time() - StartTime
StartTime = time.time()
print('Solving lossless                      status ... ', round(SolvingTime), 's')

# activate loss constraints
mTEPES.eLineLosses1.activate()
mTEPES.eLineLosses2.activate()

if mTEPES.pIndBinGenInvest*len(mTEPES.gc):
    for gc in mTEPES.gc:
        mTEPES.vGenerationInvest[gc].fix(mTEPES.vGenerationInvest[gc]())
if mTEPES.pIndBinNetInvest*len(mTEPES.lc):
    for lc in mTEPES.lc:
        mTEPES.vNetworkInvest   [lc].fix(mTEPES.vNetworkInvest   [lc]())
