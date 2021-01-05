# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.7.26 - December 28, 2020

import time
import pandas as pd
import matplotlib.pyplot as plt

def OutputResults(CaseName, model):
    print('Output results              ****')

    StartTime = time.time()

    #%% outputting the generation operation
    OutputToFile = pd.Series(data=[model.vTotalGen[g]()*1e3                                  for g in model.g], index=list(model.g))
    OutputToFile.to_frame(name='MW ').reset_index().pivot_table(columns='index', values='MW ').rename_axis([None], axis=1).to_csv(model.path+'/out/oT_Result_GenerationOutput_'   +CaseName+'.csv', sep=',')

    # %% outputting the generation operation
    OutputToFile = pd.Series(data=[model.vTotalDem[d]() * 1e3                                for d in model.d], index=list(model.d))
    OutputToFile.to_frame(name='MW ').reset_index().pivot_table(columns='index', values='MW ').rename_axis([None], axis=1).to_csv(model.path + '/out/oT_Result_DemandOutput_' + CaseName + '.csv', sep=',')

    #%% outputting the network operation
    OutputToFile = pd.Series(data=[model.vFlow[ni,nf,cc]()*1e3 for ni,nf,cc in model.la], index= pd.MultiIndex.from_tuples(list(model.la)))
    OutputToFile.index.names = ['InitialNode','FinalNode','Circuit']
    OutputToFile = pd.pivot_table(OutputToFile.to_frame(name='MW'), values='MW', columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
    OutputToFile.index.names = [None] * len(OutputToFile.index.names)
    OutputToFile.to_csv(model.path + '/out/oT_Result_NetworkFlow_'+CaseName+'.csv', sep=',')

    OutputToFile = pd.Series(data=[abs(model.vFlow[ni,nf,cc]()/max(model.pLineNTCBck[ni,nf,cc],model.pLineNTCFrw[ni,nf,cc])) for ni,nf,cc in model.la], index= pd.MultiIndex.from_tuples(list(model.la)))
    OutputToFile.index.names = ['InitialNode','FinalNode','Circuit']
    OutputToFile = pd.pivot_table(OutputToFile.to_frame(name='pu'), values='pu', columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
    OutputToFile.index.names = [None] * len(OutputToFile.index.names)
    OutputToFile.to_csv(model.path + '/out/oT_Result_NetworkUtilization_'+CaseName+'.csv', sep=',')

    OutputToFile = pd.Series(data=[model.vTheta[nd]()                   for nd in model.nd], index=list(model.nd))
    OutputToFile.to_frame(name='rad').reset_index().pivot_table(columns='index', values='rad').rename_axis([None], axis=1).to_csv(model.path + '/out/oT_Result_NetworkAngle_'+CaseName+'.csv', sep=',')

    #%% outputting the LSRMC
    OutputToFile = pd.Series(data=[model.dual[model.eBalance[nd]]*1e3 for nd in model.nd], index=list(model.nd))
    OutputToFile.to_frame(name='LSRMC').reset_index().pivot_table(columns='index', values='LSRMC').rename_axis([None], axis=1).to_csv(model.path + '/out/oT_Result_LSRMC_'+CaseName+'.csv', sep=',')

    WritingResultsTime = time.time() - StartTime
    StartTime          = time.time()
    print('Writing output results                ... ', round(WritingResultsTime), 's')
