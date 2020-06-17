# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

OutputResults = pd.Series(data=[mTEPES.vReactiveTotalOutput[sc,p,n,g]()*1e3                        for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults.to_frame(name='MW ').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW ').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ReactiveGenerationOutput_'+CaseName+'.csv', sep=',')

OutputResults = pd.DataFrame.from_dict(mTEPES.vP.extract_values(), orient='index', columns=[str(mTEPES.vP)])
OutputResults *= 1e3
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults, values=str(mTEPES.vP), index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_ActiveNetworkFlow_'+CaseName+'.csv', sep=',')

OutputResults = pd.DataFrame.from_dict(mTEPES.vQ.extract_values(), orient='index', columns=[str(mTEPES.vQ)])
OutputResults *= 1e3
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults, values=str(mTEPES.vQ), index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_ReactiveNetworkFlow_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc]()*mTEPES.pLineR[ni,nf,cc] for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la], index= pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la)))
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults.to_frame(name='pu'), values='pu', index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults *= 1e3
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_ActiveNetworkLosses_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vCurrentFlow_sqr[sc,p,n,ni,nf,cc]()*mTEPES.pLineX[ni,nf,cc]-mTEPES.pLineBsh[ni,nf,cc]*mTEPES.vVoltageMag_sqr[sc,p,n,ni]-mTEPES.pLineBsh[ni,nf,cc]*mTEPES.vVoltageMag_sqr[sc,p,n,nf] for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la], index= pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la)))
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults.to_frame(name='pu'), values='pu', index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults *= 1e3
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_ReactiveNetworkLosses_'+CaseName+'.csv', sep=',')

