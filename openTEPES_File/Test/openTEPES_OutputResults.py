# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

import time
import pandas as pd
import matplotlib.pyplot as plt

print('Output results              ****')

StartTime = time.time()

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

if sum(mTEPES.pOperReserveUp[sc,p,n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n):
    OutputResults = pd.Series(data=[mTEPES.vReserveUp  [sc,p,n,nr]()*1e3 for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationReserveUp_'  +CaseName+'.csv', sep=',')

if sum(mTEPES.pOperReserveDw[sc,p,n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n):
    OutputResults = pd.Series(data=[mTEPES.vReserveDown[sc,p,n,nr]()*1e3 for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationReserveDown_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*1e3                        for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults.to_frame(name='MW ').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW ').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationOutput_'+CaseName+'.csv', sep=',')

if len(mTEPES.r):
    OutputResults = pd.Series(data=[(mTEPES.pMaxPower[sc,p,n,r]-mTEPES.vTotalOutput[sc,p,n,r]())*1e3 for sc,p,n,r  in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.r ], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.r )))
    OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_RESCurtailment_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*mTEPES.pDuration[n]                  for sc,p,n,g  in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g ], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g )))
OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationEnergy_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,nr]()*mTEPES.pCO2EmissionRate[nr]*1e3     for sc,p,n,nr in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nr)))
OutputResults.to_frame(name='tCO2').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='tCO2').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_GenerationEmission_'+CaseName+'.csv', sep=',')

#%% outputting the ESS operation
if len(mTEPES.es):
    OutputResults = pd.Series(data=[mTEPES.vESSCharge   [sc,p,n,es]()*1e3                           for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='MW' ).reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW' ).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSChargeOutput_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[mTEPES.vESSCharge   [sc,p,n,es]()*mTEPES.pDuration[n]           for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSChargeEnergy_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(mTEPES.t2g[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSTechnologyEnergy_'+CaseName+'.csv', sep=',')

    ESSTechnologyEnergy = pd.Series(data=[0]*len(mTEPES.gt), index=list(mTEPES.gt))
    for gt in mTEPES.gt:
        ESSTechnologyEnergy[gt] = OutputResults.loc[:,:,:,gt].sum()

    fig, fg = plt.subplots()
    def func(percentage, ESSTechnologyEnergy):
        absolute = int(percentage/100.*ESSTechnologyEnergy.sum())
        return '{:.1f}%\n({:d} GWh)'.format(percentage, absolute)
    wedges, texts, autotexts = plt.pie(ESSTechnologyEnergy, labels=list(mTEPES.gt), autopct=lambda percentage: func(percentage, ESSTechnologyEnergy), shadow=True, startangle=90, textprops=dict(color='w'))
    plt.legend(wedges, list(mTEPES.gt), title='Energy Consumption', loc='center left', bbox_to_anchor=(1, 0, 0.5, 1))
    plt.setp(autotexts, size=8, weight='bold')
    plt.axis('equal')
    plt.show()
    plt.savefig('oT_ESSTechnologyEnergy_'+CaseName+'.png', bbox_inches='tight')

    OutputResults = pd.Series(data=[mTEPES.vESSInventory[sc,p,n,es]()              for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults *= 1e3
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh', dropna=False).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSInventory_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series(data=[mTEPES.vESSSpillage [sc,p,n,es]()              for sc,p,n,es in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.es)))
    OutputResults *= 1e3
    OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh', dropna=False).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_ESSSpillage_'+CaseName+'.csv', sep=',')

    OutputResults = pd.Series({Key:OptimalSolution.value*1e3 for Key,OptimalSolution in mTEPES.vESSCharge.items()})
    OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(mTEPES.t2g[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
    OutputResults.to_frame(name='MW' ).reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW' ).rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyCharge_'+CaseName+'.csv', sep=',')

    TechnologyCharge = OutputResults.loc[:,:,:,:]

OutputResults = pd.Series({Key:OptimalSolution.value*1e3 for Key,OptimalSolution in mTEPES.vTotalOutput.items()})
OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(mTEPES.t2g[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
OutputResults.to_frame(name='MW').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='MW').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyOutput_'+CaseName+'.csv', sep=',')

TechnologyOutput = OutputResults.loc[:,:,:,:]

for sc,p in mTEPES.sc*mTEPES.p:
    fig, fg = plt.subplots()
    fg.stackplot(range(len(mTEPES.n)),  TechnologyOutput.loc[sc,p,:,:].values.reshape(len(mTEPES.n),len(mTEPES.gt)).transpose().tolist(), labels=list(mTEPES.gt))
    # fg.stackplot(range(len(mTEPES.n)), -TechnologyCharge.loc[sc,p,:,:].values.reshape(len(mTEPES.n),len(mTEPES.gt)).transpose().tolist(), labels=list(mTEPES.gt))
    # fg.plot     (range(len(mTEPES.n)),  pDemand.sum(axis=1)*1e3, label='Demand', linewidth=0.2, color='k')
    fg.set(xlabel='Hours', ylabel='MW')
    plt.title(sc)
    fg.tick_params(axis='x', rotation=90)
    fg.legend()
    plt.tight_layout()
    plt.show()
    plt.savefig('oT_TechnologyOutput_'+sc+'_'+p+'_'+CaseName+'.png', bbox_inches='tight')

OutputResults = pd.Series(data=[mTEPES.vTotalOutput[sc,p,n,g]()*mTEPES.pDuration[n] for sc,p,n,g in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.g)))
OutputResults = pd.Series(data=[OutputResults[sc,p,n].filter(mTEPES.t2g[gt]).sum() for sc,p,n,gt in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.gt)))
OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_TechnologyEnergy_'+CaseName+'.csv', sep=',')

TechnologyEnergy = pd.Series(data=[0]*len(mTEPES.gt), index=list(mTEPES.gt))
for gt in mTEPES.gt:
    TechnologyEnergy[gt] = OutputResults.loc[:,:,:,gt].sum()

fig, fg = plt.subplots()
def func(percentage, TechnologyEnergy):
    absolute = int(percentage/100.*TechnologyEnergy.sum())
    return '{:.1f}%\n({:d} GWh)'.format(percentage, absolute)
wedges, texts, autotexts = plt.pie(TechnologyEnergy, labels=list(mTEPES.gt), autopct=lambda percentage: func(percentage, TechnologyEnergy), shadow=True, startangle=90, textprops=dict(color='w'))
plt.legend(wedges, list(mTEPES.gt), title='Energy Generation', loc='center left', bbox_to_anchor=(1, 0, 0.5, 1))
plt.setp(autotexts, size=8, weight='bold')
plt.axis('equal')
plt.show()
plt.savefig('oT_TechnologyEnergy_'+CaseName+'.png', bbox_inches='tight')

#%% outputting the network operation
OutputResults = pd.DataFrame.from_dict(mTEPES.vFlow.extract_values(), orient='index', columns=[str(mTEPES.vFlow)])
OutputResults *= 1e3
OutputResults.index = pd.MultiIndex.from_tuples(OutputResults.index)
OutputResults.index.names = ['Scenario','Period','LoadLevel','InitialNode','FinalNode','Circuit']
OutputResults = pd.pivot_table(OutputResults, values=str(mTEPES.vFlow), index=['Scenario','Period','LoadLevel'], columns=['InitialNode','FinalNode','Circuit'], fill_value=0)
OutputResults.index.names = [None] * len(OutputResults.index.names)
OutputResults.to_csv('oT_Result_NetworkFlow_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[abs(mTEPES.vFlow[sc,p,n,ni,nf,cc]()/mTEPES.pLineNTC[ni,nf,cc]) for sc,p,n,ni,nf,cc in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la], index= pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.la)))
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
pd.pivot_table(OutputResults, values=str(mTEPES.vENS), index=['Scenario','Period','LoadLevel'], columns=['Node'], fill_value=0).to_csv('oT_Result_NetworkPNS_'+CaseName+'.csv', sep=',')

OutputResults = pd.Series(data=[mTEPES.vENS[sc,p,n,nd]()*mTEPES.pDuration[n] for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd)))
OutputResults.to_frame(name='GWh').reset_index().pivot_table(index=['level_0','level_1','level_2'], columns='level_3', values='GWh').rename_axis(['Scenario','Period','LoadLevel'], axis=0).rename_axis([None], axis=1).to_csv('oT_Result_NetworkENS_'+CaseName+'.csv', sep=',')

#%% outputting the LSRMC
OutputResults = pd.Series(data=[mTEPES.dual[mTEPES.eBalance[sc,p,n,nd]]*1e3/mTEPES.pScenProb[sc]/mTEPES.pDuration[n] for sc,p,n,nd in mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd], index=pd.MultiIndex.from_tuples(list(mTEPES.sc*mTEPES.p*mTEPES.n*mTEPES.nd)))
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

#%% plotting the network in a map
import geopandas
import cartopy.crs as ccrs
from   cartopy.io import shapereader

# take data from http://www.naturalearthdata.com/
df = geopandas.read_file(shapereader.natural_earth(resolution='10m', category='cultural', name='admin_0_countries'))
polyAT = df.loc[df['ADMIN'] == 'Austria'       ]['geometry'].values[0]
polyBE = df.loc[df['ADMIN'] == 'Belgium'       ]['geometry'].values[0]
polyBG = df.loc[df['ADMIN'] == 'Bulgaria'      ]['geometry'].values[0]
polyCH = df.loc[df['ADMIN'] == 'Croatia'       ]['geometry'].values[0]
polyCY = df.loc[df['ADMIN'] == 'Cyprus'        ]['geometry'].values[0]
polyCZ = df.loc[df['ADMIN'] == 'Czechia'       ]['geometry'].values[0]
polyDK = df.loc[df['ADMIN'] == 'Denmark'       ]['geometry'].values[0]
polyEE = df.loc[df['ADMIN'] == 'Estonia'       ]['geometry'].values[0]
polyFI = df.loc[df['ADMIN'] == 'Finland'       ]['geometry'].values[0]
polyFR = df.loc[df['ADMIN'] == 'France'        ]['geometry'].values[0]
polyDE = df.loc[df['ADMIN'] == 'Germany'       ]['geometry'].values[0]
polyEL = df.loc[df['ADMIN'] == 'Greece'        ]['geometry'].values[0]
polyHU = df.loc[df['ADMIN'] == 'Hungary'       ]['geometry'].values[0]
polyIE = df.loc[df['ADMIN'] == 'Ireland'       ]['geometry'].values[0]
polyIT = df.loc[df['ADMIN'] == 'Italy'         ]['geometry'].values[0]
polyLV = df.loc[df['ADMIN'] == 'Latvia'        ]['geometry'].values[0]
polyLT = df.loc[df['ADMIN'] == 'Lithuania'     ]['geometry'].values[0]
polyLU = df.loc[df['ADMIN'] == 'Luxembourg'    ]['geometry'].values[0]
polyMT = df.loc[df['ADMIN'] == 'Malta'         ]['geometry'].values[0]
polyNL = df.loc[df['ADMIN'] == 'Netherlands'   ]['geometry'].values[0]
polyPL = df.loc[df['ADMIN'] == 'Poland'        ]['geometry'].values[0]
polyPT = df.loc[df['ADMIN'] == 'Portugal'      ]['geometry'].values[0]
polyRO = df.loc[df['ADMIN'] == 'Romania'       ]['geometry'].values[0]
polySK = df.loc[df['ADMIN'] == 'Slovakia'      ]['geometry'].values[0]
polySI = df.loc[df['ADMIN'] == 'Slovenia'      ]['geometry'].values[0]
polyES = df.loc[df['ADMIN'] == 'Spain'         ]['geometry'].values[0]
polySE = df.loc[df['ADMIN'] == 'Sweden'        ]['geometry'].values[0]
polyUK = df.loc[df['ADMIN'] == 'United Kingdom']['geometry'].values[0]

fg = plt.axes(projection=ccrs.PlateCarree())
fg.set_extent((min(mTEPES.pNodeLon.values())-2, max(mTEPES.pNodeLon.values())+2, min(mTEPES.pNodeLat.values())-2, max(mTEPES.pNodeLat.values())+2), crs=ccrs.PlateCarree())

# fg.add_geometries(polyAT, crs=ccrs.PlateCarree(), facecolor='C0', edgecolor='0.5')
# fg.add_geometries(polyBE, crs=ccrs.PlateCarree(), facecolor='C1', edgecolor='0.5')
# fg.add_geometries(polyBG, crs=ccrs.PlateCarree(), facecolor='C2', edgecolor='0.5')
# fg.add_geometries(polyCH, crs=ccrs.PlateCarree(), facecolor='C3', edgecolor='0.5')
# fg.add_geometries(polyCY, crs=ccrs.PlateCarree(), facecolor='C4', edgecolor='0.5')
# fg.add_geometries(polyCZ, crs=ccrs.PlateCarree(), facecolor='C5', edgecolor='0.5')
# fg.add_geometries(polyDK, crs=ccrs.PlateCarree(), facecolor='C6', edgecolor='0.5')
# fg.add_geometries(polyEE, crs=ccrs.PlateCarree(), facecolor='C7', edgecolor='0.5')
# fg.add_geometries(polyFI, crs=ccrs.PlateCarree(), facecolor='C8', edgecolor='0.5')
fg.add_geometries(polyFR, crs=ccrs.PlateCarree(), facecolor='C9', edgecolor='0.5')
# fg.add_geometries(polyDE, crs=ccrs.PlateCarree(), facecolor='C0', edgecolor='0.5')
# fg.add_geometries(polyEL, crs=ccrs.PlateCarree(), facecolor='C1', edgecolor='0.5')
# fg.add_geometries(polyHU, crs=ccrs.PlateCarree(), facecolor='C2', edgecolor='0.5')
# fg.add_geometries(polyIE, crs=ccrs.PlateCarree(), facecolor='C3', edgecolor='0.5')
# fg.add_geometries(polyIT, crs=ccrs.PlateCarree(), facecolor='C4', edgecolor='0.5')
# fg.add_geometries(polyLV, crs=ccrs.PlateCarree(), facecolor='C5', edgecolor='0.5')
# fg.add_geometries(polyLT, crs=ccrs.PlateCarree(), facecolor='C6', edgecolor='0.5')
# fg.add_geometries(polyLU, crs=ccrs.PlateCarree(), facecolor='C7', edgecolor='0.5')
# fg.add_geometries(polyMT, crs=ccrs.PlateCarree(), facecolor='C8', edgecolor='0.5')
# fg.add_geometries(polyNL, crs=ccrs.PlateCarree(), facecolor='C9', edgecolor='0.5')
# fg.add_geometries(polyPL, crs=ccrs.PlateCarree(), facecolor='C1', edgecolor='0.5')
fg.add_geometries(polyPT, crs=ccrs.PlateCarree(), facecolor='C2', edgecolor='0.5')
# fg.add_geometries(polyRO, crs=ccrs.PlateCarree(), facecolor='C3', edgecolor='0.5')
# fg.add_geometries(polySK, crs=ccrs.PlateCarree(), facecolor='C4', edgecolor='0.5')
# fg.add_geometries(polySI, crs=ccrs.PlateCarree(), facecolor='C5', edgecolor='0.5')
fg.add_geometries(polyES, crs=ccrs.PlateCarree(), facecolor='C0', edgecolor='0.5')
# fg.add_geometries(polySE, crs=ccrs.PlateCarree(), facecolor='C7', edgecolor='0.5')
# fg.add_geometries(polyUK, crs=ccrs.PlateCarree(), facecolor='C8', edgecolor='0.5')

# node name
for nd in mTEPES.nd:
    plt.annotate(nd, [mTEPES.pNodeLon[nd], mTEPES.pNodeLat[nd]])

#%% colors of the lines according to the ENTSO-E color code
# existing lines
for ni,nf,cc,lt in mTEPES.le*mTEPES.lt and mTEPES.pLineType:
   if lt == 'AC':
      if mTEPES.pLineVoltage[ni,nf,cc] > 700 and mTEPES.pLineVoltage[ni,nf,cc] <= 900:
          plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='blue'   , linewidth=1  , marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
      if mTEPES.pLineVoltage[ni,nf,cc] > 500 and mTEPES.pLineVoltage[ni,nf,cc] <= 700:
          plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='#ff8000', linewidth=1  , marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
      if mTEPES.pLineVoltage[ni,nf,cc] > 350 and mTEPES.pLineVoltage[ni,nf,cc] <= 500:
          plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='red'    , linewidth=0.5, marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
      if mTEPES.pLineVoltage[ni,nf,cc] > 290 and mTEPES.pLineVoltage[ni,nf,cc] <= 350:
          plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='green'  , linewidth=0.2, marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
      if mTEPES.pLineVoltage[ni,nf,cc] > 200 and mTEPES.pLineVoltage[ni,nf,cc] <= 290:
          plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='green'  , linewidth=0.2, marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
   else:
       plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='magenta', linewidth=0.5, marker='o', linestyle='solid' , transform=ccrs.PlateCarree())
# candidate lines
for ni,nf,cc,lt in mTEPES.lc*mTEPES.lt and mTEPES.pLineType:
    if lt == 'AC':
        if mTEPES.pLineVoltage[ni,nf,cc] > 700 and mTEPES.pLineVoltage[ni,nf,cc] <= 900:
            plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='blue'   , linewidth=1  , marker='o', linestyle='dashed', transform=ccrs.PlateCarree())
        if mTEPES.pLineVoltage[ni,nf,cc] > 500 and mTEPES.pLineVoltage[ni,nf,cc] <= 700:
            plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='#ff8000', linewidth=1  , marker='o', linestyle='dashed', transform=ccrs.PlateCarree())
        if mTEPES.pLineVoltage[ni,nf,cc] > 350 and mTEPES.pLineVoltage[ni,nf,cc] <= 500:
            plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='red'    , linewidth=0.5, marker='o', linestyle='dashed', transform=ccrs.PlateCarree())
        if mTEPES.pLineVoltage[ni,nf,cc] > 290 and mTEPES.pLineVoltage[ni,nf,cc] <= 350:
            plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='green'  , linewidth=0.2, marker='o', linestyle='dashed', transform=ccrs.PlateCarree())
        if mTEPES.pLineVoltage[ni,nf,cc] > 200 and mTEPES.pLineVoltage[ni,nf,cc] <= 290:
            plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='green'  , linewidth=0.2, marker='o', linestyle='dashed', transform=ccrs.PlateCarree())
    else:
        plt.plot([mTEPES.pNodeLon[ni], mTEPES.pNodeLon[nf]], [mTEPES.pNodeLat[ni], mTEPES.pNodeLat[nf]], color='magenta', linewidth=0.5, marker='o', linestyle='dashed', transform=ccrs.PlateCarree())

# line NTC
for ni,nf,cc in mTEPES.la:
    plt.annotate(round(mTEPES.pLineNTC[ni,nf,cc]*1e3), [(mTEPES.pNodeLon[ni]+mTEPES.pNodeLon[nf])/2, (mTEPES.pNodeLat[ni]+mTEPES.pNodeLat[nf])/2])

plt.title(CaseName+' Network Map')
plt.show()
plt.savefig('oT_Map_Network_'+CaseName+'.png', bbox_inches='tight')

PlottingNetMapsTime = time.time() - StartTime
print('Plotting network maps                 ... ', round(PlottingNetMapsTime), 's')
