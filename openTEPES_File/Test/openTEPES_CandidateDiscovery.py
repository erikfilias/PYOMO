# Open Generation and Transmission Operation and Expansion Planning Model with RES and ESS (openTEPES) - Version 1.6.33 - June 17, 2020

pBoundOnCandidates     = 50

import time
import math
import pandas as pd
from   pyomo.environ import DataPortal, Set

print('Candidate discovery     ****')

StartTime = time.time()

# The way to use this is to run it and wait for the file UpgradedNetwork to be generated
# Then, upgraded network is fed to the model in a subsequent iteration
# Transformers, converters could be added if we add node voltage
# We could introduce the more sophisticated version based on angles, but it seems that pMaxTheta is not defined yet
# Bound on candidates is established by hand, we could also introduce it in the options file
# In this version a linear relaxation of the problem is solved first. Normally, we would restrict investment at all, but I could not find a parameter for this. Is there a reason why it is not defined?

# Initializing

dictSets = DataPortal()
dictSets.load(filename='oT_Dict_Cable_'+CaseName+'.csv', set='cb', format='set')
mTEPES.cb = Set(initialize=dictSets['cb'], ordered=False, doc='cables'    )

# this union of sets gives a warning. To avoid this warning, use block.del_component() and block.add_component().
mTEPES.cc = mTEPES.cc | mTEPES.cb

dfNetwork = pd.read_csv('oT_Data_Network_'+CaseName+'.csv', index_col=[0,1,2])
dfCables  = pd.read_csv('oT_Data_Cable_'  +CaseName+'.csv', index_col=[0    ])
dfNetwork.fillna(0, inplace=True)
dfCables.fillna (0, inplace=True)

pLinearDistance              = pd.DataFrame(0, index=range(len(mTEPES.nd*mTEPES.nd)), columns=['ni','nf',     'LinearDistance'       ])
pPotentialBenefit            = pd.DataFrame(0, index=range(len(mTEPES.nd*mTEPES.nd)), columns=['ni','nf',     'PotentialBenefit'     ])
pPotentialBenefitRatio       = pd.DataFrame(0, index=range(len(mTEPES.nd*mTEPES.nd)), columns=['ni','nf','cc','PotentialBenefitRatio'])
pPotentialBenefitRatioSorted = pd.DataFrame(0, index=range(len(mTEPES.nd*mTEPES.nd)), columns=['ni','nf','cc','PotentialBenefitRatio'])

AuxIndex  = 0
AuxIndex2 = 0

# For every pair of nodes, calculating potential benefit and cost
for ni,nf in mTEPES.ni*mTEPES.nf:
    if ni > nf:
        ThisDistance = 2 * math.asin(math.sqrt(math.sin((mTEPES.pNodeLat[nf]-mTEPES.pNodeLat[ni])*math.pi/180/2)**2 + math.cos(mTEPES.pNodeLat[ni]*math.pi/180)*math.cos(mTEPES.pNodeLat[nf]*math.pi/180)*(math.sin((mTEPES.pNodeLon[nf] - mTEPES.pNodeLon[ni])*math.pi/180/2)**2)))
        pLinearDistance.loc[AuxIndex] = [ni,nf,ThisDistance]

        # Potential Benefit per unit capacity
        ThisBenefit = sum(abs(mTEPES.dual[mTEPES.eBalance[sc,p,n,nf]]-mTEPES.dual[mTEPES.eBalance[sc,p,n,ni]])*1e3/mTEPES.pScenProb[sc]/mTEPES.pDuration[n] for sc,p,n in mTEPES.sc*mTEPES.p*mTEPES.n)
        pPotentialBenefit.loc[AuxIndex] = [ni,nf,ThisBenefit]

        for cc in mTEPES.cc:
            if ThisDistance <= dfCables.loc[cc,'MaximumDistanceCT'].item() and ThisDistance >= dfCables.loc[cc,'MinimumDistanceCT'].item():
                # Cost of candidate
                ThisCost     = dfCables.loc[cc,'FixedCostperKmCT'].item()*ThisDistance
                # Potential Benefit ratio using the capacity
                ThisCapacity = dfCables.loc[cc,'TTCCT'].item()
                if ThisCost > 0:
                    pPotentialBenefitRatio.loc[AuxIndex2] = [ni, nf, cc, ThisBenefit * ThisCapacity/(ThisCost+1e-7)]
                    AuxIndex2 = AuxIndex2 + 1

        AuxIndex = AuxIndex +1

pLinearDistance.set_index       (['ni','nf'     ])
pPotentialBenefit.set_index     (['ni','nf'     ])
pPotentialBenefitRatio.set_index(['ni','nf','cc'])

pNetworkNew = pd.DataFrame(0, range(pBoundOnCandidates), columns=['ni','nf','cc','LineType','Voltage','LossFactor','Reactance','TTC','SecurityFactor','FixedCost','FixedChargeRate','BinaryInvestment'])

# Sorting and selecting candidates
# Evitar la escritura en el fichero csv
pPotentialBenefitRatio.to_csv                   ('oT_PotentialBenefitRatio_'+CaseName+'.csv')
pPBR                               = pd.read_csv('oT_PotentialBenefitRatio_'+CaseName+'.csv', index_col=[0, 1 ,2, 3])
pPotentialBenefitRatioSorted = pPBR.sort_values(by=['PotentialBenefitRatio'], ascending=False)
pSorted = pPotentialBenefitRatioSorted[0:pBoundOnCandidates]

listni = pSorted.index.get_level_values('ni')
listnf = pSorted.index.get_level_values('nf')
listnc = pSorted.index.get_level_values('cc')

# Adding the information to the lines - candidate lines
pNumberLines = len(dfNetwork.index)
for i in range(len(pSorted)):

    ni = listni[i]
    nf = listnf[i]
    cc = listnc[i]

    ThisDistance         = pLinearDistance.loc[(pLinearDistance.ni==ni)].loc[(pLinearDistance.nf==nf), 'LinearDistance'].item()
    ThisLineType         = dfCables.loc[cc,'LineType'        ]
    ThisVoltage          = dfCables.loc[cc,'Voltage'         ]
    ThisLossFactor       = dfCables.loc[cc,'LossFactor'      ]
    ThisReactance        = dfCables.loc[cc,'XperKmCT'        ] * ThisDistance
    ThisTTC              = dfCables.loc[cc,'TTCCT'           ]
    ThisSecurityFactor   = dfCables.loc[cc,'SecurityFactor'  ]
    ThisFixedCost        = dfCables.loc[cc,'FixedCostperKmCT'] * ThisDistance
    ThisFixedChargeRate  = dfCables.loc[cc,'FixedChargeRate' ]
    ThisBinaryInvestment = 1
    pNetworkNew.loc[i-1]      = [ni,nf,cc,ThisLineType,ThisVoltage,ThisLossFactor,ThisReactance,ThisTTC,ThisSecurityFactor,ThisFixedCost,ThisFixedChargeRate,ThisBinaryInvestment]
    dfNetwork.loc[(ni,nf,cc)] = [         ThisLineType,ThisVoltage,ThisLossFactor,ThisReactance,ThisTTC,ThisSecurityFactor,ThisFixedCost,ThisFixedChargeRate,ThisBinaryInvestment]

dfNetwork.to_csv('oT_Data_UpgradedNetwork_'+CaseName+'.csv')

CandidateTime = time.time() - StartTime
StartTime   = time.time()
print('Candidate discovery                   ... ', round(CandidateTime), 's')
