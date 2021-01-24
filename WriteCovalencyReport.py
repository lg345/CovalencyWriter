from __future__ import unicode_literals
import sys
from cclib.method import MPA,CSPA,LPA
from cclib.parser import ccData
from cclib.io import ccread
import fnmatch
import numpy as np
import re
import pandas as pd
import argparse
def WriteCovalencyReport(fileName,calculationType='LPA',element='Fe',homoLumoRange=15,printThreshold=5):
    filename=fileName
    data = ccread(filename)
    if calculationType=='LPA':
        m=LPA(data)#This is for Lowdin Analysis.
    elif calculationType=='MPA':
        m=MPA(data)#Mulliken Population Analysis
    elif calculationType=='csquared':
        m=CSPA(data)#C-squared Analysis
    else:
        Exception('Calculation type not supported.')
    print(element)
    #Section: Regular Expression Variables
    regexZ2=re.compile(element+'.*DZ2')
    regexXZ=re.compile(element+'.*DXZ')
    regexYZ=re.compile(element+'.*DYZ')
    regexX2Y2=re.compile(element+'.*DX2Y2')
    regexXY=re.compile(element+'.*DXY')
    regexPX=re.compile(element+'.*PX')
    regexPY=re.compile(element+'.*PY')
    regexPZ=re.compile(element+'.*PZ')
    regexS=re.compile(element+'.*S')
    
    #Section: Finding Atomic Orbitals
    Z2, matches1 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexZ2, string)])
    XZ, matches2 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexXZ, string)])
    YZ, matches3 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexYZ, string)])
    X2Y2, matches4 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexX2Y2, string)])
    XY, matches5 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexXY, string)])
    PX, matches6 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexPX, string)])
    PY, matches7 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexPY, string)])
    PZ, matches8 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexPZ, string)])
    S, matches9 = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if re.match(regexS, string)])
    #OxoOrbs, matches = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if  re.match('O', string)])#For oxo orbitals
    otherOrbs, matches = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if not (re.match(element, string) )])
    
    #Section: Population Calculation
    m.calculate([S,Z2,XZ,YZ,X2Y2,XY,PX,PY,PZ,otherOrbs])
    Spin=0#Alpha orbitals
    homo=data.homos[Spin]
    character=m.fragresults[Spin].T[:]*100
    indices=np.arange(0,len(data.aonames),1)
    orbital=[data.aonames[i][-1] for i in np.argmax(data.mocoeffs[Spin],1)]
    aos=[re.split(('_'),str(data.aonames[i])) for i in np.argmax(data.mocoeffs[Spin]**2,1)]
    AO=np.array(aos).T
    dominantAtom=np.array([data.aonames[i].split('_')[0] for i in np.argmax(data.mocoeffs[Spin],1)]).T
    dominantAO=np.array([data.aonames[i].split('_')[1] for i in np.argmax(data.mocoeffs[Spin],1)]).T
    energies=data.moenergies[Spin]
    energyArray=np.column_stack((indices,energies,dominantAtom,dominantAO))
    
    #Section: Building Alpha Dataframe.
    alp=pd.DataFrame(energyArray,columns=['Orbital Number','Energy (eV)','Dominant Atom', 'Dominant AO'])
    homo=data.homos[0]
    alp.loc[0:homo+1,'Occupancy']=1
    alp.loc[homo+1:,'Occupancy']=0
    alp['Spin']='Alpha'
    alp[element+'_S']=character[0]
    alp[element+'_DZ2']=character[1]
    alp[element+'_DXZ']=character[2]
    alp[element+'_DYZ']=character[3]
    alp[element+'_DX2Y2']=character[4]
    alp[element+'_DXY']=character[5]
    alp[element+'_PX']=character[6]
    alp[element+'_PY']=character[7]
    alp[element+'_PZ']=character[8]
    alp['Other_Character']=character[9]
    alp['Total_Metal_Character']=np.sum(character[0:6],axis=0)
    alp['Total_p_Character']=np.sum(character[6:8],axis=0)
    outDF=alp.loc[(alp.index > homo-homoLumoRange) & (alp.index < homo+homoLumoRange)& (alp['Total_Metal_Character']>printThreshold)].round(3)
    
    #Section: Final Alpha Output
    outDFalp=outDF[['Spin','Orbital Number', 'Energy (eV)','Occupancy',element+'_S',
             element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',
             element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character']]
    SumAlpha=outDFalp[['Spin',element+'_S',element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character']].loc[outDFalp['Occupancy']==0].sum(axis=0)
    
    Spin=1#beta orbitals
    homo=data.homos[Spin]
    character=m.fragresults[Spin].T[:]*100
    indices=np.arange(0,len(data.aonames),1)
    orbital=[data.aonames[i][-1] for i in np.argmax(data.mocoeffs[Spin],1)]
    aos=[re.split(('_'),str(data.aonames[i])) for i in np.argmax(data.mocoeffs[Spin]**2,1)]
    AO=np.array(aos).T
    dominantAtom=np.array([data.aonames[i].split('_')[0] for i in np.argmax(data.mocoeffs[Spin],1)]).T
    dominantAO=np.array([data.aonames[i].split('_')[1] for i in np.argmax(data.mocoeffs[Spin],1)]).T
    energies=data.moenergies[Spin]
    energyArray=np.column_stack((indices,energies,dominantAtom,dominantAO))
    
    #Section: Building Beta Dataframe.    
    bet=pd.DataFrame(energyArray,columns=['Orbital Number','Energy (eV)','Dominant Atom', 'Dominant AO'])
    bet.loc[0:homo+1,'Occupancy']=1
    bet.loc[homo+1:,'Occupancy']=0
    bet['Spin']="Beta"
    bet[element+'_S']=character[0]
    bet[element+'_DZ2']=character[1]
    bet[element+'_DXZ']=character[2]
    bet[element+'_DYZ']=character[3]
    bet[element+'_DX2Y2']=character[4]
    bet[element+'_DXY']=character[5]
    bet[element+'_PX']=character[6]
    bet[element+'_PY']=character[7]
    bet[element+'_PZ']=character[8]
    #bet['O_Character']=character[8]
    bet['Other_Character']=character[9]
    bet['Total_Metal_Character']=np.sum(character[0:6],axis=0)
    bet['Total_p_Character']=np.sum(character[6:8],axis=0)
    outDF=bet.loc[(bet.index > homo-homoLumoRange) & (bet.index < homo+homoLumoRange)& (bet['Total_Metal_Character']>printThreshold)].round(3)
    #Section: Final Beta Output
    outDFbet=outDF[['Spin','Orbital Number', 'Energy (eV)','Occupancy',element+'_S',
             element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',
             element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character']]
    SumBeta=outDFbet[['Spin',element+'_S',element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character']].loc[outDFbet['Occupancy']==0].sum(axis=0)
    
    #Section: Write Final Ouput File
    TotalDF=outDFalp.append(outDFbet)
    SumTotal=(SumAlpha+SumBeta)[1:]
    f = open(fileName.split('.')[0]+'.'+calculationType+'.report', 'w')
    f.write('Alpha Hole Character Sum:\n'+str(SumAlpha[1:])[:-13]+'\nBeta Hole Character Sum:\n'+str(SumBeta[1:])[:-13]+'\nTotal Hole Sum:\n'+str(SumTotal[1:])[:-13])
    TotalDF.to_csv(f,sep='\t')
    return(TotalDF)
    f.close()
    
if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Script to write metal-based covalency report from QM calculations. Contact: Leland Gee (lbgee@stanford.edu)')
    parser.add_argument('--file',help='Full path of the file to be analyzed.')
    parser.add_argument('--calc',help='Calculation type to be used: LPA=Lowdin, MPA=Mulliken, csquared=C-squared. Default=LPA.',type=str,default='LPA')
    parser.add_argument('--elem',help='Element to be analyzed. Default=Fe.',type=str,default='Fe')
    parser.add_argument('--gap',help='Number of MOs to include around the HOMO-LUMO gap. Default=10.',type=int,default=10)
    parser.add_argument('--printLevel',help='Minimum metal character percent to be printed in the final report. Default=5.',type=float,default=5)
    args = parser.parse_args()
    WriteCovalencyReport(args.file,args.calc,args.elem,args.gap,args.printLevel) 
