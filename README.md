## Table of contents
* [General info](#general-info)
* [Technologies](#technologies)
* [Setup](#setup)

## General info
This script is to generate covalency from DFT population calculations. Specifically for transition metal atoms in calculations. 
	
## Technologies
Project is created with:
Python 3.6 (pandas, numpy)

cclib (https://cclib.github.io/)

Compatible with any QM calculations that cclib supports (e.g. Orca, Gaussian, molpro).

	
## Setup
Usage:
```
python WriteCovalencyReport.py --file "Drive:\Folder\where\log\file\is\CuI-Cl_1cluster_xes.out" --calc "LPA" --elem "Cu"

python WriteCovalencyReport.py --help 
```
-will show variable descriptions and default values.

Note that the --file should include the path and file name.

  -h, --help            show this help message and exit
  --file FILE           Full path of the file to be analyzed.
  --calc CALC           Calculation type to be used: LPA=Lowdin, MPA=Mulliken,
                        csquared=C-squared. LPA is the default.
  --elem ELEM           Element to be analyzed. Default=Fe.
  --gap GAP             Number of MOs to include around the HOMO-LUMO gap.
                        Default=10.
  --printLevel PRINTLEVEL
                        Minimum metal character percent to be printed in the
                        final report. Default=5.


Adding specific fragments:
1) Add a new regular expression variable in the section "Regular Expression Variables"
#Example adding for an Oxo:
```
regexO=re.compile('O.*S')
```
#if there was more than one oxygen then I must add the atom index I'm interested in:
```
regexO=re.compile('O5.*S')
```
#Note here 5 means it is the fifth element in the xyz part of the input for the calculation (starting from 1).

#Alternatively we can look at any atomic Oxo orbital (instead of just S):
```
regexO=re.compile('O5.*')
```
#Note: If you want to add very complicated statements (e.g. All nitrogens and some carbons) please look up more information on regular expressions for how to do this.

2) Store the atomic orbitals that match your search in section "Finding Atomic Orbitals":
```
OxoOrbs, matches = zip(*[(idx, string) for idx, string in enumerate(data.aonames) if  re.match('O', string)])
```
3) Add your atomic orbitals variable (here: OxoOrbs) to the Population Calculation in section Population Calculation:
```
m.calculate([S,Z2,XZ,YZ,X2Y2,XY,PX,PY,PZ,otherOrbs,OxoOrbs])
```
#note: keep track of what (zero-based) index your new orbitals are for later. Here OxoOrbs is number 10 in the list starting from zero.

4) Add your new orbitals to the alpha spin dataframe. Also in Section "Building Alpha Dataframe".
```
alp['Oxo_S']=character[10]
```
#10 because it was number 10 in the m.calculate list

5) Add this orbital to the final alpha spin dataframe. In section "Final Alpha Output".
```
    outDFalp=outDF[['Spin','Orbital Number', 'Energy (eV)','Occupancy',element+'_S',
             element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',
             element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character','Oxo_S']]
```
#Note: Here you can change the column orders by moving the names around. The only requirement is that the names match the names in Section: Building Alpha Dataframe. I've added the Oxo S orbital to the end.

7) Apply the changes to the Beta part of the calculation in Section: Building Beta dataframe.
bet['Oxo_S']=character[10]	
#note here the dataframe is now called "bet" for the beta orbitals.

8) Apply the new changes to the Section: "Final Beta Output"
outDFbet=outDF[['Spin','Orbital Number', 'Energy (eV)','Occupancy',element+'_S',
             element+'_DXY',element+'_DXZ',element+'_DYZ',element+'_DX2Y2',element+'_DZ2',
             element+'_PX',element+'_PY',element+'_PZ','Other_Character','Total_Metal_Character','Total_p_Character','Oxo_S]]
#Now in the final output file a column for Oxo_S orbitals will appear at the end of the list.    

```
