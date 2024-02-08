

##
'''
---Abaqus-Python script for testing with Hill'48, PolyN and Neural-Network based Yield funcs: PLANE STRESS
---Developed by Stefan C. Soare 
---Added NN-based yield funcs to the previous PolyN-version of this scripts (https://github.com/stefanSCS/PolyN)
---License: see specs on the GitHub repo of these scripts
---Released by: S.C.S
---Theory in article: "On the use of neural networks in the modeling of yield surfaces"

---Note: Abaqus uses Python 2.7
---Element tested: S4R 
---Variable thisPATH stores path (to this script)
---All length dimensions are in mmm
---All stress dimensions are in MPa
---DisplacementBC's for uniaxial tests are infered from Swift/Voce hardening law 
---Further details are in docs: 'NNys_Usage.pdf'
'''
##
from __future__ import print_function
##
abqSims=True  ##run with Abaqus (for sims)
pytPlot=True   ##run  with Python (for post-proc: plots)    
#####----------------------------------------------------------
####------------------------------------------------------------
###Import all abaqus modules (adds some extra-imports)
try:
    from abaqus import *
    from abaqusConstants import *
    from caeModules import *
    from odbAccess import *
    #import regionToolset
except ImportError:
    abqSims,pytPlot=False,True
try: 
    import numpy as znp
    from matplotlib import pyplot as plt
except ImportError:
    pytPlot=False
##    
###Import other Python modules
from math import sqrt as msqrt
from math import exp as mexp
from math import log as mlog
from math import cos as mcos
from math import sin as msin
from math import pi as mpi
from os import chdir as oschdir
import os.path as ospath
from os import name as osn
from os import remove as osremove
from os import makedirs as osmaked
from os import getcwd as osgetcwd
from sys import __stdout__ as sysOut
from sys import exit as sysexit
#from sys import version_info 
#verprt=version_info[0]<3
#
##Abaqus scratch and save directory 
thisPATH="./"
figDirData="./DATA/"
figDirPlot="./RESULTS/"
if(osn=='nt'):
    thisPATH=".\\"
    figDirData=".\\DATA\\"
    figDirPlot=".\\RESULTS\\"    
#
umatHill='abqUMAT_Hill48.f'
umatPoly='abqUMAT_PolyN.f'
umatNNYS='abqUMAT_nnYS.f'
if(osn=='nt'):
    umatHill='abqUMAT_Hill48.for'
    umatPoly='abqUMAT_PolyN.for'
    umatNNYS='abqUMAT_nnYS.for'


caeOnly=False
uaxSim=True
nPlotPoints=201
##Abaqus scratch and save directory 
oschdir(thisPATH)
##this script for Tensile-Test 
thisName="abqPolyN"
printThisName=thisName+':'
modelName=thisPATH+thisName+"uax.cae"
mdbModelName='Model-1'
partName='Part-1'
partInstanceName='Part-1-1'
rad=mpi/180.0
noGUI=True
if(pytPlot):noGUI=False
dtMax=0.01
maxNInc=1000
HRATIO=6.0
nElemMultiplier=1
## Shell element option (for 2D, this is the only interesting option)
shell_S4R=True
#
aMAT={'name':'Material-1','descrp':'','hThick':0.0,'E':0.0, 'nuPoisson':0.0, 'muG':0.0,
     'hardLaw':'','hardParam':[0.0,0.0,0.0],'rVals':[0.0,0.0,0.0],'rExp':[],
     'blankDiam':0.0,'dieOpeningDiam':0.0,'dieShoulderRad':0.0,'punchDiam':0.0,
     'punchNoseRad':0.0,'holderClearance':0.0,'frictionCoeff':0.0,
     'punchTravel':0.0, 'overrideDefault':False,
     'seedScaleC':1.0, 'seedScaleR':1.0, 'pSoft': [1.0,1.0,1.0], 'scaleP': 1.05, 'deltaPunch': 3.0, 
     'readOdb': False, 'nCPUs': int(1),
     'dtInitial':0,'dtMax':0,'maxNInc':0,
     'subDir':''}
hhh=0.0 ##aMAT['hThick']
transvK=0.0 ##(5.0*aMAT['muG']*aMAT['hThick'])/6.0 ##transverse shear shell   
## if loadType==1/0 then BC=displacement/force 
loadType=True
##sample dimensions (length, width) 
WW=10.0  ##;WW=13.2773
LL=10.0*WW
### zero and geometric tolerance
zro,gTol=0.000, 1.0e-5
#
if(not shell_S4R):##S8R uses Green-Lagrange strain (L)
    for kk in range(len(vStrains)):
        vStrains[kk]=0.5*(mexp(2.0*vStrains[kk])-1.0)

###Select element type (if shell_S4R == 0, then S8R is used) 
####shell_S4R=True
if(shell_S4R):
    elmType='S4R'
else:     
    elmType='S8R'    
#
##strUMAT=''

dtInitial=dtMax/20.0
JOB_DESCRIPTION=''
JOB_NAME = 'JOB-Tensile-Test'
vvTests=[]
#
# 
def printMsg(msg):
    if(noGUI):
        #print >> sysOut, printThisName+str(msg)
        print(str(msg),file=sysOut)
    #elif(verprt):
    #    #print msg
    #    pass
    else:    
        print(msg)    
#
def setUaxCae(uax=True,ZcaeOnly=False):
    global uaxSim,caeOnly
    uaxSim,caeOnly=uax,ZcaeOnly
 
def readData(bParam,inpData,ZHRATIO=6.0,ZnElemMultiplier=1):
    global aMAT,hhh,transvK,vvTests,dtMax,dtInitial,maxNInc,HRATIO,nElemMultiplier,figDirData,figDirPlot
    ZsubDir=bParam['subDir']
    if(ZsubDir):
        if(osn=='nt'):figDirData+=ZsubDir+'\\'
        else:figDirData+=ZsubDir+'/'    
    if(not (ospath.exists(figDirData) and ospath.isdir(figDirData))):
        printMsg("The local data folder was not found")
        printMsg("Make sure a folder {} exists at your script location".format(figDirData))
        printMsg("Calculations aborted");sysexit()    
    if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
        printMsg("The local folder for saving reports and plots was not found")
        printMsg("Make sure a folder \'RESULTS\' exists at your script location")
        printMsg("Calculations aborted");sysexit()
    if(ZsubDir):
        if(osn=='nt'):figDirPlot+=ZsubDir+'\\'
        else:figDirPlot+=ZsubDir+'/'
        if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
            try:
                osmaked(figDirPlot)
            except OSError as err:
                printMsg('readData:');printMsg(err);sysexit()            
    HRATIO=ZHRATIO
    nElemMultiplier=ZnElemMultiplier
    if(uaxSim):dtMax,maxNInc=bParam['dtMax'],bParam['maxNInc']
    if(bParam['hLaw'] in ['Swift','Voce']):
        aMAT['hardLaw']=bParam['hLaw']
        if(bParam['hLaw']=='Voce'):dtInitial=dtMax/100.0
        else:dtInitial=dtMax/50.0
    else:
        printMsg("Only 'Swift' and 'Voce' hardening laws can be used (at the moment)")
        printMsg('Calculations aborted');sysexit()        
    aMAT['hThick']=bParam['hThick']
    aMAT['E'],aMAT['nuPoisson'],aMAT['muG']=bParam['eParam']
    aMAT['hardParam'][0],aMAT['hardParam'][1],aMAT['hardParam'][2]=bParam['hParam']
    aMAT['rVals'][0],aMAT['rVals'][1],aMAT['rVals'][2]=bParam['rParam']
    hhh=bParam['hThick']
    transvK=(5.0*aMAT['muG']*aMAT['hThick'])/6.0 ##transverse shear shell
    ##printMsg(str(transvK));sysexit()
    if(uaxSim):
        for rrr in bParam['rExp']:
            aMAT['rExp'].append(rrr)
    else:
        try:
            aMAT['blankDiam']=float(bParam['blankDiam'])
            aMAT['dieOpeningDiam']=float(bParam['dieOpeningDiam']) 
            aMAT['dieShoulderRad']=float(bParam['dieShoulderRad'])
            aMAT['punchDiam']=float(bParam['punchDiam'])
            aMAT['punchNoseRad']=float(bParam['punchNoseRad'])
            aMAT['holderClearance']=float(bParam['holderClearance'])
            aMAT['frictionCoeff']=float(bParam['frictionCoeff']) 
            ##aMAT['punchTravel']=float(bParam['punchTravel'])
            ##aMAT['overrideDefault']=bool(bParam['overrideDefault'])
        except Exception as erx:
            printMsg(erx);sysexit()            
    for ddd in inpData:
        ufile,dfile='',ddd['inpFile']
        if(not dfile):printMsg("field 'inpFile' cannot be empty\nCalculations aborted");sysexit()
        if(ddd['UMAT']):
            ufile=thisPATH+ddd['UMAT']
            if(not(ospath.exists(ufile) and ospath.isfile(ufile))):
                printMsg('File '+ufile+' not found\nCalculations aborted')
                sysexit()
        if(ddd['PolyN']):
            zdfile=figDirData+dfile
            if(not(ospath.exists(zdfile) and ospath.isfile(zdfile))):
                printMsg('File '+dfile+' not found\nCalculations aborted')
                sysexit()
        tdd={'UMAT':ddd['UMAT'],'ufile':ufile,'PolyN':ddd['PolyN'],'dfile':dfile,'angles':[]}        
        if(uaxSim):
            tdd['angles']=ddd['angles']    
        vvTests.append(tdd)
    ###print(vvTests)        
    printMsg('----------------------------------------------------------------input data: OK')
    ###sysexit()    
#

def readInpData(fName,subDir):  ##(bParam,inpData,ZHRATIO=6.0,ZnElemMultiplier=1):
    global aMAT,hhh,transvK,vvTests,dtMax,dtInitial,maxNInc,HRATIO,nElemMultiplier,figDirData,figDirPlot
    ZsubDir=subDir
    if(ZsubDir):
        if(osn=='nt'):figDirData+=ZsubDir+'\\'
        else:figDirData+=ZsubDir+'/'    
    if(not (ospath.exists(figDirData) and ospath.isdir(figDirData))):
        printMsg("The local data folder was not found")
        printMsg("Make sure a folder {} exists at your script location".format(figDirData))
        printMsg("Calculations aborted");sysexit()    
    if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
        printMsg("The local folder for saving reports and plots was not found")
        printMsg("Make sure a folder \'RESULTS\' exists at your script location")
        printMsg("Calculations aborted");sysexit()
    if(ZsubDir):
        if(osn=='nt'):figDirPlot+=ZsubDir+'\\'
        else:figDirPlot+=ZsubDir+'/'
        if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
            try:
                osmaked(figDirPlot)
            except OSError as err:
                printMsg('readInpData:');printMsg(err);sysexit()
    try:
        fp=open(figDirData+fName)
    except IOError:printMsg('readInpData: input file not found');sysexit()
    bParam={'uax':True,'caeOnly':True,'umatHill':'', 'umatPoly':'','sim':[],
           'dtInitial':1.0e-5,'dtMax':0.001,'maxNInc':10000,'hLaw':'','eParam':[],
           'hParam':[],'rParam':[],'rExp':[[],[]],'hThick':0,'HRATIO':0,
           'blankDiam':0, 'dieOpeningDiam':0, 'dieShoulderRad':0, 'punchDiam':0,
           'punchNoseRad':0, 'holderClearance':0, 'frictionCoeff':0,
           'seedScaleC':1.0, 'seedScaleR':1.0, 'pSoft': [], 
           'scaleP': 1.05, 'deltaPunch': 2.0, 'readOdb': False, 'nCPUs': int(1)}       
    bKeys=list(bParam.keys())
    dbool={'False':False,'True':True, '0':False, '1':True};kbool=dbool.keys()    
    msg='Calculations aborted'
    inpData,tuax,tcae=[],True,True    
    for line in fp:
        line=line.strip()
        if((line=='') or (line[0]=='#')):continue
        if(':' in line):
            line=line.split(':');lnz=line[0].strip()
        else:printMsg(fName+': Unknown format on line');printMsg(line);printMsg(msg);sysexit()
        if(lnz not in bKeys):printMsg(fName+': Unknown option on line');printMsg(line);printMsg(msg);sysexit()
        tpval=line[1].strip()
        if(lnz=='uax'):
            if(tpval in kbool):tuax=dbool[tpval]
            else:printMsg(line);printMsg('Unknown option');printMsg(msg);sysexit()
            continue
        if(lnz=='caeOnly'):
            if(tpval in kbool):tcae=dbool[tpval]
            else:printMsg(line);printMsg('Unknown option');printMsg(msg);sysexit()
            continue            
        if(lnz in ['dtInitial','dtMax']):
            try:bParam[lnz]=float(tpval)
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz=='maxNInc'):
            try:bParam[lnz]=int(float(tpval))
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz=='hLaw'):
            if(tpval in ['Swift','Voce']):bParam[lnz]=tpval
            else:printMsg('Unknown hardening option');printMsg(msg);sysexit()
            continue
        if(lnz=='eParam'):
            vtp=tpval.split(',')
            try:
                for item in vtp:bParam[lnz].append(float(item))
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            if(len(bParam[lnz]) <3):printMsg('Must specify 3 elastic constants');printMsg(msg);sysexit()
            continue            
        if(lnz=='hParam'):
            vtp=tpval.split(',')
            try:
                for item in vtp:bParam[lnz].append(float(item))
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            if(len(bParam[lnz])<3):printMsg('Must specify 3 hardening constants');printMsg(msg);sysexit()
            continue
        if(lnz=='rParam'):
            vtp=tpval.split(',')
            try:
                for item in vtp:bParam[lnz].append(float(item))
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            if(len(bParam[lnz])<3):printMsg('Must specify 3 experimental r-values (0,45,90)');printMsg(msg);sysexit()
            continue
        if(lnz=='hThick'):
            try:bParam[lnz]=float(tpval)
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz=='HRATIO'):
            try:bParam[lnz]=float(tpval)
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz in ['blankDiam','dieOpeningDiam','dieShoulderRad','punchDiam','punchNoseRad','holderClearance','frictionCoeff']):
            try:bParam[lnz]=float(tpval)
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz in ['seedScaleC','seedScaleR','scaleP','deltaPunch']):
            try:bParam[lnz]=float(tpval)
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            continue
        if(lnz=='pSoft'):
            vtp=tpval.split(',');vtp=[item.strip() for item in vtp]
            if(len(vtp)<3):printMsg(tpval);printMsg('Incorrect number of data items');printMsg(msg);sysexit()
            try:
                for item in vtp:bParam['pSoft'].append(float(item))
            except ValueError as err:printMsg(tpval);printMsg(err);sysexit()
            continue
        if(lnz=='readOdb'):
            #if(tpval in kbool):bParam[lnz]=dbool[tpval]
            if(tpval in ['0','3','4']):bParam[lnz]=int(tpval)
            else:printMsg(line);printMsg('Unknown option');printMsg(msg);sysexit()
            continue
        if(lnz=='nCPUs'):
            try:bParam[lnz]=int(float(tpval))
            except ValueError as err:printMsg(lnz);printMsg(err);printMsg(msg);sysexit()
            if(bParam[lnz]<1):bParam[lnz]=1            
        if(lnz=='sim'):
            vtp=tpval.split('|')
            vtp=[item.strip() for item in vtp]
            if(tuax and (len(vtp)<4)):printMsg(tpval);printMsg('Incorrect number of data items');printMsg(msg);sysexit()
            elif(len(vtp)<3):printMsg(tpval);printMsg('Incorrect number of data items:');printMsg(msg);sysexit()
            dtemp={'inpFile':vtp[0],'UMAT':False,'PolyN':False,'angles':[],'digProfile':[],'digFileSep':[]}
            if(vtp[1]=='False'):dtemp['UMAT']=False
            elif(vtp[1]=='umatHill'):dtemp['UMAT']=umatHill
            elif(vtp[1]=='umatPoly'):dtemp['UMAT']=umatPoly
            elif(vtp[1]=='umatNNYS'):dtemp['UMAT']=umatNNYS
            else:printMsg(tpval);printMsg('Unknown option for UMAT');printMsg(msg);sysexit()
            if(vtp[2] in kbool):dtemp['PolyN']=dbool[vtp[2]]
            else:printMsg(tpval);printMsg('Unknown option for PolyN');printMsg(msg);sysexit()
            if(tuax):
                vtp=vtp[3].strip(',').split(',')
                try:
                    for item in vtp: dtemp['angles'].append(float(item))
                except ValueError as err:printMsg(tpval);printMsg(err);printMsg(msg);sysexit()    
            inpData.append(dtemp)                            
    fp.close()
    setUaxCae(tuax,tcae)    
    HRATIO=bParam['HRATIO']
    nElemMultiplier=int(1) ##ZnElemMultiplier
    dtMax,maxNInc=bParam['dtMax'],bParam['maxNInc']
    ##aMAT['subDir']=ZsubDir
    if(bParam['hLaw'] in ['Swift','Voce']):
        aMAT['hardLaw']=bParam['hLaw']
        if(bParam['hLaw']=='Voce'):dtInitial=dtMax/100.0
        else:dtInitial=dtMax/50.0
    elif(bParam['hLaw']):
        printMsg("Only 'Swift' and 'Voce' hardening laws can be used (at the moment)")
        printMsg('Calculations aborted');sysexit()        
    aMAT['hThick']=bParam['hThick']
    if(bParam['eParam']):aMAT['E'],aMAT['nuPoisson'],aMAT['muG']=bParam['eParam']
    if(bParam['hParam']):aMAT['hardParam'][0],aMAT['hardParam'][1],aMAT['hardParam'][2]=bParam['hParam']
    if(bParam['rParam']):aMAT['rVals'][0],aMAT['rVals'][1],aMAT['rVals'][2]=bParam['rParam']
    hhh=bParam['hThick']
    transvK=(5.0*aMAT['muG']*aMAT['hThick'])/6.0 ##transverse shear shell
    ##printMsg(str(transvK));sysexit()
    if(uaxSim):
        for rrr in bParam['rExp']:
            aMAT['rExp'].append(rrr)
    else:
        try:
            aMAT['blankDiam']=float(bParam['blankDiam'])
            aMAT['dieOpeningDiam']=float(bParam['dieOpeningDiam']) 
            aMAT['dieShoulderRad']=float(bParam['dieShoulderRad'])
            aMAT['punchDiam']=float(bParam['punchDiam'])
            aMAT['punchNoseRad']=float(bParam['punchNoseRad'])
            aMAT['holderClearance']=float(bParam['holderClearance'])
            aMAT['frictionCoeff']=float(bParam['frictionCoeff']) 
            ##aMAT['punchTravel']=float(bParam['punchTravel'])
            ##aMAT['overrideDefault']=bool(bParam['overrideDefault'])
            aMAT['seedScaleC']=bParam['seedScaleC']
            aMAT['seedScaleR']=bParam['seedScaleR']
            aMAT['pSoft']=bParam['pSoft'] 
            aMAT['scaleP']=bParam['scaleP']
            aMAT['deltaPunch']=bParam['deltaPunch']
            aMAT['readOdb']=bParam['readOdb']
            aMAT['nCPUs']=bParam['nCPUs']
            aMAT['dtInitial']=bParam['dtInitial']
            aMAT['dtMax']=bParam['dtMax']
            aMAT['maxNInc']=bParam['maxNInc']
        except Exception as erx:
            printMsg(erx);sysexit()            
    for ddd in inpData:
        ufile,dfile='',ddd['inpFile']
        if(not dfile):printMsg("field 'inpFile' cannot be empty\nCalculations aborted");sysexit()
        if(ddd['UMAT']):
            ufile=thisPATH+ddd['UMAT']
            if(not(ospath.exists(ufile) and ospath.isfile(ufile))):
                printMsg('File '+ufile+' not found\nCalculations aborted')
                sysexit()
        if(ddd['PolyN']):
            zdfile=figDirData+dfile
            if(not(ospath.exists(zdfile) and ospath.isfile(zdfile))):
                printMsg('File '+dfile+' not found\nCalculations aborted')
                sysexit()
        tdd={'UMAT':ddd['UMAT'],'ufile':ufile,'PolyN':ddd['PolyN'],'dfile':dfile,'caseNNys':False,'angles':[]}
        if(ddd['UMAT']==umatNNYS):tdd['caseNNys']=True        
        if(uaxSim):tdd['angles']=ddd['angles']    
        vvTests.append(tdd) 
    ####print(vvTests)        
    printMsg('----------------------------------------------------------------input data: OK')
    return inpData



def ffHardFunc(hLaw,hParam,veps):
    if(hLaw=='Swift'):###Swift: sigma=A*(B+eps)^C
        return hParam[0]*(hParam[1]+veps)**hParam[2]
    elif(hLaw=='Voce'):###Voce: A-B*exp(-C*eps)
        #return hParam[0]-hParam[1]*uax.znp.exp(-hParam[2]*veps)
        return hParam[0]-hParam[1]*znp.exp(-hParam[2]*veps)
    else:
        print("hardFunc: Hardening option must be one of ['Swift','Voce']")
        #uax.sysexit()
        sysexit()        

def ffHRat(hLaw,hParam,eps):
    if(hLaw=='Swift'):
        return ((hParam[1]+eps)/hParam[1])**(1.0-hParam[2])
    elif(hLaw=='Voce'):
        #return uax.znp.exp(hParam[2]*eps)  
        return znp.exp(hParam[2]*eps)          

def ffdHRat(hLaw,hParam,eps):
    if(hLaw=='Swift'):
        return ((1.0-hParam[2])/hParam[1])*(hParam[1]/(hParam[1]+eps))**(hParam[2])
    elif(hLaw=='Voce'):
        #return hParam[2]*uax.znp.exp(hParam[2]*eps)
        return hParam[2]*znp.exp(hParam[2]*eps)

def ffplotHrat(hLaw,hParam,eB=0.5):
    if(not pytPlot):
        print('plotHrat: No plotting in sim-mode')
        #uax.sysexit()
        sysexit()
    #veps=uax.znp.linspace(0.0,eB,1+10**4)
    veps=znp.linspace(0.0,eB,1+10**4)
    vsig=ffHardFunc(hLaw,hParam,veps)
    vhrat=ffHRat(hLaw,hParam,veps)
    vdhrat=ffdHRat(hLaw,hParam,veps)
    swiftBound=34.0
    voceBound=50.2
    km=0
    if(hLaw=='Voce'):
        for kk in range(vdhrat.shape[0]):
            if(vdhrat[kk]>voceBound):
                km=kk;print('km = ',km) 
                print('dhrat = {}'.format(vdhrat[kk]),
                ', eps = ',veps[kk], ', hratio = ',ffHRat(hLaw,hParam,veps[kk]))
                break
    elif(hLaw=='Swift'):
        for kk in range(vdhrat.shape[0]):
            if(vdhrat[kk]<swiftBound):
                km=kk;print('km = ',km) 
                print('dhrat = {}'.format(vdhrat[kk]),
                ', eps = ',veps[kk], ', HRATIO = ',ffHRat(hLaw,hParam,veps[kk]))
                break           
    if((km==0) or (km>=vdhrat.shape[0]-2)):
        print("ffplotHrat: No suitable 'km' was found")
        if(hLaw=='Swift'):print("ffplotHrat: Adjust value of variable 'swiftBound'")
        if(hLaw=='Voce'):print("ffplotHrat: Adjust value of variable 'voceBound'")        
    #fg=uax.plt.figure()
    fg=plt.figure(figsize=(12,6))
    ax1=fg.add_subplot(1,3,1)
    ax2=fg.add_subplot(1,3,2)
    ax3=fg.add_subplot(1,3,3)
    ax1.plot(veps,vsig);ax1.plot([veps[km]],[vsig[km]],marker='*',markersize=10)
    ax1.set_title(r'$\sigma = H(\epsilon)$')
    ax2.plot(veps,vhrat);ax2.plot([veps[km]],[vhrat[km]],marker='*',markersize=10)
    ax2.set_title(r"$h_r(\epsilon)=H'(0)/H'(\epsilon)$")
    ax3.plot(veps,vdhrat);ax3.plot([veps[km]],[vdhrat[km]],marker='*',markersize=10)
    ax3.set_title(r"$dHRATIO: h_r'(\epsilon)$")
    #uax.plt.show();uax.sysexit() 
    plt.show();sysexit()     
######-----------------------------------------------------------------------------------------------------
#
#
def abqHillRatios(rVals):
    r0,r45,r90=rVals
    R11=1
    R22=msqrt((r90*(r0+1))/(r0*(r90+1)))
    R33=msqrt((r90*(r0+1))/(r0+r90))
    R12=msqrt((3*r90*(r0+1))/((2*r45+1)*(r0+r90)))
    rr=5
    return [round(R11,rr),round(R22,rr),round(R33,rr),round(R12,rr),round(1.0,rr),round(1.0,rr)]


def strainsHill_old(rVals,hardParam,vangle,hratio=HRATIO):
    r0,r45,r90=rVals
    G=1/(1+r0)
    H=r0/(1+r0)
    F=H/r90
    N2=(2*r45+1)*(F+G)
    GH,FH=G+H,F+H
    e0,nn=hardParam[1:]
    e0n=e0**(1+nn)
    ##small hratio for testing
    #hratio=20 ## max strain eMax corresponds to H'(eMax)=H'(0)/hratio
    #hratio=HRATIO
    rr=5
    erd=e0*(hratio**(1/(1-nn))-1.0)
    ##vangle=[0,15, 30, 45, 60, 75, 90]
    vemax=[round(erd,rr)]
    rad=mpi/180.0
    for angle in vangle[1:]:
        tt=angle*rad
        ct,st=mcos(tt),msin(tt)
        ct,st=ct*ct,st*st
        ft=msqrt(GH*ct*ct+FH*st*st+(N2-2*H)*ct*st)
        em=(ft*((e0+erd)**(1+nn)-e0n)+e0n)**(1.0/(1+nn))-e0
        vemax.append(round(em,rr))
    return vemax


def strainsHill(hLaw,hardParam,vangle,hratio=HRATIO):
    #F,G,H,N2=yParam
    if(hLaw=='Voce'):
        em=mlog(hratio)/hardParam[2]
    elif(hLaw=='Swift'):
        em=hardParam[1]*(hratio**(1.0/(1.0-hardParam[2]))-1.0)    
    return [round(em,5) for angle in vangle]
    
def vStrainsHill(hLaw,hardParam,vangle,vcf=[],hratio=HRATIO):
    if(hLaw=='Swift'):
        e0,nn=hardParam[1:]
        e0n=e0**(1+nn)
        erd=e0*(hratio**(1.0/(1.0-nn))-1.0)
        #if(0.0 in vangle):????
        #    vemax=[erd]?????
        if(vcf):
            for angle in vangle[1:]:
                tt=angle*rad
                ct,st=mcos(tt),msin(tt)
                ft=fYF(ct*ct,st*st,ct*st,vcf)
                em=(ft*((e0+erd)**(1+nn)-e0n)+e0n)**(1.0/(1+nn))-e0
                vemax.append(em)
        else:
            r0,r45,r90=aMAT['rVals']
            G=1/(1+r0)
            H=r0/(1+r0)
            F=H/r90
            N2=(2*r45+1)*(F+G)
            GH,FH=G+H,F+H        
            for angle in vangle[1:]:
                tt=angle*rad
                ct,st=mcos(tt),msin(tt)
                ct,st=ct*ct,st*st
                ft=msqrt(GH*ct*ct+FH*st*st+(N2-2*H)*ct*st)
                em=(ft*((e0+erd)**(1+nn)-e0n)+e0n)**(1.0/(1+nn))-e0
                vemax.append(em)
        return vemax                
    elif(hLaw=='Voce'):
        a,b,c=hardParam[0:]
        erd=mlog(hratio)/c
        vemax=[erd]
        if(vcf):
            for angle in vangle[1:]:
                tt=angle*rad
                ct,st=mcos(tt),msin(tt)
                ft=fYF(ct*ct,st*st,ct*st,vcf)
                KT=ft*(a*erd-b*(1.0-1.0/mexp(c*erd))/c)+b/c
                em=0.0
                for kIter in range(100):
                    mcm=mexp(c*em)
                    phi=a*em+b/(c*mcm)-KT
                    if(abs(phi)<1.0e-7):break
                    dphi=a-b/mcm
                    em-=phi/dphi
                vemax.append(em)
        else:
            r0,r45,r90=aMAT['rVals']
            G=1/(1+r0)
            H=r0/(1+r0)
            F=H/r90
            N2=(2*r45+1)*(F+G)
            GH,FH=G+H,F+H        
            for angle in vangle[1:]:
                tt=angle*rad
                ct,st=mcos(tt),msin(tt)
                ct,st=ct*ct,st*st
                ft=msqrt(GH*ct*ct+FH*st*st+(N2-2*H)*ct*st)
                KT=ft*(a*erd-b*(1.0-1.0/mexp(c*erd))/c)+b/c
                em=0.0
                for kIter in range(100):
                    mcm=mexp(c*em)
                    phi=a*em+b/(c*mcm)-KT
                    if(abs(phi)<1.0e-7):break
                    dphi=a-b/mcm
                    em-=phi/dphi
                vemax.append(em)
        return vemax                
    else:printMsg("vStrainsHill: hardening must be in ['Swift','Voce']");sysexit()    
                    

def hardFunc(hLaw,hParam,epMax):
    deps,eps=5.0e-4,0.0
    tb,nn=[],int(epMax/deps)+2
    if(hLaw=='Swift'):###Swift: sigma=K*(eps0+eps)^n
        for k in range(nn):
            tb.append((round(hParam[0]*(hParam[1]+eps)**hParam[2],4),round(eps,5)))
            eps+=deps
    elif(hLaw=='Voce'):
        for k in range(nn):
            tb.append((round(hParam[0]-hParam[1]*mexp(-hParam[2]*eps),4),round(eps,5)))
            eps+=deps    
    return tb

#def hardFunc2(hParam,logStrain=0.4):                
#    return round(hParam[0]*(hParam[1]+logStrain)**hParam[2],3)
    
def paramHill(rVals):
    r0,r45,r90=rVals
    G=1/(1+r0)
    H=r0/(1+r0)
    F=H/r90
    N2=(2*r45+1)*(F+G)
    return [F,G,H,N2]
    
##reset everything
##Mdb()

def genABQmodel():
    session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
    mmd=mdb.Model(name=mdbModelName)
    ######
    LL2,WW2=0.5*LL,0.5*WW
    sRR=3*WW2  ;sRR=6*WW2
    ptheta=30*mpi/180.0
    cptheta,sptheta=mcos(ptheta),msin(ptheta)
    dRRx,dRRy=sRR*cptheta,sRR*sptheta
    sLLmax,sWWmax=LL+2.0*sRR,WW2+sRR
    sk=mmd.ConstrainedSketch(name='pProfile', sheetSize=2*(LL+WW))
    skP=[(-LL,-WW2),(LL,-WW2,LL,-sWWmax),(LL+sRR,-sWWmax),(sLLmax,-sWWmax),(sLLmax,sWWmax),(LL+sRR,sWWmax,LL,sWWmax),
    (LL,WW2),(-LL,WW2,-LL,sWWmax),(-LL-sRR,sWWmax),(-sLLmax,sWWmax),(-sLLmax,-sWWmax),(-LL-sRR,-sWWmax,-LL,-sWWmax),
    (-LL,-WW2)]      
    skDrw=['pp','cc','pp','pp','pp','cc','pp','cc','pp','pp','pp','cc']     
    for ks in range(len(skDrw)):
        drw=skDrw[ks]
        if(drw=='pp'):
            sk.Line(point1=(skP[ks][0],skP[ks][1]), point2=(skP[ks+1][0],skP[ks+1][1]))
        if(drw=='cc'):    
            sk.ArcByCenterEnds(center=(skP[ks][2],skP[ks][3]), 
            point1=(skP[ks][0],skP[ks][1]), point2=(skP[ks+1][0],skP[ks+1][1]), direction=CLOCKWISE)
    ###sk.Line(point1=skP[-1],point2=skP[0])

    pPart=mmd.Part(dimensionality=THREE_D, name=partName, type=DEFORMABLE_BODY)
    pPart.BaseShell(sketch=sk)
    del mmd.sketches['pProfile']
    session.viewports['Viewport: 1'].setValues(displayedObject=pPart)
    #
    W75=0.75*WW2
    setWhole=pPart.Set(faces=pPart.faces.findAt(((LL,zro,zro),)), name='Set-1-Part-1')
    #edgeLeftA=pPart.edges.findAt(((zro,W75,zro),))
    #edgeLeftB=pPart.edges.findAt(((zro,-W75,zro),))
    #setEdgeLeft=pPart.Set(edges=(edgeLeftA,edgeLeftB,),name='EdgeLeft')
    setEdgeLeft=pPart.Set(edges=pPart.edges.findAt(((-sLLmax,zro,zro),)),name='EdgeLeft')
    edgeBottomA=pPart.edges.findAt(((-LL-1.5*sRR,-sWWmax,zro),))
    edgeBottomB=pPart.edges.findAt(((LL+1.5*sRR,-sWWmax,zro),))
    setEdgeBottom=pPart.Set(edges=(edgeBottomA,edgeBottomB,),name='EdgeBottom')
    edgeRight=pPart.edges.findAt(((sLLmax,zro,zro),)) ##
    surfRight=pPart.Surface(side1Edges=edgeRight, name='Surf-Load') ##this is a Surface object used for load BC's (pressure)
    setEdgeRight=pPart.Set(edges=edgeRight,name='EdgeRight') ##this is used for displacement BC's
    ##setVertex=pPart.Set(name='VertexLL', vertices=pPart.vertices.findAt(((orgn,orgn,zro),)))
    matCSYS=pPart.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='csys-1-Material', origin=(zro,zro,zro), 
            point1=(LL,zro,zro), point2=(zro,WW,zro))    
    ##
    #mmat=mmd.Material(description=aMAT['descrp'], name=aMAT['name'])
    ###
    aa=mmd.rootAssembly
    aa.DatumCsysByDefault(CARTESIAN)
    aa.Instance(dependent=OFF, name=partInstanceName, part=pPart)
    ipPart=aa.instances[partInstanceName] ##Note: this ipPart points to the instance of pPart
    ####
    mmd.StaticStep(description='Apply tensile load', initialInc=dtInitial, 
        matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=dtInitial/100.0,maxInc=dtMax, maxNumInc=maxNInc, 
        name='Step-1', nlgeom=ON, previous='Initial')  
    #mmd.steps['Step-1'].setValues(adaptiveDampingRatio=0.001, continueDampingFactors=False, stabilizationMagnitude=0.0002, 
    #    stabilizationMethod=DISSIPATED_ENERGY_FRACTION)
    mmd.FieldOutputRequest(name='F-Output-1',createStepName='Step-1')
    for kk in mmd.historyOutputRequests.keys():
        del mmd.historyOutputRequests[kk]        
    ###
    ###mmd.XsymmBC(name='BC-1-SymmX', createStepName='Initial', region=ipPart.sets['EdgeLeft'], localCsys=None)
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-1-FixXZ',region=ipPart.sets['EdgeLeft'], 
        u1=SET, u2=UNSET, u3=SET, ur1=UNSET, ur2=UNSET, ur3=UNSET)    
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-3-FixY',region=ipPart.sets['EdgeBottom'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Step-1',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='BC-4-MoveX', region=ipPart.sets['EdgeRight'], 
        u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)        
    #
    fac = ipPart.faces.findAt(coordinates=(zro, zro, zro))
    aa.PartitionFaceByShortestPath(point1=(-sLLmax,zro,zro), point2=(sLLmax,zro,zro), faces=fac)
    fac = ipPart.faces.findAt(coordinates=(zro, W75, zro))
    aa.PartitionFaceByShortestPath(point1=(-LL,WW2,zro), point2=(-LL,-WW2,zro), faces=fac)
    fac = ipPart.faces.findAt(coordinates=(zro, -W75, zro))
    aa.PartitionFaceByShortestPath(point1=(-LL,WW2,zro), point2=(-LL,-WW2,zro), faces=fac)
    fac = ipPart.faces.findAt(coordinates=(zro, W75, zro))
    aa.PartitionFaceByShortestPath(point1=(LL,WW2,zro), point2=(LL,-WW2,zro), faces=fac)
    fac = ipPart.faces.findAt(coordinates=(zro, -W75, zro))
    aa.PartitionFaceByShortestPath(point1=(LL,WW2,zro), point2=(LL,-WW2,zro), faces=fac)
    #
    xDatum=LL+(2*sRR+WW2)*((1.0-sptheta)/cptheta)
    fac = ipPart.faces.findAt(coordinates=(LL+sRR,W75, zro))
    vedg=ipPart.edges
    edgeZ=vedg.getByBoundingBox(xMin=LL-gTol,xMax=LL+sRR+gTol,yMin=WW2-gTol,yMax=sWWmax+gTol,zMin=-gTol,zMax=gTol)[0]
    aa.PartitionFaceByCurvedPathEdgePoints(face=fac, 
        edge1=vedg.findAt(coordinates=(LL+sRR,zro,zro)), point1=(xDatum,zro,zro), 
        edge2=edgeZ, point2=(LL+dRRx,sWWmax-dRRy,zro))
    fac = ipPart.faces.findAt(coordinates=(LL+sRR,-W75, zro))
    vedg=ipPart.edges
    edgeZ=vedg.getByBoundingBox(xMin=LL-gTol,xMax=LL+sRR+gTol,yMax=-WW2+gTol,yMin=-sWWmax-gTol,zMin=-gTol,zMax=gTol)[0]
    aa.PartitionFaceByCurvedPathEdgePoints(face=fac, 
        edge1=vedg.findAt(coordinates=(LL+sRR/2.0,zro,zro)), point1=(xDatum,zro,zro), 
        edge2=edgeZ, point2=(LL+dRRx,-sWWmax+dRRy,zro))  
    fac = ipPart.faces.findAt(coordinates=(-LL-sRR,W75, zro))
    vedg=ipPart.edges
    edgeZ=vedg.getByBoundingBox(xMax=-LL+gTol,xMin=-LL-sRR-gTol,yMin=WW2-gTol,yMax=sWWmax+gTol,zMin=-gTol,zMax=gTol)[0]
    aa.PartitionFaceByCurvedPathEdgePoints(face=fac, 
        edge1=vedg.findAt(coordinates=(-LL-sRR,zro,zro)), point1=(-xDatum,zro,zro), 
        edge2=edgeZ, point2=(-LL-dRRx,sWWmax-dRRy,zro))
    fac = ipPart.faces.findAt(coordinates=(-LL-sRR/2.0,-W75, zro))
    vedg=ipPart.edges
    edgeZ=vedg.getByBoundingBox(xMax=-LL+gTol,xMin=-LL-sRR-gTol,yMax=-WW2+gTol,yMin=-sWWmax-gTol,zMin=-gTol,zMax=gTol)[0]
    aa.PartitionFaceByCurvedPathEdgePoints(face=fac, 
        edge1=vedg.findAt(coordinates=(-LL-sRR/2.0,zro,zro)), point1=(-xDatum,zro,zro), 
        edge2=edgeZ, point2=(-LL-dRRx,-sWWmax+dRRy,zro))  
    #    
    fac=ipPart.faces.findAt(((zro,W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((zro,-W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((LL+sRR/2.0,W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((LL+sRR/2.0,-W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((-LL-sRR/2.0,W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((-LL-sRR/2.0,-W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=ipPart.faces.findAt(((LL+1.5*sRR,W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac)
    fac=ipPart.faces.findAt(((LL+1.5*sRR,-W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac)
    fac=ipPart.faces.findAt(((-LL-1.5*sRR,W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac)
    fac=ipPart.faces.findAt(((-LL-1.5*sRR,-W75,zro),))
    aa.setMeshControls(elemShape=QUAD, regions=fac)
    adEL=0.0
    nY=nElemMultiplier*int(WW2+adEL)
    nX=nElemMultiplier*int(LL*(1.0+adEL/WW2))
    #n2X=int(nX*(sRR*(0.5*mpi-ptheta))/LL)
    n2X=int(nX*((xDatum-LL)/LL))
    nERatioX,nERatioY=float(nX)/float(nElemMultiplier),float(nY)/float(nElemMultiplier)
    n3X,n2Y=int(nERatioX*((sLLmax-xDatum)/LL)+7),int(nERatioY*(WW2+sRR)/WW2+7)
    vedg=ipPart.edges
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL,W75,zro),)), number=nY)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL,-W75,zro),)), number=nY)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL,W75,zro),)), number=nY)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL,-W75,zro),)), number=nY)
    ##aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt((skC,skD)), number=nY)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((zro,zro,zro),)), number=2*nX)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL+0.5*sRR,zro,zro),)), number=n2X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL-0.5*sRR,zro,zro),)), number=n2X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL+1.5*sRR,-sWWmax,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL+1.5*sRR,zro,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((LL+1.5*sRR,sWWmax,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL-1.5*sRR,-sWWmax,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL-1.5*sRR,zro,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-LL-1.5*sRR,sWWmax,zro),)), number=n3X)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((sLLmax,W75,zro),)), number=n2Y)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((sLLmax,-W75,zro),)), number=n2Y)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-sLLmax,W75,zro),)), number=n2Y)
    aa.seedEdgeByNumber(constraint=FINER, edges=ipPart.edges.findAt(((-sLLmax,-W75,zro),)), number=n2Y)
    dX,dY=float(LL)/float(nX),float(WW2)/float(nY)
    #
    if(shell_S4R):
        aa.setElementType(elemTypes=(mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=ON, hourglassControl=ENHANCED),),
                          regions=ipPart.sets['Set-1-Part-1']) 
    else:
        aa.setElementType(elemTypes=(mesh.ElemType(elemCode=S8R, elemLibrary=STANDARD,secondOrderAccuracy=ON),),regions=ipPart.sets['Set-1-Part-1'])
    aa.generateMesh(regions=ipPart.faces.findAt(((zro,W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((zro,-W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((LL+sRR/2.0,W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((LL+sRR/2.0,-W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((-LL-sRR/2.0,W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((-LL-sRR/2.0,-W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((LL+1.5*sRR,W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((LL+1.5*sRR,-W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((-LL-1.5*sRR,W75,zro),)))
    aa.generateMesh(regions=ipPart.faces.findAt(((-LL-1.5*sRR,-W75,zro),)))
    ####
    mdb.Job(atTime=None, contactPrint=OFF, description=JOB_DESCRIPTION, echoPrint=OFF, explicitPrecision=SINGLE, 
        getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, 
        model=mmd, modelPrint=OFF, multiprocessingMode=DEFAULT, 
        name=JOB_NAME, nodalOutputPrecision=FULL, numCpus=1, numGPUs=0, 
        queue=None, resultsFormat=ODB, scratch=thisPATH, 
        type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    ##if(UMAT):
    ##    mdb.jobs[JOB_NAME].setValues(userSubroutine=UMATpath)
    mdb.saveAs(modelName)
    #
    epsX,epsY=dX/3.0,dY/3.0
    #epsX,epsY=dX/2.0,dY/2.0
    nodX,nodY=(0.75*nX)*dX,zro
    ndA=ipPart.nodes.getByBoundingBox(xMin=nodX-epsX,xMax=nodX+epsX,yMin=nodY-epsY,yMax=nodY+epsY,zMin=-gTol,zMax=gTol)
    ndB=ipPart.nodes.getByBoundingBox(xMin=-nodX-epsX,xMax=-nodX+epsX,yMin=nodY-epsY,yMax=nodY+epsY,zMin=-gTol,zMax=gTol)
    nodX,nodY=zro,WW2
    ndC=ipPart.nodes.getByBoundingBox(xMin=nodX-epsX,xMax=nodX+epsX,yMin=nodY-epsY,yMax=nodY+epsY,zMin=-gTol,zMax=gTol)
    ndD=ipPart.nodes.getByBoundingBox(xMin=nodX-epsX,xMax=nodX+epsX,yMin=-nodY-epsY,yMax=-nodY+epsY,zMin=-gTol,zMax=gTol)
    #print ndA[0].coordinates
    #print ndB[0].coordinates
    eX,eY=dX/2.0,dY/2.0
    #print dX,dY,eX,eY
    elm=ipPart.elements.getByBoundingBox(xMin=eX-dX-epsX,xMax=eX+dX+epsX,yMin=eY-dY-epsY,yMax=eY+dY+epsY,zMin=-gTol,zMax=gTol)
    #print ipPart.nodes[elm[0].connectivity[0]].coordinates
    lenA,lenB,lenC,lenD,lenE=len(ndA),len(ndB),len(ndC),len(ndD),len(elm)
    if((lenA==1) and (lenB==1) and (lenC==1) and (lenD==1) and (lenE==1)):
        ndAlabel,ndBlabel,ndClabel,ndDlabel,elmLabel=ndA[0].label,ndB[0].label,ndC[0].label,ndD[0].label,elm[0].label
    else:
        printMsg('------- lenA, lenB, lenC, lenD, lenE: {}, {}, {}, {}, {}'.format(lenA,lenB,lenC,lenD,lenE))
        printMsg('lenA, lenB, lenC, lenD and lenE must = 1')
        printMsg('Check your dimensions, sampling node/element, seeds and/or dX, dY definitions...')
        printMsg('Calculations aborted');sysexit()
    return mdb,matCSYS.id,ndAlabel,ndBlabel,ndClabel,ndDlabel,elmLabel
    
    
####-------------------------------------------------------------------------------------------------------------------------------
####def runTests(mdb,datumID,ndAlabel,ndBlabel,ndClabel,ndDlabel,elmLabel):
def runTests(nCPUs=1):
    oschdir(thisPATH)
    mdb,datumID,ndAlabel,ndBlabel,ndClabel,ndDlabel,elmLabel=genABQmodel()
    ##errStatusFile=figDirPlot+'_000_allJobsStatus.txt'
    ##errf=open(errStatusFile,'w');errf.close()
    for ttest in vvTests:
        vAngles=ttest['angles']
        ff=figDirPlot+ttest['dfile'].strip('.txt')
        vStrains=strainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles,hratio=HRATIO)
        vDeltaL=[LL*(mexp(eep)-1.0) for eep in vStrains]
        #printMsg(str(vStrains));printMsg(str(vDeltaL));sysexit()
        mmd=mdb.models[mdbModelName]
        mmat=mmd.Material(description=aMAT['descrp'], name=aMAT['name'])
        polyParam=[]
        if(ttest['PolyN']):
            try:
                zff=open(figDirData+ttest['dfile'],'r')
                for line in zff:
                    polyParam.append(float(line.strip()))
                zff.close()    
            except IOError:
                printMsg('Input file for Poly_UMAT not found: '+ttest['dfile'])
                printMsg('Calculations aborted');sysexit()
            mdb.jobs[JOB_NAME].setValues(userSubroutine=ttest['ufile'])    
            mmat.UserMaterial(mechanicalConstants=tuple(polyParam))
            mmat.Depvar(n=1)
            ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles,vcf=polyParam[5:])
        elif(ttest['UMAT']):
            mdb.jobs[JOB_NAME].setValues(userSubroutine=ttest['ufile'])
            hardP,matP=aMAT['hardParam'],paramHill(aMAT['rVals'])
            #tmp=(aMAT['E'],aMAT['nuPoisson'],hardP[0],hardP[1],hardP[2],matP[0],matP[1],matP[2],matP[3])
            #tmp=(len(tmp)*'{:.4f}, ').format(*tmp);printMsg(tmp);sysexit()
            mmat.UserMaterial(mechanicalConstants=(aMAT['E'],aMAT['nuPoisson'], 
            hardP[0],hardP[1],hardP[2],matP[0],matP[1],matP[2],matP[3]))
            mmat.Depvar(n=1)
            ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles)
        else:
            mmat.Elastic(table=((aMAT['E'], aMAT['nuPoisson']), ))
            epMax=max(vStrains)+0.001 ###max strain for *MATERIAL hardening curve 
            mmat.Plastic(table=hardFunc(aMAT['hardLaw'],aMAT['hardParam'],epMax))
            mmat.plastic.Potential(table=(abqHillRatios(aMAT['rVals']),))
            ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles)            
        #else:
        #    pass
    ###
        mmd.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
        integrationRule=SIMPSON, material=aMAT['name'], name='Section-1', 
        numIntPts=5, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
        thickness=aMAT['hThick'], thicknessField='', thicknessModulus=None, 
        thicknessType=UNIFORM, useDensity=OFF)
        aa=mmd.rootAssembly
        #pPart=mdb.parts[partName]
        pPart=aa.instances[partInstanceName].part
        mmd.sections['Section-1'].TransverseShearShell(k11=transvK,k12=0.0, k22=transvK)
        pPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, region=pPart.sets['Set-1-Part-1'], 
        sectionName='Section-1', thicknessAssignment=FROM_SECTION)
        if(ttest['UMAT']):    
            mmd.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'U', 'SDV'))
        else:
            mmd.fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE', 'U', 'PE', 'PEEQ'))  
        aa.regenerate()
        ##vDeltaL=[LL*(mexp(eep)-1.0) for eep in vStrains]
        vFileStress,vFileNode,vFileJob=[],[],[]
        casePoly=ttest['PolyN']
        mdb.jobs[JOB_NAME].setValues(numCpus=nCPUs)
    ##############################Loop on angles -------------------------------------------------------------
        for idxAngle in range(len(vAngles)):
            ANGLE_FROM_RD=vAngles[idxAngle]        
            if ospath.exists(JOB_NAME+'.lck'):
                osremove(JOB_NAME+'.lck')
            strAngle=str(int(ANGLE_FROM_RD))
            fileStress=ff+'_ES_DATA_'+strAngle+'.txt';vFileStress.append(fileStress)
            fileNode=ff+'_Node_'+strAngle+'.txt';vFileNode.append(fileNode)
            fileJob=ff+'_JOB_'+strAngle+'.txt';vFileJob.append(fileJob)
            fileScreen=ff+'_Screen_'+strAngle+'.png'
            printMsg('---------- Testing with {} at angle from RD: {}'.format(ttest['dfile'], ANGLE_FROM_RD))
            #printMsg('---------------- Max increment: dtMax = {}'.format(dtMax))
            #if(UMAT):
            #    printMsg('---------------- Testing with UMAT: '+UMATname)
            pPart.MaterialOrientation(additionalRotationField='', additionalRotationType=ROTATION_ANGLE, 
                angle=ANGLE_FROM_RD, axis=AXIS_3, fieldName='', localCsys=pPart.datums[datumID], 
                orientationType=SYSTEM, region=pPart.sets['Set-1-Part-1'])
            ###printMsg(vDeltaL[idxAngle]);sysexit()    
            mmd.boundaryConditions['BC-4-MoveX'].setValues(u1=vDeltaL[idxAngle]+12.0)
            #if(loadType):
            #    mmd.DisplacementBC(amplitude=UNSET, createStepName='Step-1',distributionType=UNIFORM, fieldName='', 
            #                       fixed=OFF, localCsys=None, name='BC-4-MoveX', region=ipPart.sets['EdgeRight'], 
            #                       u1=DeltaL, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET)    
            #else:
            #    LOADpressure=-((WW*hhh)/(WW+2.0*sRR))*hardFunc2(aMAT['hardParam'],logStrain=vStrains[idxAngle]) 
            #    mmd.ShellEdgeLoad(createStepName='Step-1',distributionType=UNIFORM,field='',localCsys=None, 
            #                      magnitude=LOADpressure,name='Load-X',region=ipPart.surfaces['Surf-Load'],
            #                      resultant=ON)        
            aa.regenerate()
            if(caeOnly):
                cwd=osgetcwd()
                mdb.jobs[JOB_NAME].setValues(scratch=cwd)
                #printMsg('ufile:    '+ttest['ufile'])
                #mdb.jobs[JOB_NAME].setValues(userSubroutine=cwd+ttest['ufile'].strip(thisPATH))
                mdb.saveAs(modelName);printMsg("CAE file created: "+modelName);sysexit()
            #
            ##monitorManager.addMessageCallback(jobName=JOB_NAME,messageType=ANY_MESSAGE_TYPE,callback=monitMessage,
            ##userData=[errStatusFile,ttest['dfile'].strip('.txt')+'__'+strAngle])
            printMsg('Job submitted.........')            
            try:   
                mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
            except:
                pass    
            mdb.jobs[JOB_NAME].waitForCompletion()
            ##monitorManager.removeMessageCallback(jobName=JOB_NAME,messageType=ANY_MESSAGE_TYPE,callback=monitMessage,
            ##userData=[errStatusFile,ttest['dfile'].strip('.txt')+'__'+strAngle]) 
            printMsg('------------------------------ Opening Odb.....')
            odb=openOdb(path=JOB_NAME+'.odb')
            #podb=odb.rootAssembly.parts['PART-1-1']
            podb=odb.rootAssembly.instances['PART-1-1'] 
            velm=podb.elements[elmLabel-1]
            vnds=podb.nodes
            vndsA,vndsB,vndsC,vndsD=vnds[ndAlabel-1],vnds[ndBlabel-1],vnds[ndClabel-1],vnds[ndDlabel-1]  
            printMsg('Exporting stress-strain data for mesh element: '+str(elmLabel))
            printMsg('Exporting node data for node(s): '+str(ndAlabel)+', '+str(ndBlabel)+', '+str(ndClabel)+', '+str(ndDlabel))
            ff2=open(JOB_NAME+'.dat','r')
            wTime=0 ##Wall Clock 
            for line in ff2:
                if('WALLCLOCK TIME' in line):
                    wTime+=int(line.strip().split('=')[1].strip())
                    printMsg(line.strip())
            ff2.close()
            ff2=open(fileJob,'w')
            ff2.write('Angle | Wall Clock (sec)\n')
            ff2.write(strAngle+' | '+str(wTime)+'\n')
            ff2.close()
            #vStress,vStrain=[],[]
            ss11,ss22,ss33,ss12=0.0,0.0,0.0,0.0
            bss11,bss22,bss33,bss12=0.0,0.0,0.0,0.0
            ee11,ee22,ee33,ee12=0.0,0.0,0.0,0.0
            bee11,bee22,bee33,bee12=0.0,0.0,0.0,0.0
            ff1=open(fileStress,'w')
            ff1.write('Test angle: '+strAngle+'\n')
            ff1.write('Mesh element_'+str(elmLabel)+' : S11 | S22 | S33 | S12 | E11 | E22 | E33 | 2*E12 | EPbar \n')
            lineOut=4*"{:.3f} | " + 5*"{:.4f} |" + "\n"  
            rrS,rrE=3,4
            lineOut2=3*"{:.4f} | " + "\n"
            ff2=open(fileNode,'w')
            ff2.write('Test angle: '+strAngle+'\n')
            ff2.write('Node_A_{} coordinates: initialX | initialY | initialZ \n'.format(str(ndAlabel)))
            ff2.write(lineOut2.format(*vndsA.coordinates))
            ff2.write('Node_B_{} coordinates: initialX | initialY | initialZ \n'.format(str(ndBlabel)))
            ff2.write(lineOut2.format(*vndsB.coordinates))
            ff2.write('Node_C_{} coordinates: initialX | initialY | initialZ \n'.format(str(ndClabel)))
            ff2.write(lineOut2.format(*vndsC.coordinates))
            ff2.write('Node_D_{} coordinates: initialX | initialY | initialZ \n'.format(str(ndDlabel)))
            ff2.write(lineOut2.format(*vndsD.coordinates))
            ff2.write('  UA1  |  UA2  |  UA3 | UB1 | UB2 | UB3 | UC1 | UC2 | UC3 | UD1 | UD2 | UD3 \n')
            lineOut2=12*"{:.4f} | " + "\n"
            ct,st=mcos(rad*ANGLE_FROM_RD),msin(rad*ANGLE_FROM_RD)
            ctst,ct2,st2=ct*st,ct*ct,st*st
            zprec=odb.steps['Step-1'].frames[0].fieldOutputs['S'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values[0].precision
            if(zprec==SINGLE_PRECISION):
                for kframe in odb.steps['Step-1'].frames:
                    ss=kframe.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    ss11,ss22,ss33,ss12=0.0,0.0,0.0,0.0
                    try:
                        for itm in ss:
                            ss11+=itm.data[0]
                            ss22+=itm.data[1]
                            ss33+=itm.data[2]
                            ss12+=itm.data[3]
                    except:
                        for itm in ss:
                            ss11+=itm.dataDouble[0]
                            ss22+=itm.dataDouble[1]
                            ss33+=itm.dataDouble[2]
                            ss12+=itm.dataDouble[3]                
                    nss=len(ss)    
                    ss11,ss22,ss33,ss12=ss11/nss,ss22/nss,ss33/nss,ss12/nss ##components wrt local (material) frame 
                    bss11=round(ss11*ct2+ss22*st2-2.0*ss12*ctst,rrS)###components wrt global (spatial) frame 
                    bss22=round(ss11*st2+ss22*ct2+2.0*ss12*ctst,rrS)
                    bss12=round((ss11-ss22)*ctst+ss12*(ct2-st2),rrS) 
                    #if(UMAT):    
                    #    ss=kframe.fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    #else:
                    #    ss=kframe.fieldOutputs['PE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    ss=kframe.fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values    
                    ee11,ee22,ee33,ee12=0.0,0.0,0.0,0.0
                    try:
                        for itm in ss:
                            ee11+=itm.data[0]
                            ee22+=itm.data[1]
                            ee33+=itm.data[2]
                            ee12+=itm.data[3]
                    except:
                        for itm in ss:
                            ee11+=itm.dataDouble[0]
                            ee22+=itm.dataDouble[1]
                            ee33+=itm.dataDouble[2]
                            ee12+=itm.dataDouble[3]                    
                    #nss=len(ss)    
                    ee11,ee22,ee33,ee12=ee11/nss,ee22/nss,ee33/nss,ee12/nss ##components wrt local (material) frame 
                    bee11=round(ee11*ct2+ee22*st2-ee12*ctst,rrE)###components wrt global (spatial) frame 
                    bee22=round(ee11*st2+ee22*ct2+ee12*ctst,rrE)
                    bee12=round(2.0*(ee11-ee22)*ctst+ee12*(ct2-st2),rrE) ##Abaqus reports engineering shear strain
                    if(ttest['UMAT']):    
                        ss=kframe.fieldOutputs['SDV1'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    else:
                        ss=kframe.fieldOutputs['PEEQ'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values    
                    ee11=0.0
                    for itm in ss:
                        ee11+=itm.data
                    ee11/=nss
                    ff1.write(lineOut.format(bss11,bss22,bss33,bss12,bee11,bee22,bee33,bee12,ee11))
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsA).values
                    dUA=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsB).values
                    dUB=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsC).values
                    dUC=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsD).values
                    dUD=ss[0].dataDouble    
                    ff2.write(lineOut2.format(dUA[0],dUA[1],dUA[2],dUB[0],dUB[1],dUB[2],dUC[0],dUC[1],dUC[2],dUD[0],dUD[1],dUD[2]))
            elif(zprec==DOUBLE_PRECISION):
                for kframe in odb.steps['Step-1'].frames:
                    ss=kframe.fieldOutputs['S'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    ss11,ss22,ss33,ss12=0.0,0.0,0.0,0.0
                    for itm in ss:
                        ss11+=itm.dataDouble[0]
                        ss22+=itm.dataDouble[1]
                        ss33+=itm.dataDouble[2]
                        ss12+=itm.dataDouble[3]
                    nss=len(ss)    
                    ss11,ss22,ss33,ss12=ss11/nss,ss22/nss,ss33/nss,ss12/nss ##components wrt local (material) frame 
                    bss11=round(ss11*ct2+ss22*st2-2.0*ss12*ctst,rrS)###components wrt global (spatial) frame 
                    bss22=round(ss11*st2+ss22*ct2+2.0*ss12*ctst,rrS)
                    bss12=round((ss11-ss22)*ctst+ss12*(ct2-st2),rrS) 
                    #if(UMAT):    
                    #    ss=kframe.fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    #else:
                    #    ss=kframe.fieldOutputs['PE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    ss=kframe.fieldOutputs['LE'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values    
                    ee11,ee22,ee33,ee12=0.0,0.0,0.0,0.0
                    for itm in ss:
                        ee11+=itm.dataDouble[0]
                        ee22+=itm.dataDouble[1]
                        ee33+=itm.dataDouble[2]
                        ee12+=itm.dataDouble[3]
                    #nss=len(ss)    
                    ee11,ee22,ee33,ee12=ee11/nss,ee22/nss,ee33/nss,ee12/nss ##components wrt local (material) frame 
                    bee11=round(ee11*ct2+ee22*st2-ee12*ctst,rrE)###components wrt global (spatial) frame 
                    bee22=round(ee11*st2+ee22*ct2+ee12*ctst,rrE)
                    bee12=round(2.0*(ee11-ee22)*ctst+ee12*(ct2-st2),rrE) ##Abaqus reports engineering shear strain
                    if(ttest['UMAT']):    
                        ss=kframe.fieldOutputs['SDV1'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values
                    else:
                        ss=kframe.fieldOutputs['PEEQ'].getSubset(position=INTEGRATION_POINT,elementType=elmType,region=velm).values    
                    ee11=0.0
                    for itm in ss:
                        ee11+=itm.data
                    ee11/=nss
                    ff1.write(lineOut.format(bss11,bss22,bss33,bss12,bee11,bee22,bee33,bee12,ee11))
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsA).values
                    dUA=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsB).values
                    dUB=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsC).values
                    dUC=ss[0].dataDouble
                    ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vndsD).values
                    dUD=ss[0].dataDouble    
                    ff2.write(lineOut2.format(dUA[0],dUA[1],dUA[2],dUB[0],dUB[1],dUB[2],dUC[0],dUC[1],dUC[2],dUD[0],dUD[1],dUD[2]))               
            ff1.close()
            ff2.close()
            ######session.viewports['Viewport: 1'].setValues(displayedObject=odb)
            #####odb.close()
            ##
            session.printOptions.setValues(vpDecorations=OFF)
            session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
            #session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(datumCoordSystems=ON)
            #session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(datumCoordSystems=OFF)
            session.printOptions.setValues(vpBackground=ON) 
            vv=session.viewports['Viewport: 1']
            vv.setValues(displayedObject=odb)
            vv.view.fitView(drawImmediately=True)
            vv.odbDisplay.display.setValues(plotState=(UNDEFORMED, CONTOURS_ON_DEF, ))
            vv.view.setViewpoint(viewVector=(0.0,0,1.0),cameraUpVector=(0,1,0),drawImmediately=True)
            vv.odbDisplay.superimposeOptions.setValues(visibleEdges=FEATURE)
            vv.odbDisplay.superimposeOptions.setValues(renderStyle=FILLED)
            vv.odbDisplay.superimposeOptions.setValues(deformedOffsetMode=NONUNIFORM)
            vv.odbDisplay.superimposeOptions.setValues(nonuniformOffset=(0, 0, 1.))
            #vv.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,'S11'))
            vv.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT,'Mises'))
            vv.view.fitView()
            vv.view.zoom(0.999)
            session.printToFile(fileName=fileScreen, format=PNG, canvasObjects=(vv,))
            odb.close()
            session.graphicsOptions.setValues(backgroundStyle=GRADIENT,backgroundColor='#1B2D46')
            printMsg('---------- '+ff+'/'+strAngle+' : DONE--------------------------------\n\n')
        if(pytPlot and False):
            plotOneTest(baseName=ff,vAngle=vAngles,vfStress=vFileStress,vfNode=vFileNode,vfJob=vFileJob,zvcf=polyParam,casePoly=casePoly)
        ###del mmd.fieldOutputRequests['F-Output-1']
        del mmd.sections['Section-1']
        del mmd.materials[aMAT['name']]        
    printMsg('------------ ALL Abaqus SIM TESTS: DONE ---------------------------------------------------')    
    
####----------------------------------------------------------------------------------------------------------------
def paramPlotHill(rVals):
    r0,r45,r90=rVals
    G=1/(1+r0)
    H=r0/(1+r0)
    F=H/r90
    N2=(2*r45+1)*(F+G)
    GH,FH=G+H,F+H
    return GH,FH,H,N2

def fHill(GH,FH,H,N2,theta):
    ct,st=znp.cos(theta),znp.sin(theta)
    ct,st=ct*ct,st*st
    return znp.sqrt(GH*ct*ct+FH*st*st+(N2-2*H)*ct*st)
    
def OBSrvalHill(GH,FH,H,N2):
    vtheta=znp.linspace(0,znp.pi/2,nPlotPoints)
    ct,st=znp.cos(vtheta),znp.sin(vtheta)
    ctst=ct*st
    ct,st=ct*ct,st*st
    dx=GH*ct-H*st
    dy=FH*st-H*ct
    dxy=N2*ctst
    return vtheta,(dxy*ctst-(dx*st+dy*ct))/(dx+dy)
    
def rvalHill2(theta,GH,FH,H,N2):
    ct,st=znp.cos(theta),znp.sin(theta)
    ctst=ct*st
    ct,st=ct*ct,st*st
    dx=GH*ct-H*st
    dy=FH*st-H*ct
    dxy=N2*ctst
    return (dxy*ctst-(dx*st+dy*ct))/(dx+dy)

def fYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        #KKo+=2
        #tau=abs(sxy)
        #yf=vcf[KKo]*tau       
        return (vcf[KKo+2]**(1.0/float(deg)))*abs(sxy)    
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        #kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[2:nCoeff+2]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]  
    yf=cc[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
        yf=bb+yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    #zyf=yf**odeg
    #zyf=att*zyf    
    return abs(tt)*(yf**(1.0/float(deg)))
    
def fGYF(sx,sy,sxy,vcf):
    deg=int(vcf[0])
    odeg=1.0/float(deg)
    nMon=int(deg/2)
    nCoeff=int(vcf[1])
    KKo=nCoeff
    KKo-=1
    msx,msy=abs(sx),abs(sy)
    if(msx+msy<gTol):
        KKo+=2
        tau=abs(sxy)
        yf=(vcf[KKo]**odeg)*tau 
        tau=sxy/tau ##( = +/- 1)
        ##gyf=[d/dsx,d/dsy,d/dsxy]     
        return yf,[0.0,0.0,tau]
    if(msx>=msy):
        tt=sx;rho=sy/sx;gamm=sxy/sx;gamm=gamm*gamm
        kcfA=2;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    else: 
        tt=sy;rho=sx/sy;gamm=sxy/sy;gamm=gamm*gamm
        kcfA=2+6*nCoeff;kcfB=kcfA+nCoeff 
        cc=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d1=vcf[kcfA:kcfB]
        kcfA=kcfB;kcfB+=nCoeff
        d2=vcf[kcfA:kcfB]
    #d3=znp.zeros(nMon)
    d3=[0.0 for pq in range(nMon)]    
    yf,g1yf,g2yf=cc[KKo],d1[KKo],d2[KKo]
    jj,jj2=nMon,nMon
    while(jj>0):
        KKo-=1
        bb=cc[KKo]
        d1bb,d2bb=d1[KKo],d2[KKo]
        jj-=1
        mm=deg-2*jj ##degree of multiplying polynomial in sx and sy
        while(mm>0):
            mm-=1
            KKo-=1
            bb=cc[KKo]+bb*rho
            d1bb=d1[KKo]+d1bb*rho
            d2bb=d2[KKo]+d2bb*rho           
        d3[jj-1]=jj*bb
        yf=bb+yf*gamm
        g1yf=d1bb+g1yf*gamm
        g2yf=d2bb+g2yf*gamm  
    #zyf=(vcf[1]*yf)**odeg  ##in case of coefficients normalization 
    zyf=yf**odeg
    yval,att,y2val=zyf/(deg*yf),abs(tt),float(deg-1)/zyf
    sgntt=(tt/att)*yval
    g1yf=sgntt*g1yf
    g2yf=sgntt*g2yf  
    jj2-=1    
    g3yf=float(nMon)*cc[-1]
    while(jj2>0):
        jj2-=1
        g3yf=d3[jj2]+g3yf*gamm
    tt=(2.0*sxy)/tt    
    g3yf=sgntt*tt*g3yf
    zyf=att*zyf    
    return zyf,[g1yf,g2yf,g3yf]



def nnfYF(vx,dictParam):###
    nnLayers=dictParam['nnLayers']
    nAllUnits=dictParam['nAllUnits']
    firstWeight=dictParam['firstWeight']
    vUnits=dictParam['vUnits']
    vDegree=dictParam['vDegree']
    vw=dictParam['vNNparam']
    ###
    ssum,Fval=0.0,0.0
    KKpreviousLayer=0
    KKcurrentLayer=vUnits[1]
    vAllUnits=znp.zeros(nAllUnits+1)
    vAllUnits[1:KKcurrentLayer+1]=vx[:]
    KKw=firstWeight
    for KKlayer in range(2,nnLayers+1): 
        for KKunit in range(1,vUnits[KKlayer]+1):
            ssum=0.0
            for JJ in range(1,vUnits[KKlayer-1]+1):
                ssum=ssum+vAllUnits[KKpreviousLayer+JJ]*vw[KKw]
                KKw=KKw+1  
            vAllUnits[KKcurrentLayer+KKunit]=ssum**vDegree[KKlayer]
        KKpreviousLayer=KKcurrentLayer
        KKcurrentLayer=KKcurrentLayer+vUnits[KKlayer]   
    Fval=0.0
    for JJ in range(1,vUnits[nnLayers]+1): ### Calculate output
        Fval=Fval+vAllUnits[KKpreviousLayer+JJ]*vw[KKw]
        KKw=KKw+1
    return  Fval**dictParam['homDegReciproc']

def nnfGYF(vx,dictParam):
    nnLayers=dictParam['nnLayers']
    nAllUnits=dictParam['nAllUnits']
    firstWeight=dictParam['firstWeight']
    vUnits=dictParam['vUnits']
    vDegree=dictParam['vDegree']
    vw=dictParam['vNNparam']
    KKpreviousLayer=0
    KKcurrentLayer=vUnits[1]
    vAllUnits=znp.zeros(nAllUnits+1)
    vAllUnits[1:KKcurrentLayer+1]=vx[:]
    gDegree=znp.zeros(nnLayers+1)
    gDegree[:]=vDegree[:]-1.0
    gAllUnits=znp.zeros((nAllUnits+1,4))
    gAllUnits[1,1],gAllUnits[2,2],gAllUnits[3,3]=1.0,1.0,1.0
    FGval=znp.zeros(3)
    KKw=firstWeight
    ssum,gsum=0.0,0.0
    for KKlayer in range(2,nnLayers+1): 
        for KKunit in range(1,vUnits[KKlayer]+1):
            ssum, FGval[:]=0.0,0.0
            for JJ in range(1,vUnits[KKlayer-1]+1):
                ssum=ssum+vAllUnits[KKpreviousLayer+JJ]*vw[KKw]
                FGval[0]+=gAllUnits[KKpreviousLayer+JJ,1]*vw[KKw]
                FGval[1]+=gAllUnits[KKpreviousLayer+JJ,2]*vw[KKw]
                FGval[2]+=gAllUnits[KKpreviousLayer+JJ,3]*vw[KKw]
                KKw=KKw+1  
            vAllUnits[KKcurrentLayer+KKunit]=ssum**vDegree[KKlayer]
            gsum=vDegree[KKlayer]*ssum**gDegree[KKlayer]
            gAllUnits[KKcurrentLayer+KKunit,1:]=gsum*FGval[:]
        KKpreviousLayer=KKcurrentLayer
        KKcurrentLayer=KKcurrentLayer+vUnits[KKlayer]   
    ssum,FGval[:]=0.0,0.0
    for JJ in range(1,vUnits[nnLayers]+1): ### Calculate output
        ssum=ssum+vAllUnits[KKpreviousLayer+JJ]*vw[KKw]
        FGval[0]+=gAllUnits[KKpreviousLayer+JJ,1]*vw[KKw]
        FGval[1]+=gAllUnits[KKpreviousLayer+JJ,2]*vw[KKw]
        FGval[2]+=gAllUnits[KKpreviousLayer+JJ,3]*vw[KKw]
        KKw=KKw+1
    homDegReciproc=dictParam['homDegReciproc']    
    Fval=ssum**homDegReciproc ## output value
    ssum=homDegReciproc*Fval/ssum
    ##FGval[:]=ssum*FGval[:]	
    return  Fval, ssum*FGval[:]     
    
def plotPointsPolyN(zvcf,vphi):
    degree=int(zvcf[5])
    vcos=znp.cos(vphi);vsin=znp.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    #vsx,vsy,vsxy=np.zeros(Nphi),np.zeros(Nphi),np.zeros(Nphi)
    nPhi=vphi.shape[0]
    vR = znp.zeros(nPhi)
    for k in range(nPhi):
        vsx,vsy,vsxy=vcos2[k],vsin2[k],vsc[k]
        yf,[vDX,vDY,vDXY]=fGYF(vsx,vsy,vsxy,zvcf[5:])            
        vR[k]=(vDXY*vsc[k]-vDX*vsin2[k]-vDY*vcos2[k])/(vDX+vDY)
    return degree,vR    


def plotPointsNNys(ddParam,vphi):
    vcos=znp.cos(vphi);vsin=znp.sin(vphi)
    vcos2,vsin2,vsc=vcos*vcos,vsin*vsin,vcos*vsin
    nPhi=vphi.shape[0]
    vR = znp.zeros(nPhi)
    for k in range(nPhi):
        vsx,vsy,vsxy=vcos2[k],vsin2[k],vsc[k]
        yf,[vDX,vDY,vDXY]=nnfGYF([vsx,vsy,vsxy],ddParam)            
        vR[k]=(vDXY*vsc[k]-vDX*vsin2[k]-vDY*vcos2[k])/(vDX+vDY)
    return vR

def plotOneTest(baseName='',vAngle=[],vfStress=[],vfNode=[],vfJob=[],zvcf=[],casePoly=False,ddParam={},caseNNys=False):
    sep='|'
    epsTol=1.0e-6
    matE=aMAT['E']
    matN=aMAT['nuPoisson']
    vMaxStrains=strainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngle,hratio=HRATIO)
    if(casePoly):
        if('symb' in baseName):
            labelPoly='Poly{}-Hill fitted r-value'.format(int(zvcf[5]))
        else:
            labelPoly='Poly{} fitted r-value'.format(int(zvcf[5]))
        ##vMaxStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngle,vcf=zvcf[5:])
    elif(caseNNys):
        labelPoly='nnYS fitted r-value'    
    else:    
        labelPoly='Hill fitted r-value'
        ##vMaxStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngle)
        GH,FH,H,N2=paramPlotHill(aMAT['rVals'])
    vrv={}
    fg1,fg2,fg3=plt.figure(),plt.figure(),plt.figure()
    for idxAngle in range(len(vAngle)):
        ANGLE_FROM_RD=int(vAngle[idxAngle])
        strAngle=str(int(ANGLE_FROM_RD))
        fileStress,fileNode,fileJob=vfStress[idxAngle],vfNode[idxAngle],vfJob[idxAngle]
        fileRval=baseName+'_Rvals.txt'
        fileJob2=baseName+'_JOB.txt'
        ss,ss2,ee11,ee22,eep,skipL=[],[],[],[],[],0
        try:
            ff=open(fileStress,'r')
            sAngle=int(float(ff.readline().strip().split(':')[-1].strip()))
            ff.readline()
            for line in ff:
                vline=line.strip().split(sep)
                eps=float(vline[4])
                if(float(eps<epsTol)):skipL+=1;continue
                ss.append(float(vline[0]))
                ss2.append(float(vline[1]))
                ee11.append(eps)
                ee22.append(float(vline[5]))
                eep.append(float(vline[-2]))        
            ff.close()
            ee11,ee22,ss,ss2,eep=znp.array(ee11),znp.array(ee22),znp.array(ss),znp.array(ss2),znp.array(eep)
            ze11=(ss-matN*ss2)/matE ##elastic strains 
            ze22=(ss2-matN*ss)/matE
            ee11-=ze11
            ee22-=ze22
        except IOError as err:
            printMsg(err) 
            printMsg('Skipping test angle '+strAngle);continue
        #
        if(not shell_S4R):##Transform Green-Lagrange strains to log-strains 
            ee11=0.5*znp.log(1.0+2.0*ee11)
            ee22=0.5*znp.log(1.0+2.0*ee22)
        #    
        nee11,nee22,LL,WW=[],[],0.0,0.0
        d11,d22=[],[]
        nAe11,nBe11,nCe22,nDe22=[],[],[],[]
        nAx,nBx,nCy,nDy=zro,zro,zro,zro
        try:
            ff=open(fileNode,'r')
            nAngle=int(float(ff.readline().strip().split(':')[-1].strip()))
            ff.readline()
            #line=ff.readline()
            #vline=line.strip().split(sep)
            #LL,WW=float(vline[0]),float(vline[1])
            nAx=float(ff.readline().strip().split(sep)[0])
            ff.readline()
            nBx=float(ff.readline().strip().split(sep)[0])
            ff.readline()
            nCy=float(ff.readline().strip().split(sep)[1])
            ff.readline()
            nDy=float(ff.readline().strip().split(sep)[1])
            LL,WW=nAx-nBx,nCy-nDy
            for kk in range(skipL+1):
                ff.readline()
            kk=0    
            for line in ff:
                vline=line.strip().split(sep)
                dAX,dBX,dCY,dDY=float(vline[0]),float(vline[3]),float(vline[7]),float(vline[10])
                nee11.append(znp.log((LL+dAX-dBX)/LL)-ze11[kk])
                nee22.append(znp.log((WW+dCY-dDY)/WW)-ze22[kk])
                kk+=1
                #d11.append(DX)
                #d22.append(DY)
            ff.close()    
        except IOError as err:
            printMsg(err)
            printMsg('Skipping test angle '+strAngle);continue            
        #   
        vrv={}
        try:
            ff=open(fileRval,'r')
            for line in ff:
                vline=line.strip().split(sep)
                vrv[str(int(float(vline[0])))]=[float(vline[1]),float(vline[2])]
            ff.close()
        except IOError:
            pass
        #
        if(casePoly):
            tt=vAngle[idxAngle]*rad
            ct,st=mcos(tt),msin(tt)
            ft=fYF(ct*ct,st*st,ct*st,zvcf[5:])
        elif(caseNNys):
            tt=vAngle[idxAngle]*rad
            ct,st=mcos(tt),msin(tt)
            ft=nnfYF([ct*ct,st*st,ct*st],ddParam)        
        else:        
            ft=fHill(GH,FH,H,N2,sAngle*znp.pi/180.0)  
        #
        djb={}
        try:
            ff=open(fileJob2,'r')
            ff.readline()
            for line in ff:
                if('-' not in line):
                    vline=line.strip().split(sep)
                    djb[str(int(float(vline[0])))]=int(float(vline[1]))
            ff.close()
        except IOError:
            pass
        try:
            ff=open(fileJob,'r')
            ff.readline()
            vline=ff.readline().strip().split(sep)
            #if(ANGLE_FROM_RD != int(float(vline[0]))):
            #    print('Test angle in job file and this script are not equal')
            #    exit()
            djb[str(ANGLE_FROM_RD)]=int(float(vline[1]))    
            ff.close()
        except IOError:
            pass 
        ff=open(fileJob2,'w')
        ff.write('ANGLE_FROM_RD | WallClock (secs)\n')
        timeTotal=0
        for key in djb:
            timeTotal+=djb[key]
            ff.write('{:} | {:}\n'.format(str(key),str(djb[key])))
        ff.write('---TOTAL: '+str(timeTotal))    
        ff.close()
        #
        ###find the last stable iuncrement (if the case) 
        nRval=0
        while(nRval<len(nee11)-1):
            if(ss[nRval]>=ss[nRval+1]):break
            nRval+=1
        ##print('nRval: ',nRval)  
        nRval-=3 
        #
        vrv[str(sAngle)]=[-nee22[nRval]/(nee11[nRval]+nee22[nRval]),-ee22[nRval]/(ee11[nRval]+ee22[nRval])]
        ff=open(fileRval,'w')
        for key in vrv:
            ff.write(key+"|{:.3f}|{:.3f}\n".format(*vrv[key]))
        ff.close()
        ##print('angle from RD: ',sAngle)
        ##print('rValue:  ',rVnode)
        ##lStrain=np.log(LLdef/LL)
        ##print('rValue2: ', -wStrain/(lStrain+wStrain))
        #
        vRval,vnRval=[],[]
        for k in range(1,nRval+1):
            vRval.append(-ee22[k]/(ee11[k]+ee22[k]))
            vnRval.append(-nee22[k]/(nee11[k]+nee22[k]))
        #            
        vse=znp.array(hardFunc(aMAT['hardLaw'],aMAT['hardParam'],vMaxStrains[idxAngle]+0.025))
        #    
        ###fg1=plt.figure()
        ax=fg1.add_subplot()
        ax.plot([],[],linestyle='',label=r'$Test\,\,\, angle: \theta = $'+strAngle)
        ##ax.plot(vse[:,1],vse[:,0]/ft,color='k',linewidth=1,linestyle='dashed',label='$\sigma = H(\epsilon) = K(\epsilon_0+\epsilon)^n$')
        label=r'$\sigma = k_{'+strAngle+r'}H\left(\overline{\epsilon}\right)$'
        ax.plot(vse[:,1],vse[:,0]/ft,color='k',linewidth=2,linestyle='dashed',label=label)
        ###ax.plot(nee11,ss,linestyle='',marker='+',color='b',markersize=6,label='FE (strain from node)')
        #if(UMAT):
        #    ax.plot(eep,ss,linestyle='',marker='+',color='r',markersize=6,label='FE (strain from element)')
        #else:
        #    ax.plot(ee11,ss,linestyle='',marker='+',color='r',markersize=6,label='FE (strain from element)')    
        ax.plot(eep,ss,linestyle='',marker='+',color='r',markersize=6,label='FE (strain from element)')
        ax.legend(loc='lower right',fontsize=15)
        #ax.set_xlabel('$\epsilon_{11}$')
        #ax.set_ylabel('$\sigma_{11}$')   
        ax.text(eep[0],znp.max(vse[:,0])-20,r'$\sigma_{11}$',fontsize=16)
        ax.text(0.3*eep[-1],ss[0],r'$\overline{\epsilon}$',fontsize=16)
        ax.tick_params(axis='both', which='major', labelsize=12)
        #
        ax3=fg3.add_subplot()
        ax3.plot([],[],linestyle='',label=r'$Test\,\,\, angle: \theta = $'+strAngle)
        #ax3.plot(nee11[1:nRval+1],vnRval,linestyle='',marker='+',markersize=6,color='b',label='FE (strain from node)')
        #ax3.plot(ee11[1:nRval+1],vRval,linestyle='',marker='+',markersize=6,color='r',label='FE (strain from element)')
        ax3.plot(eep[1:nRval+1],vnRval,linestyle='',marker='+',markersize=6,color='b',label='FE (strain from nodes)')
        ax3.plot(eep[1:nRval+1],vRval,linestyle='',marker='+',markersize=6,color='r',label='FE (strain from element)')
        if(casePoly):
            degree,targetRval=plotPointsPolyN(zvcf,znp.array([sAngle*znp.pi/180,]))
            targetRval=targetRval[0]
        elif(caseNNys):
            targetRval=plotPointsNNys(ddParam,znp.array([sAngle*znp.pi/180,]))[0]        
        else:    
            targetRval=rvalHill2(sAngle*znp.pi/180.0,GH,FH,H,N2)
        ax3.plot([0,ee11[-1]+0.025],[targetRval,targetRval],linewidth=2,linestyle='dashed',color='k',label=labelPoly)
        #ax3.text(0.3*ee11[-1],max([targetRval]+vRval+vnRval),'$r_{'+str(sAngle)+'}$',fontsize=16)
        #ax3.text(0.3*ee11[-1],min([targetRval]+vRval+vnRval),r'$\overline{\epsilon}$',fontsize=14)
        if(0):
            ax3.text(0.12*ee11[-1],0.58,'$r_{'+str(sAngle)+'}$',fontsize=17)
            ax3.text(0.31*ee11[-1],0.115,r'$\overline{\epsilon}$',fontsize=16)
            ax3.tick_params(axis='both', which='major', labelsize=13)
            ax3.set_ylim(0.1,0.65)
            ax3.legend(loc='lower right',fontsize=15)            
        else:
            ax3.text(0.3*ee11[-1],max([targetRval]+vRval+vnRval),r'$r_{'+str(sAngle)+'}$',fontsize=16)
            ax3.text(0.3*ee11[-1],min([targetRval]+vRval+vnRval),r'$\overline{\epsilon}$',fontsize=14)
            ax3.legend()            
        #
        if(True):
            fg1.savefig(baseName+'_Hcurve_{}.png'.format(str(sAngle)),dpi=300,bbox_inches='tight')
            #fg2.savefig(ffz+'_RVals.png',dpi=300,bbox_inches='tight')
            fg3.savefig(baseName+'_RVals_{}.png'.format(str(sAngle)),dpi=300,bbox_inches='tight')
            fg1.clf()
            #fg2.clf()
            fg3.clf()
        #if(False):
        #    plt.show()
    ax2=fg2.add_subplot()
    vtheta=znp.linspace(0,znp.pi/2,nPlotPoints)
    if(casePoly):        
        degree,rtheta=plotPointsPolyN(zvcf,vtheta)
    elif(caseNNys):
        rtheta=plotPointsNNys(ddParam,vtheta)    
    else:    
        rtheta=rvalHill2(vtheta,GH,FH,H,N2)
    #ax2.plot(aMAT['rExp'][0],aMAT['rExp'][1],linestyle='',marker='o',markersize=10,color='g',markerfacecolor=(1,1,1,0),markeredgecolor='g',label='Exp. data')    
    ax2.plot((180/znp.pi)*vtheta,rtheta,linewidth=2,linestyle='dashed',color='k',label=labelPoly)
    for key in vrv:
        ax2.plot(float(key),vrv[key][0],linestyle='',marker='+',markersize=10,color='b',markeredgewidth=2)
        ax2.plot(float(key),vrv[key][1],linestyle='',marker='x',markersize=10,color='r',markeredgewidth=2)
    if(vrv):    
        ax2.plot(float(key),vrv[key][0],linestyle='',marker='+',markersize=10,color='b',markeredgewidth=2,label='FE (strain from nodes)')
        ax2.plot(float(key),vrv[key][1],linestyle='',marker='x',markersize=10,color='r',markeredgewidth=2,label='FE (strain from element)')
    ax2.set_xticks([0,15,30,45,60,75,90])
    ax2.tick_params(axis='both', which='major', labelsize=13)
    ax2.legend(fontsize=15,borderaxespad=0,borderpad=0.2)
    ax2.set_xlim(-2,92)
    ymin,ymax=ax2.get_ylim()
    ax2.text(87,ymin+0.01*(ymax-ymin),r'$\theta$',fontsize=16)
    ax2.text(1,0.5*(ymax+ymin),r'$r_{\theta}$',fontsize=16)
    fg2.savefig(baseName+'_RVals.png',dpi=300,bbox_inches='tight');plt.close(fg1);plt.close(fg2);plt.close(fg3)  
    printMsg('----------Reports and plots for {}: DONE -------------------------------------------'.format(baseName))



def read_NN_InpFile(wFile):
    NTENS=int(3)
    vNNparam=[0.0,] ### Added a first element to start counting from 1
    try:         
        fileNNparam=open(wFile,'r')
        for line in fileNNparam:
            vNNparam.append(float(line.strip()))
        fileNNparam.close()
        nnLayers=int(vNNparam[6])
        vUnits,vDegree=znp.zeros(nnLayers+1,dtype=znp.int32),znp.zeros(nnLayers+1) ### Added a first element to start counting from 1
        vUnits[1]= NTENS  ###!!! input dimension
        nAllUnits=NTENS
        nnHom=1.0
        JJ=7
        for KK in range(2,nnLayers+1):
            nUnits=int(vNNparam[JJ])
            vUnits[KK]=nUnits
            nAllUnits+=nUnits
            vDegree[KK]=vNNparam[JJ+1]
            nnHom=nnHom*vDegree[KK]
            JJ+=2
        firstWeight=JJ
        homDegReciproc=1.0/nnHom 
        return {'nnLayers':int(nnLayers),'nAllUnits':int(nAllUnits),'homDegReciproc':homDegReciproc,
                'firstWeight':int(firstWeight),'vUnits':vUnits,'vDegree':vDegree,'vNNparam':vNNparam}
        fileNNparam.close()        
    except IOError:
        printMsg('Input file for abqUMAT_nnYS not found: '+wFile)
        return False        

#
def readPlotData(bParam,inpData,ZHRATIO=6.0,ZsubDir=''):
    global aMAT,hhh,vvTests,noGUI,HRATIO,figDirData,figDirPlot
    if(not uaxSim):printMsg('readPlotData: Only for uniaxial sims\nPost-proc aborted');sysexit()
    noGUI,HRATIO=False,ZHRATIO
    ZsubDir=bParam['subDir']
    if(ZsubDir):
        if(osn=='nt'):figDirData+=ZsubDir+'\\'
        else:figDirData+=ZsubDir+'/'    
    if(not (ospath.exists(figDirData) and ospath.isdir(figDirData))):
        printMsg("The local data folder was not found")
        printMsg("Make sure a folder {} exists at your script location".format(figDirData))
        printMsg("Calculations aborted");sysexit()    
    if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
        printMsg("The local folder for saving reports and plots was not found")
        printMsg("Make sure a folder \'RESULTS\' exists at your script location")
        printMsg("Calculations aborted");sysexit()
    if(ZsubDir):
        if(osn=='nt'):figDirPlot+=ZsubDir+'\\'
        else:figDirPlot+=ZsubDir+'/'
        if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
            try:
                osmaked(figDirPlot)
            except OSError as err:
                printMsg('readPlotData:');printMsg(err);sysexit()
    if(bParam['hLaw'] in ['Swift','Voce']):
        aMAT['hardLaw']=bParam['hLaw']
    else:
        printMsg("Only 'Swift' and 'Voce' hardening laws can be used (at the moment)")
        printMsg('Calculations aborted');sysexit()        
    aMAT['hThick']=bParam['hThick']
    aMAT['E'],aMAT['nuPoisson'],aMAT['muG']=bParam['eParam']
    aMAT['hardParam'][0],aMAT['hardParam'][1],aMAT['hardParam'][2]=bParam['hParam']
    aMAT['rVals'][0],aMAT['rVals'][1],aMAT['rVals'][2]=bParam['rParam']
    hhh=bParam['hThick']
    for rrr in bParam['rExp']:
        aMAT['rExp'].append(rrr)
    ##printMsg(str(transvK));sysexit()
    for ddd in inpData:
        dfile=ddd['inpFile']
        if(not dfile):printMsg("field 'inpFile' cannot be empty\nCalculations aborted");sysexit()
        if(ddd['PolyN']):
            zdfile=figDirData+dfile
            if(not(ospath.exists(zdfile) and ospath.isfile(zdfile))):
                printMsg('File '+dfile+' not found\nCalculations aborted')
                sysexit()
        tdd={'PolyN':ddd['PolyN'],'dfile':dfile,'angles':ddd['angles']}
        vvTests.append(tdd)        
    printMsg('------------------------------------------------------------Plot input data: OK')
    ####sysexit()
#    
def plotTests():
    if(not uaxSim):printMsg('plotTests: Only for uniaxial sims\nPost-proc aborted');sysexit()
    for ttest in vvTests:
        vAngles=ttest['angles']
        ff=figDirPlot+ttest['dfile'].strip('.txt')
        vFileJob,vFileNode,vFileStress=[],[],[]
        for idxAngle in range(len(vAngles)):
            ANGLE_FROM_RD=vAngles[idxAngle]        
            strAngle=str(int(ANGLE_FROM_RD))
            vFileStress.append(ff+'_ES_DATA_'+strAngle+'.txt')
            vFileNode.append(ff+'_Node_'+strAngle+'.txt')
            vFileJob.append(ff+'_JOB_'+strAngle+'.txt')
        ##vStrains=strainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles,hratio=HRATIO)
        zvcf,casePoly,ddParam,caseNNys=[],ttest['PolyN'],{},ttest['caseNNys']
        if(caseNNys):
            casePoly=False
            ddParam=read_NN_InpFile(figDirData+ttest['dfile'])
            if(not ddParam):continue  
        if(casePoly):
            try:
                zff=open(figDirData+ttest['dfile'],'r')
                for line in zff:
                    zvcf.append(float(line.strip()))
                zff.close()    
            except IOError:
                printMsg('Input file for Poly_UMAT not found: '+ttest['dfile'])
                continue          
        plotOneTest(baseName=ff,vAngle=vAngles,vfStress=vFileStress,vfNode=vFileNode,vfJob=vFileJob,zvcf=zvcf,casePoly=casePoly,ddParam=ddParam,caseNNys=caseNNys)
#
###set these as properties of bParam and transfer them via readData to aMAT or global vars 
thisDDname='abqPolyNcpdrwShell'
modelDDname=thisPATH+thisDDname+".cae" 
blankName,dieName,punchName,holderName='Blank','Die','Punch','Holder'
blankInstName,dieInstName,punchInstName,holderInstName='Blank-1','Die-1','Punch-1','Holder-1'
def genDDabqModel(jbName,pTravel,zdtInitial,zdtMax,zmaxNInc,zpSoft,zdeltaPunch=0):
    session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
    ##Cup drawing global params (must be set via var 'bParam' in the 'main' driver)
    blankDiam=aMAT['blankDiam']
    dieOpeningDiam=aMAT['dieOpeningDiam']
    dieShoulderRad=aMAT['dieShoulderRad']
    punchDiam=aMAT['punchDiam']
    punchNoseRad=aMAT['punchNoseRad']
    holderClearance=aMAT['holderClearance']
    frictionCoeff=aMAT['frictionCoeff']
    blankRad=blankDiam/2.0
    flangeRad=blankRad+20.0
    #dieOpeningRad=dieOpeningDiam/2.0-aMAT['hThick']
    punchRad=punchDiam/2.0
    holderRad=punchRad+1
    holderShoulderRad=dieShoulderRad/2.0
    holderClearance-=aMAT['hThick']
    #epsContact=0.0005*holderClearance
    #zBlank=epsContact+0.5*holderClearance
    epsContact=0.001
    dieOpeningRad=punchRad+9.0*epsContact
    zBlank=epsContact
    dzPunch=3.0*epsContact
    #
    ctctSoft,ctctSoft2,ctctSoft3,ctctHard='contactSoft','ctctSoft2','ctctSoft3','contactHard'
    mmd=mdb.Model(name=mdbModelName)
    sk=mmd.ConstrainedSketch(name='pProfile', sheetSize=1.5*blankDiam)
    sk.Line(point1=(0.0,0.0),point2=(0.0,blankRad))
    sk.ArcByCenterEnds(center=(0.0, 0.0), point1=(0.0, blankRad), point2=(blankRad,0.0), 
    direction=CLOCKWISE)    
    sk.Line(point1=(blankRad, 0.0),point2=(0.0,0.0))
    pPart=mmd.Part(name=blankName, dimensionality=THREE_D, type=DEFORMABLE_BODY)
    pPart.BaseShell(sketch=sk)
    ##pPart.BaseSolidExtrude(sketch=sk, depth=aMAT['hThick'])
    pPart.Set(faces=pPart.faces,name='wholeBlankSection')
    del mmd.sketches['pProfile']
    sk=mmd.ConstrainedSketch(name='pProfile',sheetSize=2.0*flangeRad)
    sk.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    sk.Line(point1=(dieOpeningRad,-max(dieShoulderRad,punchNoseRad)-dieShoulderRad),
    point2=(dieOpeningRad,-dieShoulderRad))
    sk.ArcByCenterEnds(center=(dieOpeningRad+dieShoulderRad,-dieShoulderRad),
    point1=(dieOpeningRad,-dieShoulderRad),point2=(dieOpeningRad+dieShoulderRad,0.0),
    direction=CLOCKWISE)
    sk.Line(point1=(dieOpeningRad+dieShoulderRad,0.0),point2=(flangeRad,0.0))
    pPart=mmd.Part(name=dieName, dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
    pPart.AnalyticRigidSurfRevolve(sketch=sk)
    del mmd.sketches['pProfile']
    pPart.ReferencePoint(point=(0.0, 0.0, 1.0))
    sk=mmd.ConstrainedSketch(name='pProfile',sheetSize=2.0*flangeRad)
    sk.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
    sk.Line(point1=(punchRad,blankRad+punchNoseRad),point2=(punchRad,punchNoseRad))
    sk.ArcByCenterEnds(center=(punchRad-punchNoseRad,punchNoseRad),
    point1=(punchRad,punchNoseRad),point2=(punchRad-punchNoseRad,0.0),
    direction=CLOCKWISE)
    sk.Line(point1=(punchRad-punchNoseRad,0.0),point2=(0.0,0.0))
    pPart=mmd.Part(name=punchName, dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
    pPart.AnalyticRigidSurfRevolve(sketch=sk)
    del mmd.sketches['pProfile']
    pPart.ReferencePoint(point=(0.0, 0.0, 0.0))
    ##sk=mmd.ConstrainedSketch(name='pProfile',sheetSize=2.0*flangeRad,
    ##transform=(1.0,0.0,0.0, 0.0,0.0,1.0, 0.0,-1.0,0.0, 0.0,0.0,0.0))
    ##sk.assignCenterline(sk.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0)))
    sk=mmd.ConstrainedSketch(name='pProfile',sheetSize=2.0*flangeRad)
    sk.ConstructionLine(point2=(0.0, -100.0), point1=(0.0, 100.0))
    sk.Line(point1=(flangeRad,0.0),point2=(holderRad+holderShoulderRad,0.0))
    sk.ArcByCenterEnds(center=(holderRad+holderShoulderRad,holderShoulderRad),
    point1=(holderRad+holderShoulderRad,0.0),point2=(holderRad,holderShoulderRad),
    direction=CLOCKWISE)
    sk.Line(point1=(holderRad,holderShoulderRad),point2=(holderRad,holderShoulderRad+5.0))
    ##sk.ArcByCenterEnds(center=(holderRad+holderShoulderRad,holderShoulderRad),
    ##point2=(holderRad+hoCOUNTERCLOCKWISE)lderShoulderRad,0.0),point1=(holderRad,holderShoulderRad),
    pPart=mmd.Part(name=holderName, dimensionality=THREE_D, type=ANALYTIC_RIGID_SURFACE)
    pPart.AnalyticRigidSurfRevolve(sketch=sk)
    del mmd.sketches['pProfile']
    pPart.ReferencePoint(point=(0.0, 0.0, 1.0))
    ara=mmd.rootAssembly
    ara.DatumCsysByDefault(CARTESIAN)
    ara.Instance(name=blankInstName, part=mmd.parts[blankName],dependent=OFF)
    ara.Instance(name=dieInstName,part=mmd.parts[dieName],dependent=OFF)
    ara.Instance(name=punchInstName,part=mmd.parts[punchName],dependent=OFF)
    ara.Instance(name=holderInstName,part=mmd.parts[holderName],dependent=OFF)
    ara.rotate(instanceList=(dieInstName,punchInstName,holderInstName), 
    axisPoint=(0.0, 0.0, 0.0), axisDirection=(10.0, 0.0, 0.0), angle=90.0)
    ##ara.rotate(instanceList=(holderInstName,), 
    ##axisPoint=(0.0, 0.0, 0.0), axisDirection=(10.0, 0.0, 0.0), angle=-90.0)
    #ara.translate(instanceList=(holderInstName,),vector=(0.0,0.0,holderClearance+2*epsContact))
    ara.translate(instanceList=(holderInstName,),vector=(0.0,0.0,3.0*epsContact))
    ara.translate(instanceList=(blankInstName,),vector=(0.0,0.0,zBlank))
    #ara.translate(instanceList=(punchInstName,),vector=(0.0,0.0,zBlank+dzPunch))
    ara.translate(instanceList=(punchInstName,),vector=(0.0,0.0,dzPunch))
    pPart=ara.instances[blankInstName].part
    ##wholeBlank=pPart.cells.getByBoundingBox(-gTol,-gTol,-gTol,blankRad+gTol,blankRad+gTol,aMAT['hThick']+gTol)
    ##wholeBlank=pPart.cells.getByBoundingBox(-gTol,-gTol,-gTol,blankRad+gTol,blankRad+gTol,aMAT['hThick']+gTol)
    matCSYS=pPart.DatumCsysByThreePoints(coordSysType=CARTESIAN, name='csys-1-Material', origin=(zro,zro,zro), 
            point1=(blankRad,zro,zro), point2=(zro,blankRad,zro))    
    pPart.MaterialOrientation(additionalRotationField='', additionalRotationType=ROTATION_ANGLE, 
                angle=zro, axis=AXIS_3, fieldName='', localCsys=pPart.datums[matCSYS.id], orientationType=SYSTEM, 
                region=pPart.faces)
    ####
    x3Punch2=punchNoseRad
    if(pTravel>0.0):
        x3Punch=pTravel+dieShoulderRad
    else:printMsg('punchTravel must be >0');sysexit()
    dtInitial,dtMax,maxNInc=zdtInitial,zdtMax,zmaxNInc
    if(zdtInitial<0.0):printMsg('dtInitial must be >0');sysexit()
    if(zdtMax<0.0):printMsg('dtMax must be >0');sysexit()
    if(zmaxNInc<1):printMsg('maxNInc must be >1');sysexit() 
    mmd.StaticStep(description='Make H_contact', initialInc=1.0e-2, 
        matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-16,maxInc=0.025, maxNumInc=1000, 
        name='Step-1', nlgeom=ON, previous='Initial')   
    mmd.StaticStep(description='Make P_contact', initialInc=1.0e-6, 
        matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-16,maxInc=0.025, maxNumInc=1000, 
        name='Step-2', nlgeom=ON, previous='Step-1')           
    mmd.StaticStep(description='Move_1 punch', initialInc=dtInitial, 
        matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-16,maxInc=dtMax, maxNumInc=maxNInc, 
        name='Step-3', nlgeom=ON, previous='Step-2')
    mmd.StaticStep(description='Remove tools: Die', initialInc=0.001, 
        matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-16,maxInc=0.005, maxNumInc=maxNInc, 
        name='Step-4', nlgeom=ON, previous='Step-3')    
    #mmd.StaticStep(description='Remove tools: Die', initialInc=0.001, 
    #    matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-16,maxInc=0.001, maxNumInc=maxNInc, 
    #    name='Step-4', nlgeom=ON, previous='Step-3')        
    #mmd.StaticStep(description='Remove tools: Punch', initialInc=0.01, 
    #    matrixSolver=DIRECT, matrixStorage=SYMMETRIC, minInc=1.0e-9,maxInc=0.05, maxNumInc=maxNInc, 
    #    name='Step-5', nlgeom=ON, previous='Step-4') 
    #
    mmd.FieldOutputRequest(name='F-Output-3',createStepName='Step-3')
    mmd.FieldOutputRequest(name='F-Output-4',createStepName='Step-4')
    #mmd.FieldOutputRequest(name='F-Output-4',createStepName='Step-5')
    ##mmd.steps['Step-3'].setValues(matrixSolver=DIRECT,matrixStorage=UNSYMMETRIC)
    if(0):
        mmd.steps['Step-3'].control.setValues(allowPropagation=OFF,resetDefaultValues=OFF,
        discontinuous=ON, 
        displacementField=(0.005, 0.05, 0.0, 0.0, 0.02, 1e-05, 0.001, 1e-08, 1.0, 1e-05, 1e-08), 
        rotationField=(0.005, 0.05, 0.0, 0.0, 0.02, 1e-05, 0.001, 1e-08, 1.0, 1e-05),
        electricalPotentialField=DEFAULT,hydrostaticFluidPressureField=DEFAULT, 
        lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1))
        mmd.steps['Step-4'].control.setValues(allowPropagation=OFF,resetDefaultValues=OFF,
        discontinuous=ON, 
        displacementField=(0.005, 0.05, 0.0, 0.0, 0.02, 1e-05, 0.001, 1e-08, 1.0, 1e-05, 1e-08),
        rotationField=(0.005, 0.05, 0.0, 0.0, 0.02, 1e-05, 0.001, 1e-08, 1.0, 1e-05),    
        electricalPotentialField=DEFAULT,hydrostaticFluidPressureField=DEFAULT, 
        lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1))
    for kk in mmd.historyOutputRequests.keys():
        del mmd.historyOutputRequests[kk]        
    ###
    mmd.steps['Step-3'].Restart(frequency=0,numberIntervals=1,overlay=ON,timeMarks=ON)
    #mmd.steps['Step-4'].Restart(frequency=0,numberIntervals=1,overlay=ON,timeMarks=ON)
    pSoftH,pSoftD,pSoftP=zpSoft[:]
    yldZero=ffHardFunc(aMAT['hardLaw'],aMAT['hardParam'],zro)
    mmd.ContactProperty(ctctSoft)
    mmd.interactionProperties[ctctSoft].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((frictionCoeff, ), ), 
    shearStressLimit=None, maximumElasticSlip=FRACTION,fraction=0.005, elasticSlipStiffness=None)
    mmd.interactionProperties[ctctSoft].NormalBehavior(pressureOverclosure=EXPONENTIAL, 
    table=((pSoftH*yldZero,0.0), (0.0, zBlank)),
    maxStiffness=None, constraintEnforcementMethod=DEFAULT)
    mmd.ContactProperty(ctctSoft2)
    mmd.interactionProperties[ctctSoft2].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((frictionCoeff, ), ), 
    shearStressLimit=None, maximumElasticSlip=FRACTION,fraction=0.005, elasticSlipStiffness=None)
    mmd.interactionProperties[ctctSoft2].NormalBehavior(pressureOverclosure=EXPONENTIAL, 
    table=((pSoftD*yldZero, 0.0), (0.0, zBlank)),
    maxStiffness=None, constraintEnforcementMethod=DEFAULT)
    mmd.ContactProperty(ctctSoft3)
    mmd.interactionProperties[ctctSoft3].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((frictionCoeff, ), ), 
    shearStressLimit=None, maximumElasticSlip=FRACTION,fraction=0.005, elasticSlipStiffness=None)
    mmd.interactionProperties[ctctSoft3].NormalBehavior(pressureOverclosure=EXPONENTIAL, 
    table=((pSoftP*yldZero, 0.0), (0.0, zBlank)),
    maxStiffness=None, constraintEnforcementMethod=DEFAULT)
    mmd.ContactProperty(ctctHard)
    mmd.interactionProperties[ctctHard].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, table=((
        frictionCoeff, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    mmd.interactionProperties[ctctHard].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)   
    regionDie=ara.Surface(side1Faces=ara.instances[dieInstName].faces,name='master_Surf_Die')
    regionHolder=ara.Surface(side1Faces=ara.instances[holderInstName].faces,name='master_Surf_Holder')
    regionPunch=ara.Surface(side1Faces=ara.instances[punchInstName].faces,name='master_Surf_Punch')
    regionBlankTop=ara.Surface(side1Faces=ara.instances[blankInstName].faces,name='slave_Surf_BlankTop')
    regionBlankBottom=ara.Surface(side2Faces=ara.instances[blankInstName].faces,name='slave_Surf_BlankBottom')
    #mmd.SurfaceToSurfaceContactStd(name='Die-Blank', 
    #createStepName='Initial', master=regionDie, slave=regionBlankBottom, sliding=FINITE, 
    #thickness=ON, interactionProperty=ctctSoft, adjustMethod=NONE, enforcement=SURFACE_TO_SURFACE,
    #initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF)
    try:
        mmd.SurfaceToSurfaceContactStd(name='Die-Blank', 
        createStepName='Initial', master=regionDie, slave=regionBlankBottom, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft2, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)
        mmd.SurfaceToSurfaceContactStd(name='Holder-Blank', 
        createStepName='Initial', master=regionHolder, slave=regionBlankTop, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)
        mmd.SurfaceToSurfaceContactStd(name='Punch-Blank', 
        createStepName='Step-2', master=regionPunch, slave=regionBlankTop, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft3, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)
    except:
        mmd.SurfaceToSurfaceContactStd(name='Die-Blank', 
        createStepName='Initial', main=regionDie, secondary=regionBlankBottom, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft2, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)
        mmd.SurfaceToSurfaceContactStd(name='Holder-Blank', 
        createStepName='Initial', main=regionHolder, secondary=regionBlankTop, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)
        mmd.SurfaceToSurfaceContactStd(name='Punch-Blank', 
        createStepName='Step-2', main=regionPunch, secondary=regionBlankTop, sliding=FINITE, 
        thickness=OFF, interactionProperty=ctctSoft3, adjustMethod=NONE, enforcement=NODE_TO_SURFACE,
        initialClearance=OMIT, datumAxis=None, clearanceRegion=None, tied=OFF,smooth=0.1, bondingSet=None)    
    if(1):
        mmd.StdContactControl(name='ContCtrl-1', stabilizeChoice=AUTOMATIC, 
        stiffnessScaleFactor=10000.0,penetrationTolChoice=ABSOLUTE, absolutePenetrationTolerance=0.0001) 
        mmd.interactions['Die-Blank'].setValuesInStep(stepName='Step-1',contactControls='ContCtrl-1')
        mmd.interactions['Holder-Blank'].setValuesInStep(stepName='Step-1',contactControls='ContCtrl-1')
        mmd.interactions['Punch-Blank'].setValuesInStep(stepName='Step-2',contactControls='ContCtrl-1')
    ###mmd.interactions['Punch-Blank'].deactivate('Step-5')    
    refP=ara.instances[punchInstName].referencePoints
    ara.Set(referencePoints=(refP[2],),name='punchRP')
    ara.Set(referencePoints=(ara.instances[dieInstName].referencePoints[2],),name='dieRP')
    ara.Set(referencePoints=(ara.instances[holderInstName].referencePoints[2],),name='holderRP')
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_Die',
        region=ara.sets['dieRP'], 
        u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)    
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_Holder',
        region=ara.sets['holderRP'], 
        u1=SET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET)
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Fix_Punch', 
        region=ara.sets['punchRP'], 
        u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET)
    #mmd.boundaryConditions['Fix_Holder'].setValuesInStep(stepName='Step-1',u3=-2.0*epsContact)
    mmd.boundaryConditions['Fix_Holder'].setValuesInStep(stepName='Step-1',u3=-epsContact)    
    #
    mmd.boundaryConditions['Fix_Punch'].setValuesInStep(stepName='Step-2',u3=-(dzPunch+epsContact))     
    mmd.boundaryConditions['Fix_Punch'].setValuesInStep(stepName='Step-3',u3=-(dzPunch+epsContact+x3Punch)) 
    mmd.boundaryConditions['Fix_Die'].setValuesInStep(stepName='Step-4',u3=zdeltaPunch+punchNoseRad/2.0) 
    mmd.boundaryConditions['Fix_Punch'].setValuesInStep(stepName='Step-4',u3=-(dzPunch+epsContact+x3Punch+zdeltaPunch+punchNoseRad/2.0))
    #mmd.boundaryConditions['Fix_Punch'].setValuesInStep(stepName='Step-4',u3=-(2*dzPunch+x3Punch+x3Punch2))      
    #
    ara.Set(edges=ara.instances[blankInstName].edges.findAt(((0.66*blankRad,0.0,zBlank),)),name='blankAxisRD')
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Symm_RD', 
        region=ara.sets['blankAxisRD'], 
        u1=UNSET, u2=SET, u3=UNSET, ur1=SET, ur2=UNSET, ur3=SET)
    ara.Set(edges=ara.instances[blankInstName].edges.findAt(((0.0,0.66*blankRad,zBlank),)),name='blankAxisTD')
    mmd.DisplacementBC(amplitude=UNSET, createStepName='Initial',
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name='Symm_TD', 
        region=ara.sets['blankAxisTD'], 
        u1=SET, u2=UNSET, u3=UNSET, ur1=UNSET, ur2=SET, ur3=SET)
    ara.Set(vertices=ara.instances[blankInstName].vertices.findAt(((zro,zro,zBlank),)),name='blankCenter')        
    #mmd.DisplacementBC(name='Fix_cBlank', createStepName='Step-1', 
    #region=ara.sets['blankCenter'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    #amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)
    #mmd.boundaryConditions['Fix_cBlank'].deactivate(stepName='Step-2')    
    #
    ##mmd.DisplacementBC(name='Fix_cBlank', createStepName='Step-3', 
    ##region=ara.sets['blankCenter'], u1=SET, u2=SET, u3=SET, ur1=SET, ur2=SET, ur3=SET, 
    ##amplitude=UNSET, fixed=ON, distributionType=UNIFORM, fieldName='', localCsys=None) 
    edges=ara.instances[blankInstName].edges
    edge=edges.findAt(((blankRad*mcos(mpi/6.0),blankRad*msin(mpi/6.0),zBlank),))
    mmd.ShellEdgeLoad(name='StretchZero', createStepName='Step-1', 
    region=ara.Surface(side1Edges=edge, name='SStretchZero'), 
    magnitude=-0.1, distributionType=UNIFORM, field='', localCsys=None)
    mmd.loads['StretchZero'].deactivate(stepName='Step-2')    
    faces=ara.instances[blankInstName].faces[0]
    edges=ara.instances[blankInstName].edges
    partRad=0.92*(punchRad-punchNoseRad)
    ara.PartitionFaceByCurvedPathEdgePoints(face=faces, 
        edge2=edges.findAt(coordinates=(0.66*blankRad,zro,zBlank)), point2=(partRad,zro,zBlank), 
        edge1=edges.findAt(coordinates=(zro,0.66*blankRad,zBlank)), point1=(zro,partRad,zBlank))
    faces=ara.instances[blankInstName].faces
    fac=faces.findAt(((0.8*blankRad*mcos(mpi/6.0),0.8*blankRad*msin(mpi/6.0),zBlank),))
    ara.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    fac=faces.findAt(((0.6*partRad*mcos(mpi/6.0),0.6*partRad*msin(mpi/6.0),zBlank),))
    ara.setMeshControls(elemShape=QUAD, regions=fac, technique=STRUCTURED)
    nsCirc=int((0.5*mpi*blankRad)/2)
    nsRadA=int((blankRad-partRad)/0.8)
    nsRadB=int((nsRadA*partRad/(blankRad-partRad))/2)
    #seedScaleC,seedScaleR=1.3,1.2;nsCirc=int(seedScaleC*nsCirc);nsRadA=int(seedScaleR*nsRadA);nsRadB=int(seedScaleR*nsRadB)
    seedScaleC,seedScaleR=aMAT['seedScaleC'],aMAT['seedScaleR'] ###1.12,1.15
    nsCirc=int(seedScaleC*nsCirc);nsRadA=int(seedScaleR*nsRadA);nsRadB=int(seedScaleR*nsRadB)
    ##printMsg('nsCirc = '+str(nsCirc)+'; nsRadA = '+str(nsRadA)+'; nsRadB = '+str(nsRadB))
    edges=ara.instances[blankInstName].edges
    edge=edges.findAt(((blankRad*mcos(mpi/6.0),blankRad*msin(mpi/6.0),zBlank),))
    ara.seedEdgeByNumber(constraint=FINER, edges=edge, number=nsCirc)
    for edge in edges:
        vp=edge.pointOn[0]
        if(vp[0]>0.0 and vp[1]>0.0 and msqrt(vp[0]**2+vp[1]**2)<0.9*blankRad):break
    ###edge=edges.findAt(coordinates=(partRad*mcos(mpi/6.0),partRad*msin(mpi/6.0),zBlank));print(edge)
    ara.seedEdgeByNumber(constraint=FINER, edges=(edges[edge.index],), number=nsCirc)
    edge=edges.findAt(((0.66*partRad,0.0,zBlank),))
    ara.seedEdgeByNumber(constraint=FINER, edges=edge, number=nsRadB)
    edge=edges.findAt(((0.0,0.66*partRad,zBlank),))
    ara.seedEdgeByNumber(constraint=FINER, edges=edge, number=nsRadB)  
    edge=edges.findAt(((0.9*blankRad,0.0,zBlank),))
    ara.seedEdgeByNumber(constraint=FINER, edges=edge, number=nsRadA)
    ara.setElementType(elemTypes=(mesh.ElemType(elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=ON, hourglassControl=ENHANCED),),
    regions=ara.Set(faces=faces,name='wholeBlank'))
    #ara.setElementType(elemTypes=(mesh.ElemType(elemCode=S4, elemLibrary=STANDARD, secondOrderAccuracy=ON, hourglassControl=DEFAULT),),
    #regions=ara.Set(faces=faces,name='wholeBlank'))
    ara.generateMesh(regions=faces.findAt(((0.8*blankRad*mcos(mpi/6.0),0.8*blankRad*msin(mpi/6.0),zBlank),)))
    ara.generateMesh(regions=faces.findAt(((0.6*partRad*mcos(mpi/6.0),0.6*partRad*msin(mpi/6.0),zBlank),)))    
    ara.regenerate()
    mdb.Job(atTime=None, contactPrint=OFF, description=JOB_DESCRIPTION, echoPrint=OFF, explicitPrecision=SINGLE, 
        getMemoryFromAnalysis=True, historyPrint=OFF, memory=90, memoryUnits=PERCENTAGE, 
        model=mmd, modelPrint=OFF, multiprocessingMode=DEFAULT, 
        name=jbName, nodalOutputPrecision=FULL, numCpus=1, numGPUs=0, 
        queue=None, resultsFormat=ODB, scratch=thisPATH, 
        type=ANALYSIS, userSubroutine='', waitHours=0, waitMinutes=0)
    vNodes,vnCircLabels,vnEdgeRDLabels,vnEdgeTDLabels=ara.instances[blankInstName].nodes,[],[],[]
    for node in vNodes:
        vp=node.coordinates        
        if(abs(msqrt(vp[0]*vp[0]+vp[1]*vp[1])-blankRad)<gTol):
            vnCircLabels.append(node.label)
    ncBlank=0
    for node in vNodes:
        vp=node.coordinates
        if(abs(vp[0])+abs(vp[1])<gTol/100.0):ncBlank=node.label;break
    if(ncBlank<0.1):printMsg('Warning: no center found')        
    ara.SetFromNodeLabels(nodeLabels=((blankInstName,tuple(vnCircLabels)),),name='nodeCirc')     
    ara.regenerate() 
    ##mdb.saveAs(modelDDname)    
    return mdb,vnCircLabels,vnEdgeRDLabels,vnEdgeTDLabels,ncBlank

###def runDDtests(punchTravel,dtInitial,dtMax,maxNInc,pSoft,deltaPunch=0,nCPUs=1):
def runDDtests():
    oschdir(thisPATH)
    ##JOB_NAME='JOB-CpDrw-Shell'
    JOB_NAME=thisDDname+'-JOB'
    pSoft=aMAT['pSoft']
    if(len(pSoft)!=3):printMsg('pSoft must be a 3-elem sequence');sysexit()
    deltaPunch=aMAT['deltaPunch']
    punchTravel=isoProfile(aMAT,scale=aMAT['scaleP'])
    dtInitial,dtMax,maxNInc,nCPUs=aMAT['dtInitial'],aMAT['dtMax'],aMAT['maxNInc'],aMAT['nCPUs']
    mdb,vnCircLabels,vnEdgeRDLabels,vnEdgeTDLabels,ncBlank=genDDabqModel(JOB_NAME,
    punchTravel,dtInitial,dtMax,maxNInc,pSoft,deltaPunch)
    mmd=mdb.models[mdbModelName]
    aa=mmd.rootAssembly
    staALL=figDirPlot+'cpdrwEndStatus.txt'
    staFF=open(staALL,'w');staFF.close()
    for ttest in vvTests:
        ff=figDirPlot+ttest['dfile'].strip('.txt')
        fileJob,fileScreen,fileStatus=ff+'_cpdrwJOB.txt',ff+'_cpdrwScreen.png',ff+'_cpdrwSTATUS.txt'
        fileProfile,fileStrain,filePunch=ff+'_cpdrwPROFILE.txt',ff+'_cpdrwSTRAIN.txt',ff+'_cpdrwPUNCHf.txt'
        #
        if(not aMAT['readOdb']):
            mmat=mmd.Material(description=aMAT['descrp'], name=aMAT['name'])
            polyParam=[]
            if(ttest['PolyN']):
                try:
                    zff=open(figDirData+ttest['dfile'],'r')
                    for line in zff:
                        polyParam.append(float(line.strip()))
                    zff.close()    
                except IOError:
                    printMsg('Input file for Poly_UMAT not found: '+ttest['dfile'])
                    printMsg('Calculations aborted');sysexit()
                mdb.jobs[JOB_NAME].setValues(userSubroutine=ttest['ufile'])    
                mmat.UserMaterial(mechanicalConstants=tuple(polyParam))
                mmat.Depvar(n=1)
                ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles,vcf=polyParam[5:])
            elif(ttest['UMAT']):
                mdb.jobs[JOB_NAME].setValues(userSubroutine=ttest['ufile'])
                hardP,matP=aMAT['hardParam'],paramHill(aMAT['rVals'])
                #tmp=(aMAT['E'],aMAT['nuPoisson'],hardP[0],hardP[1],hardP[2],matP[0],matP[1],matP[2],matP[3])
                #tmp=(len(tmp)*'{:.4f}, ').format(*tmp);printMsg(tmp);sysexit()
                mmat.UserMaterial(mechanicalConstants=(aMAT['E'],aMAT['nuPoisson'], 
                hardP[0],hardP[1],hardP[2],matP[0],matP[1],matP[2],matP[3]))
                mmat.Depvar(n=1)
                ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles)
            else:
                mmat.Elastic(table=((aMAT['E'], aMAT['nuPoisson']), ))
                epMax=0.651 ###max strain for *MATERIAL hardening curve 
                mmat.Plastic(table=hardFunc(aMAT['hardLaw'],aMAT['hardParam'],epMax))
                mmat.plastic.Potential(table=(abqHillRatios(aMAT['rVals']),))
                ##vStrains=vStrainsHill(aMAT['hardLaw'],aMAT['hardParam'],vAngles)            
            #else:
            #    pass
        ###
            mmd.HomogeneousShellSection(idealization=NO_IDEALIZATION, 
            integrationRule=SIMPSON, material=aMAT['name'], name='Section-1', 
            numIntPts=11, poissonDefinition=DEFAULT, preIntegrate=OFF, temperature=GRADIENT, 
            thickness=aMAT['hThick'], thicknessField='', thicknessModulus=None, 
            thicknessType=UNIFORM, useDensity=OFF)
            #aa=mmd.rootAssembly
            #pPart=mdb.parts[partName]
            pPart=aa.instances[blankInstName].part
            mmd.sections['Section-1'].TransverseShearShell(k11=transvK,k12=0.0, k22=transvK)
            pPart.SectionAssignment(offset=0.0, offsetField='', offsetType=MIDDLE_SURFACE, 
            region=pPart.sets['wholeBlankSection'], 
            sectionName='Section-1', thicknessAssignment=FROM_SECTION)
            if(ttest['UMAT']):    
                #mmd.fieldOutputRequests['F-Output-2'].setValues(variables=('S', 'LE', 'U', 'SDV'))
                mmd.fieldOutputRequests['F-Output-3'].setValues(variables=('S', 'LE', 'U', 'SDV'))
                mmd.fieldOutputRequests['F-Output-4'].setValues(variables=('S', 'LE', 'U', 'SDV'))
            else:
                #mmd.fieldOutputRequests['F-Output-2'].setValues(variables=('S', 'LE', 'U', 'PEEQ'))
                mmd.fieldOutputRequests['F-Output-3'].setValues(variables=('S', 'LE', 'U', 'PEEQ'))
                mmd.fieldOutputRequests['F-Output-4'].setValues(variables=('S', 'LE', 'U', 'PEEQ'))            
            aa.regenerate()
            ##
            vFileStress,vFileNode,vFileJob=[],[],[]
            #
            mdb.jobs[JOB_NAME].setValues(numCpus=nCPUs)
            printMsg('---------- Cup drawing with {}'.format(ttest['dfile']))      
            aa.regenerate()
            if(caeOnly):
                cwd=osgetcwd()
                mdb.jobs[JOB_NAME].setValues(scratch=cwd)
                if(osn=='nt'):
                    mdb.jobs[JOB_NAME].setValues(userSubroutine=cwd+'\\'+ttest['ufile'].strip(thisPATH))
                else:
                    mdb.jobs[JOB_NAME].setValues(userSubroutine=cwd+'/'+ttest['ufile'].strip(thisPATH))                
                mdb.saveAs(modelDDname);printMsg("CAE file created: "+modelDDname);sysexit()
            mdb.saveAs(modelDDname)    
            printMsg('Job submitted.........')            
        ##if(0):
            try:   
                mdb.jobs[JOB_NAME].submit(consistencyChecking=OFF)
            except:
                pass    
            mdb.jobs[JOB_NAME].waitForCompletion()
        ##mdb.saveAs(modelDDname)
        printMsg('------------------------------ Post-Proc Odb.....')
        try:
            ff2=open(JOB_NAME+'.dat','r')
            wTime=0 ##Wall Clock 
            for line in ff2:
                if('WALLCLOCK TIME' in line):
                    wTime+=int(line.strip().split('=')[1].strip())
                    printMsg(line.strip())
            ff2.close()
            ff2=open(fileJob,'w')
            ff2.write('Wall Clock (sec)\n')
            ff2.write(str(wTime)+'\n')
            ff2.close()
        except IOError as zerr:
            printMsg(zerr)
        try:
            ff2=open(JOB_NAME+'.sta','r')
            ff3=open(fileStatus,'w')
            wTime=''
            for line in ff2:
                ff3.write(line)
                if('ANALYSIS' in line):
                    wTime=line.strip()
            ff2.close();ff3.close()
            ff2=open(staALL,'a')
            ff2.write(ttest['dfile']+':  '+wTime+'\n')
            ff2.close()
        except IOError as zerr:
            printMsg(zerr)
        #        
        odb=openOdb(path=JOB_NAME+'.odb')
        #odb=openOdb(path='E:\\ABAQUStemp\\JOB-CpDrw-Test.odb')
        #for name in odb.rootAssembly.instances.keys():
        #    printMsg(name)
        podb=odb.rootAssembly.instances[blankInstName.upper()] 
        vnds,vndsA=podb.nodes,[]
        vnBlank=aa.instances[blankInstName].nodes
        if(len(odb.steps['Step-3'].frames)<0.1):
            printMsg('No frames in Step-3: post-proc aborted');odb.close();sysexit()
        if(aMAT['readOdb']==3):strFrame='Step-3'
        else:
            strFrame='Step-4'
            if(len(odb.steps['Step-4'].frames)<0.1):
                printMsg('No frames in Step-4: stepping-back post-processing to Step-3')
                strFrame='Step-3'
        ##        
        kframe=odb.steps[strFrame].frames[-1]
        #
        ff2=open(fileProfile,'w')
        ff2.write("X0 | Y0 | Z0 | UX | UY | UZ\n")
        for kk in vnCircLabels:
            #vndsA.append(vnds[kk-1])
            ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vnds[kk-1]).values
            bp=ss[0].dataDouble
            vpz=vnBlank[kk-1].coordinates            
            ff2.write("{}|{}|{}|{}|{}|{}\n".format(vpz[0],vpz[1],vpz[2],bp[0],bp[1],bp[2]))
        vpz=vnBlank[ncBlank-1].coordinates
        ss=kframe.fieldOutputs['U'].getSubset(position=NODAL,region=vnds[ncBlank-1]).values
        bp=ss[0].dataDouble
        ff2.write("{}|{}|{}|{}|{}|{}\n".format(vpz[0],vpz[1],vpz[2],bp[0],bp[1],bp[2]))        
        ff2.close()
        session.printOptions.setValues(vpDecorations=OFF)
        session.graphicsOptions.setValues(backgroundStyle=SOLID, backgroundColor='#FFFFFF')
        #session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(datumCoordSystems=ON)
        #session.viewports['Viewport: 1'].assemblyDisplay.geometryOptions.setValues(datumCoordSystems=OFF)
        session.printOptions.setValues(vpBackground=ON) 
        vv=session.viewports['Viewport: 1']
        vv.odbDisplay.setFrame(step=2, frame=len(odb.steps['Step-2'].frames)-1)
        vv.setValues(displayedObject=odb)
        vv.view.fitView(drawImmediately=True)
        vv.odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))
        leaf = dgo.LeafFromPartInstance(partInstanceName=(punchInstName.upper(), ))
        vv.odbDisplay.displayGroup.remove(leaf=leaf)
        leaf = dgo.LeafFromPartInstance(partInstanceName=(holderInstName.upper(), ))
        vv.odbDisplay.displayGroup.remove(leaf=leaf)
        vv.view.setValues(nearPlane=621.834,farPlane=1005.57, width=601.13, height=294.743, 
        cameraPosition=(-241.385,-527.764, 525.246), cameraUpVector=(0.373725, 0.812994, 0.44651), 
        cameraTarget=(-19.8594, -42.6011, 10.0452))
        vv.view.fitView(drawImmediately=True)
        #vv.view.setValues(nearPlane=492.05,farPlane=844.722, width=409.741, height=200.902, 
        #viewOffsetX=56.4518, viewOffsetY=6.71708)
        vv.odbDisplay.setPrimaryVariable(variableLabel='S', outputPosition=INTEGRATION_POINT, refinement=(INVARIANT,'Mises'))
        session.printToFile(fileName=fileScreen, format=PNG, canvasObjects=(vv,))
        session.graphicsOptions.setValues(backgroundStyle=GRADIENT,backgroundColor='#1B2D46')
        odb.close()
        printMsg('---------- '+ttest['dfile'].strip('.txt')+': DONE--------------------------------\n\n')
    printMsg('-----------ALL TESTS: DONE------------------------------------------------------')

    
def plotDDtests(bParam,inpData):
    global figDirData,figDirPlot
    ZsubDir=bParam['subDir']
    if(ZsubDir):
        if(osn=='nt'):figDirData+=ZsubDir+'\\'
        else:figDirData+=ZsubDir+'/'    
    if(not (ospath.exists(figDirData) and ospath.isdir(figDirData))):
        printMsg("The local data folder was not found")
        printMsg("Make sure a folder {} exists at your script location".format(figDirData))
        printMsg("Calculations aborted");sysexit()    
    if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
        printMsg("The local folder for saving reports and plots was not found")
        printMsg("Make sure a folder \'RESULTS\' exists at your script location")
        printMsg("Calculations aborted");sysexit()
    if(ZsubDir):
        if(osn=='nt'):figDirPlot+=ZsubDir+'\\'
        else:figDirPlot+=ZsubDir+'/'
        if(not(ospath.exists(figDirPlot) and ospath.isdir(figDirPlot))):
            try:
                osmaked(figDirPlot)
            except OSError as err:
                printMsg('plotDDtests:');printMsg(err);sysexit()
    fgg=plt.figure()            
    for ttest in inpData:
        ff=figDirPlot+ttest['inpFile'].strip('.txt')
        fileProfile=ff+'_cpdrwPROFILE'
        try:
            fHandle=open(fileProfile+'.txt','r')
        except IOError:
            printMsg('plotDDtests: file '+fileProfile+'.txt not found; Skipping');continue
        fHandle.readline()
        zData=[]
        for line in fHandle:
            line=line.strip().split('|')
            if(not line[0]):continue
            zData.append([float(kk) for kk in line])
        fHandle.close()
        cBlank=zData[-1][-1]
        nData=len(zData)-1
        zData=znp.array(zData[0:nData])
        zData=zData[zData[:,1].argsort()]
        zProfile,zTheta=znp.zeros(4*nData),znp.zeros(4*nData)
        ##print(zData.shape,zProfile.shape)
        zProfile[0:nData]=bParam['hThick']+zData[:,-1]-cBlank
        zTheta[0:nData]=180.0*(znp.arctan2(zData[:,1],zData[:,0])/znp.pi)
        zProfile[nData:2*nData]=zProfile[nData-1::-1]
        zTheta[nData:2*nData]=zTheta[nData-1]+zTheta[nData-1]-zTheta[nData-1::-1]
        zProfile[2*nData:]=zProfile[2*nData-1::-1]
        zTheta[2*nData:]=zTheta[2*nData-1]+zTheta[2*nData-1]-zTheta[2*nData-1::-1]
        fileTheta=open(ff+'_thetaPROFILE.txt','w')
        for kk in range(0,len(zTheta)):
            fileTheta.write("{}, {}\n".format(zTheta[kk],zProfile[kk]))
        fileTheta.close()            
        zData,kFile=[],0
        if(ttest['digProfile']):
            try:
                fHandle=open(figDirData+ttest['digProfile'][kFile],'r')
                for line in fHandle:
                    line=line.strip().split(ttest['digFileSep'][kFile])
                    if(not line[0]):continue
                    try:
                        zData.append([float(line[0]),float(line[1])])
                    except ValueError:
                        print(line)
                        printMsg('plotDDtests: Cannot convert to float data from file '+ttest['digProfile'][kFile])
                        break
                fHandle.close()                    
            except (IndexError,IOError):
                printMsg('plotDDtests: file '+ttest['digProfile'][0]+' not found; Skipping')
        ax=fgg.add_subplot()
        ax.plot(zTheta,zProfile,linestyle='-',linewidth=1,color='k',label='')
        if(zData):
            nData=len(zData)
            if(zData[-1][0]<100.0):
                zdProfile,zdTheta=znp.zeros(4*nData),znp.zeros(4*nData)
                zdTheta[:nData]=[yy[0] for yy in zData]
                zdProfile[:nData]=[yy[1] for yy in zData]
                zdProfile[nData:2*nData]=zdProfile[nData-1::-1]
                zdTheta[nData:2*nData]=zdTheta[nData-1]+zdTheta[nData-1]-zdTheta[nData-1::-1]
                zdProfile[2*nData:]=zdProfile[2*nData-1::-1]
                zdTheta[2*nData:]=zdTheta[2*nData-1]+zdTheta[2*nData-1]-zdTheta[2*nData-1::-1]
                ax.plot(zdTheta,zdProfile,linestyle='',marker='o',
                markersize=5,markeredgecolor='b',markerfacecolor='b',label='')
            else:    
                ax.plot([yy[0] for yy in zData],[yy[1] for yy in zData],linestyle='',marker='o',
                markersize=5,markeredgecolor='b',markerfacecolor='b',label='')
        ax.set_xlim(-2,362)
        ax.set_xticks([0,45,90,135,180,225,270,315,360])        
        fgg.savefig(fileProfile+'.png',dpi=300,bbox_inches='tight')
        fgg.clf()
    printMsg('--------------------All plots done')

def isoProfile(bParam,scale=1.075):
    hB=bParam['hThick']    
    blankDiam=bParam['blankDiam']
    dieOpeningDiam=bParam['dieOpeningDiam']
    dieShoulderRad=bParam['dieShoulderRad']
    punchDiam=bParam['punchDiam']
    punchNoseRad=bParam['punchNoseRad']
    blankRad=blankDiam/2.0
    dieOpeningRad=dieOpeningDiam/2.0
    punchRad=punchDiam/2.0
    x3Punch=punchRad-punchNoseRad    
    x3Punch=blankRad*blankRad-x3Punch*x3Punch-mpi*(x3Punch+(punchNoseRad+0.5*hB)/msqrt(2.0))*(punchNoseRad+hB)    
    x3Punch=(x3Punch/(dieOpeningRad+punchRad)+hB+punchNoseRad)    
    printMsg('Plane strain isotropic estimate of the cup height  = {}'.format(round(x3Punch,3)))
    x3Punch*=scale
    printMsg('Plane strain isotropic estimate of the punch travel  = scale*height = {}'.format(round(x3Punch,3)))
    return x3Punch

def zexec(inpData):
    if(abqSims):###---------STEP_1: simulation
        if(uaxSim):
            ##Run sims (or create cae only)
            runTests()
        else:
            runDDtests()
            #
    elif(pytPlot):###-------- STEP_2: post-processing
        if(uaxSim):
            ##readPlotData(bParam,zINP,ZHRATIO=HRATIO)
            plotTests()
        else:
            plotDDtests(aMAT,inpData)    
    else:
        msg='STEP_1 requires Abaqus\n'
        msg+='STEP_2 requires Python3[Numpy,Matplotlib]\nCalculations aborted'
        printMsg(msg);sysexit()

    
#if(__name__=='__main__'):
if(__name__=='__main__'):
    printMsg('\n'+3*'Use the driver script, Luke !\n')

