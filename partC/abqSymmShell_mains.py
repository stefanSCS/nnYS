
### Execute this script in a command shell with the command: 
# abaqus cae noGUI=abqSymmShell_mains.py 

import abqSymmShell as abqp

##------- input data block 

####  uniaxial tests simulations 
fName='abq_uniaxial_nnYS_AA6016T4.txt'

####  or cup drawing simulations 
#fName='abq_cpDraw_nnYS_AA6016T4.txt'

### specify the subDir of DATA for the input file 'fName'
subDir='AA6016T4'

####------------end of input data block

###---execution block
inpD=abqp.readInpData(fName,subDir)
abqp.zexec(inpD)


