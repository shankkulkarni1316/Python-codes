##    ############################################################################
##    ##             Submitting INP file
##    ############################################################################
# Shank Kulkarni
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from optimization import *
from job import *
from sketch import *
from visualization import *
from connectorBehavior import *


PP = 51   # number of iterations
FF = 12 # change in frequency steps


for w in range (1,FF):
    for q in range (1,PP):
        INPfilename='JobPBC-%d-%d' %(w,q)
        path='./'
        INPpath=path+INPfilename+'.inp'    

        mdb.ModelFromInputFile(inputFileName=INPpath
            , name='Model-1' )
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name='Job-%d-%d' %(w,q), nodalOutputPrecision=SINGLE, 
            numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
            waitHours=0, waitMinutes=0)
            
        mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
# Save by skulka17 on Tue Oct 04 08:58:07 2016
