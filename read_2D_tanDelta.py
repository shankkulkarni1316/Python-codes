##    ############################################################################
##    ##             Reading the ODB file
##    ############################################################################
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
import random
from array import *
from odbAccess import openOdb
import odbAccess
import math
import numpy    
import os        # Operating system
import shutil    # copying or moving files

PP = 2   # number of jobs
FF = 19  # change in frequency steps

alpha=numpy.zeros(PP)
td=numpy.zeros(PP)
td_final=numpy.zeros(FF)
alpha_final=numpy.zeros(FF)
deltax=6
S22=0
dl=1
pow=1

sortie = open('Results2D.txt' , 'w')
sortie.write('\t Shank S. Kulkarni \n')
sortie.write('\t tan delta is \n')
sortie.write('\n')


for v in range (1,2):
    for w in range (1,FF):

        for q in range (0,25):
            
            z = 0
            z = q+1
            avg1S22 = 0
            odbname='Job-%d-%d' %(w,q)
            path='./'
            myodbpath=path+odbname+'.odb'    
            odb=openOdb(myodbpath)

            count = 0        
            nameOfStep = 'Step-1'

            lastFrame = odb.steps[nameOfStep].frames[-1]

            Stress = lastFrame.fieldOutputs['S']

            ##    Number of elements  #### 
            numberOfElements = len(Stress.values)
            ##########################

            #for dl in range (1,2):
    
            allIPs = odb.steps['Step-1'].historyRegions.keys()
            #print allIPs
            isFirstIP = True
            td_avg = 0
            Eequ_avg = 0
            Eequ1_avg = 0
            for integrationPoint in allIPs:
                
                if (isFirstIP == True):
                    isFirstIP = False
                    continue
                    
                tipHistories = odb.steps['Step-1'].historyRegions[integrationPoint]
    
                S22time = tipHistories.historyOutputs['S22'].data
                E22time = tipHistories.historyOutputs['E22'].data
                
                sizeS22 = numpy.shape(S22time)
                sizeE22 = numpy.shape(E22time)
                
                def column(matrix, i):
                    return [row[i] for row in matrix]
                    
                    
                S22data1=column(S22time,1)
                S22data2=column(S22time,0)
                E22data1=column(E22time,1)
                E22data2=column(E22time,0)
    
                
                maxS22 = 0
                maxE22 = 0
                for i in range (666,1000):
                    if maxS22 < S22data1[i]:
                        maxS22 = S22data1[i]
                        point1=i
                        
                        
                    if maxE22 < E22data1[i]:
                        maxE22 = E22data1[i]
                        point2=i

                #############################
                ## tan delta
                #############################                
                
                diff=point2-point1
                delta=((360*3*diff*pi)/(1000*180))
                td=math.tan(delta)
                
                if ((td < 0) or (td > 4.0)):
                    continue
                else:
                
                    td_avg += td
                    count = count + 1
                #Eequ_avg += Eequ
                #Eequ1_avg += Eequ1

            td_avg_value = td_avg/(count)

            closeOdb(odb)
            sortie.write('\t \t  %f ' %td_avg_value)

            sortie.write('\n')  
        sortie.write('\n')        
sortie.close()    

#os.chdir(pathToFolder)