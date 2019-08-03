##    ############################################################################
##    ##             Reading the ODB file 3D
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

PP = 51   # number of jobs
FF = 11  # change in frequency steps

sortie = open('Results3DAgain.txt' , 'w')
sortie.write('\t Shank S. Kulkarni \n')
sortie.write('\t tan delta is \n')
sortie.write('\n')

for w in range (1,FF):

    for q in range (1,PP):

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
        
            S33time = tipHistories.historyOutputs['S33'].data
            E33time = tipHistories.historyOutputs['E33'].data
            
            sizeS33 = numpy.shape(S33time)
            sizeE33 = numpy.shape(E33time)
            
            def column(matrix, i):
                return [row[i] for row in matrix]
                
                
            S33data1=column(S33time,1)
            S33data2=column(S33time,0)
            E33data1=column(E33time,1)
            E33data2=column(E33time,0)
        
            maxS33 = 0
            maxE33 = 0
            for i in range (333,620):
                if maxS33 < S33data1[i]:
                    maxS33 = S33data1[i]
                    point1=i
                    
                    
                if maxE33 < E33data1[i]:
                    maxE33 = E33data1[i]
                    point2=i
        

            
            #Eequ = (S33data1[point2]/E33data1[point2])
            #Eequ1 = (S33data1[point1]/E33data1[point1])
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
                #sortie.write('\t %f' %td)
                #sortie.write('\n')
        
        td_avg_value = td_avg/(count)
        #Eequ_avg_value = Eequ_avg/(len(allIPs)-1)
        #Eequ1_avg_value = Eequ1_avg/(len(allIPs)-1)
        
        closeOdb(odb)
        
        sortie.write('\t \t  %f ' %td_avg_value)
        #sortie.write('\t \t  %f ' %Eequ_avg_value)
        #sortie.write('\t \t  %f ' %Eequ1_avg_value)
        sortie.write('\n')  
    sortie.write('\n')        

sortie.close()    

