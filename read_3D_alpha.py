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


PP = 5   # number of jobs
FF = 5  # change in frequency steps

#frequency ranges
startHere = [700, 600, 550, 550]
endHere = [900, 850, 800, 750]

highEnd = [-12.0, -12.0, -12.0, -12.0]
lowEnd = [-6, -4, -2, -2]

sortie = open('Results3D_Alpha.txt' , 'w')
sortie.write('\t Shank S. Kulkarni \n')
sortie.write('\t Alpha \n')
sortie.write('\n')

received_waveMin=[]
received_waveMin_Good=[]

for w in range (1,FF):

    for q in range (1,PP):

        avg1S22 = 0
        odbname='Job-%d-%d' %(w,q) 
        
        
        path='./'
        myodbpath=path+odbname+'.odb'    
        odb=openOdb(myodbpath)
        count = 0

        
        allIPs = odb.steps['Step-1'].historyRegions.keys()

        isFirstIP = True
        
        for integrationPoint in allIPs:
            
            if (isFirstIP == True):
                isFirstIP = False
                continue
                
            tipHistories = odb.steps['Step-1'].historyRegions[integrationPoint]
        
            S33time = tipHistories.historyOutputs['S33'].data               
            
            sizeS33 = numpy.shape(S33time)
            
            def column(matrix, i):
                return [row[i] for row in matrix]
                
                
            S33data1=column(S33time,1)
            minS33 = 0
            applied = False
            for i in range (startHere[w-1],endHere[w-1]):          #------ FOR FREQUENCY 1,2,3,4
                if minS33 > S33data1[i]:
                    minS33 = S33data1[i]

                    
            if applied == False:
                received_waveMin.append(minS33)
                if ((minS33 > highEnd[w-1]) and (minS33 < lowEnd[w-1])):            
                    received_waveMin_Good.append(minS33)

        if (len(received_waveMin) != 0):
            Min_Avg_All = sum(received_waveMin) / float(len(received_waveMin))
            print Min_Avg_All
        
        if (len(received_waveMin_Good) != 0):
            Min_Avg_Good = sum(received_waveMin_Good) / float(len(received_waveMin_Good))
            print Min_Avg_Good

        # Calculate Alpha
        Alpha_1 = (log(10/abs(Min_Avg_All)))/1.5
        print Alpha_1
        
        Alpha_2 = (log(10/abs(Min_Avg_Good)))/1.5
        print Alpha_2        
        
        sortie.write('\t  %f \t\t %f' %(Alpha_1,Alpha_2))
        #for j in received_waveMin:        
        #    sortie.write('\n\t  %f ' %float(j))
        #    sortie.write('\n')  
           
        closeOdb(odb)
  
    sortie.write('\n Now frequency %f' %(w+1))        
    sortie.write('\n')
    
sortie.close()    
#workbook.close()
