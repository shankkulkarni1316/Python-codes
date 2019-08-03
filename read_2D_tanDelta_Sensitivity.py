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

W=12

allInclusions = [25]


sortie = open('Results2D.txt' , 'w')
sortie.write('\t Shank S. Kulkarni \n')
sortie.write('\t tan delta is \n')
sortie.write('\n')


for max_incl in allInclusions:
    for ik in range (25,50):
        for w in range (1,W):

            avg1S22 = 0
            odbname='Job-%d-%d-%d' %(max_incl,ik,w)
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
            isFirstIP = True
            td_avg = 0
            #Eequ_avg = 0
            #Eequ1_avg = 0
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

    
                #Eequ = (S22data1[point2]/E22data1[point2])
                #Eequ1 = (S22data1[point1]/E22data1[point1])
                #############################
                ## tan delta
                #############################                
                #print point1
                #print point2
                
                diff=point2-point1
                delta=((360*3*diff*pi)/(1000*180))
                td=math.tan(delta)
                td_avg += td
                
                #Eequ_avg += Eequ
                #Eequ1_avg += Eequ1

            td_avg_value = td_avg/(len(allIPs)-1)
            #print td_avg_value
            #Eequ_avg_value = Eequ_avg/(len(allIPs)-1)
            #Eequ1_avg_value = Eequ1_avg/(len(allIPs)-1)

            
            sortie.write('\t \t  %f ' %td_avg_value)
            #sortie.write('\t \t  %f ' %Eequ_avg_value)
            #sortie.write('\t \t  %f ' %Eequ1_avg_value)
            sortie.write('\n')  
        sortie.write('\n')        
sortie.close()    
