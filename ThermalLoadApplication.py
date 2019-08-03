##########################################################
##########################################################
## CONTAINS 3 MAIN FUNCTIONS: (FOR THERMAL ANALYSIS)    ##
##  1) APPLICATION OF USER SPECIFIC TEMP GRADIENT       ##
##  2) APPLICATION OF USER SPECIFIC HEAT FLUX           ##
##  3) UniformSurfaceFlux TYPE LOADING                  ##
##########################################################
##########################################################

from abaqus import *
from math import *
from abaqusConstants import *
import mesh
from part import *
from material import *
from section import *
from assembly import *
from step import *
from interaction import *
from load import *
from mesh import *
from job import *
from sketch import *
from visualization import *
import sys
import py_compile
import os
import string
import math
import regionToolset

import numpy

from ConstraintsUtilities import *

def applyTemperarureFromODBToCenterOfRVE(Model, historyOutputT):
    """
    applyTemperarureFromODBToCenterOfRVE:
    This function takes temperature history or XY data from ODB and apply that to center node of RVE.
    
    Inputs:
    Model-                          A model object from a Model Database
    
    historyOutputT -                temperature history
    """
    
    Amp_T1 = createAmplitudeObject(Model, 'Amp_T1', historyOutputT)
    
    T1 = 1.0
    
    applyTemperatureToRVE(Model, T1, Amp_T1.name)
    

def applyTemperatureToRVE(Model, T1, Amp_T1):
    """
    applyTemperatureToRVE:
    This function takes magnitude and amplitude of temperature and apply to RVE center.
    
    Inputs:
    Model-                          A model object from a Model Database
    
    T1 -                            magnitude
    
    Amp_T1 -                        amplitude
    
    centreNodeSet -                 center node set object 
    """
    
    stepNames = Model.steps.keys()
    centreNodeSet = findNodeAtCenterOfRVE(Model)    
    region1 = Model.rootAssembly.sets['Centre Node']    

    
    for stepName in stepNames:
        if ('LoadHist' in stepName):
            if (Amp_T1 != None):            
                Model.TemperatureBC(amplitude=Amp_T1, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=T1, name=
                    'Temp for-%s' %stepName, region=region1)
            else:
                Model.TemperatureBC(amplitude=UNSET, createStepName=
                    stepName, magnitude=T1, name='Temp for-%s' %stepName, region=region1) 

                    

def applyLoadsFromODBFluxDrivenPBC(Model, historyOutputQ1=None, historyOutputQ2=None, historyOutputQ3=None):
    """
    applyLoadsFromODBFluxDrivenPBC:
    this function will read flux history or XY data from ODB and apply loads according to them.
    
    Inputs:
    Model-                          A model object from a Model Database
    
    historyOutputQ1  -            history output object for 3 components
    """

    Amp_Q1 = createAmplitudeObject(Model, 'Amp_Q1', historyOutputQ1)
    Amp_Q2 = createAmplitudeObject(Model, 'Amp_Q2', historyOutputQ2)
    Amp_Q3 = createAmplitudeObject(Model, 'Amp_Q3', historyOutputQ3)

    
    Q1 = 1.0
    Q2 = 1.0
    Q3 = 1.0

    
    userFluxDrivenPBC(Model, Q1, Q2, Q3, Amp_Q1.name, Amp_Q2.name, Amp_Q3.name)


def applyLoadsFromODBUniformSurfaceFluxThermalBC(Model, historyOutputQ1=None, historyOutputQ2=None, historyOutputQ3=None):
    """
    applyLoadsFromODBUniformSurfaceFluxThermalBC:
    this function will read flux history or XY data from ODB and apply loads according to them for UniformSurfaceFlux BC.
    
    Inputs:
    Model-                          A model object from a Model Database
    
    historyOutputQ1  -            history output object for 3 components
    """

    Amp_Q1 = createAmplitudeObject(Model, 'Amp_Q1', historyOutputQ1)
    Amp_Q2 = createAmplitudeObject(Model, 'Amp_Q2', historyOutputQ2)
    Amp_Q3 = createAmplitudeObject(Model, 'Amp_Q3', historyOutputQ3)

    
    Q1 = 1.0
    Q2 = 1.0
    Q3 = 1.0

    
    UniformSurfaceFluxThermalBC(Model, Q1, Q2, Q3, Amp_Q1, Amp_Q2, Amp_Q3)
    


# NOT IN USE    
def userTempetatureGradientDrivenPBC(Model, dTX, dTY, dTZ, Amp_dTX, Amp_dTY, Amp_dTZ):
    """
    userTempetatureGradientDrivenPBC:
    This function defines user specified temperature gradient values to reference nodes of PBC.
    
    Inputs:
    Model-                          A model object from a Model Database

    dTX, dTY, dTZ-   temperature gradient components
    
    Amp_dTX, Amp_dTY, Amp_dTZ-       Amplitudes of temperature gradient components
    """
    """
    region1 = Model.rootAssembly.sets['RP-ThermGrad_x']
    region2 = Model.rootAssembly.sets['RP-ThermGrad_y']
    region3 = Model.rootAssembly.sets['RP-ThermGrad_z']

    stepNames = Model.steps.keys()

    isNotFirstStep = False
    isNotFirstRegularStep = False

    # Value of perturbation
    perturbation = 0.01
    #counter for name of the loads and boundary conditions
    k = 0
    
    for stepName in stepNames:
        if (('X' in stepName) and (isNotFirstRegularStep == False)):  
            if (dTX != None):
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=perturbation, name=
                    'SmallX-%d' %k, region=region1)
            else:
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=perturbation, name='SmallX-%d' %k, region=region1)            

        if (('X' in stepName) and (isNotFirstRegularStep == True)):
            if (dTX != None):
                Model.boundaryConditions['dTX-%d' %(k-1)].deactivate(stepName)
            if (dTY != None):
                Model.boundaryConditions['dTY-%d' %(k-1)].deactivate(stepName)
            if (dTZ != None):
                Model.boundaryConditions['dTZ-%d' %(k-1)].deactivate(stepName)                

            if (dTX != None):                
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=perturbation, name=
                    'SmallX-%d' %k, region=region1)
            else:
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=perturbation, name='SmallX-%d' %k, region=region1) 
        
        if ('Y' in stepName):
            if (dTY != None):        
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=perturbation, name=
                    'SmallY-%d' %k, region=region2)
            else:
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=perturbation, name='SmallY-%d' %k, region=region2)
                    
            if (dTX != None):                    
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=0.0, name=
                    'SmallX-%d' %k, region=region1)                
            else:
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=0.0, name='SmallX-%d' %k, region=region1)

            
        if ('Z' in stepName): 
            if (dTZ != None):        
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=perturbation, name=
                    'SmallZ-%d' %k, region=region3)
            else:
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=perturbation, name='SmallZ-%d' %k, region=region3)
                    
            if (dTY != None):
                Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                    distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=0.0, name=
                    'SmallY-%d' %k, region=region2)
            else:        
                Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                    stepName, magnitude=0.0, name='SmallY-%d' %k, region=region2)
                    
            isNotFirstStep = True

        if ('LoadHist' in stepName):
            if (isNotFirstStep == True):
                if (dTZ != None):            
                    Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=0.0, name=
                        'SmallZ-%d' %k, region=region3)
                else:
                    Model.ConcentratedHeatFlux(amplitude=UNSET, createStepName=
                        stepName, magnitude=0.0, name='SmallZ-%d' %k, region=region3)                
                    
            if (dTX != None):
                if (Amp_dTX != None):
                    Model.TemperatureBC(amplitude=Amp_dTX, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTX, name=
                        'dTX-%d' %k, region=region1)
                else:
                    Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTX, name=
                        'dTX-%d' %k, region=region1)      
                
        
            if (dTY != None):
                if (Amp_dTY != None):
                    Model.TemperatureBC(amplitude=Amp_dTY, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTY, name=
                        'dTY-%d' %k, region=region2)
                else:
                    Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTY, name=
                        'dTY-%d' %k, region=region2)        
        
            if (dTZ != None):
                if (Amp_dTZ != None):
                    Model.TemperatureBC(amplitude=Amp_dTZ, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTZ, name=
                        'dTZ-%d' %k, region=region3)
                else:
                    Model.TemperatureBC(amplitude=UNSET, createStepName=stepName, 
                        distributionType=UNIFORM, fieldName='', fixed=OFF, magnitude=dTZ, name=
                        'dTZ-%d' %k, region=region3)
            
            isNotFirstRegularStep = True
        k = k + 1            
        """    


    
def userFluxDrivenPBC(Model, Q1, Q2, Q3, Amp_Q1, Amp_Q2, Amp_Q3, Volume=None):
    """
    userFluxDrivenPBC:
    This function defines user specified flux values to reference nodes of PBC.
    
    Inputs:
    Model-                          A model object from a Model Database

    Q1, Q2, Q3-   flux components
    
    Amp_Q1, Amp_Q2, Amp_Q3-       Amplitudes of flux components
    
    Volume -    Volume of RVE
    """

    if (Volume == None):
        Volume = CalculateVolumeFromMDB(Model)
        
    stepNames = Model.steps.keys()
    
    region1 = Model.rootAssembly.sets['RP-ThermGrad_x']
    region2 = Model.rootAssembly.sets['RP-ThermGrad_y']
    region3 = Model.rootAssembly.sets['RP-ThermGrad_z']

    #default values
    Q1_conc = None
    Q2_conc = None
    Q3_conc = None
    
    if ((Q1 != None) and (Volume != None)):
        Q1_conc = Q1*Volume
    if ((Q2 != None) and (Volume != None)):    
        Q2_conc = Q2*Volume
    if ((Q3 != None) and (Volume != None)):
        Q3_conc = Q3*Volume


    allComponents = [Q1_conc, Q2_conc, Q3_conc]
    allComponentsNames = ['Q1','Q2','Q3']
    allAmplitudes = [Amp_Q1, Amp_Q2, Amp_Q3]
    
    #Create history loads in first history-type step (non-perturbation step)
    HistoryLoads = []
    for stepName in stepNames:
        if 'LoadHist' in stepName:  #This is a history step
            for component,componentName,amplitudeName in zip(allComponents,allComponentsNames,allAmplitudes):
                if not component == None: #If this component is constrained
                    #Determine whether to apply to the normal ref node or the shear ref node
                    if ('Q1' in componentName):
                        regionToApply = region1
                    if ('Q2' in componentName):
                        regionToApply = region2
                    if ('Q3' in componentName):
                        regionToApply = region3                        
                    #Set up the amplitude if it's specified
                    ampArg = UNSET if amplitudeName == None else amplitudeName
                    HistoryLoads.append( Model.ConcentratedHeatFlux(name='%s_' %stepName+componentName, 
                                                                    createStepName=stepName, 
                                                                    region=regionToApply,
                                                                    magnitude=component,
                                                                    amplitude=ampArg) )
            #We only want to set up loads for the first history step - they will run thru the whole analysis.
            #break
    # Now we need to deactivate the history step in every pertrubation step used for homogenization.
    IsAfterFirstHistoryStep = False
    for stepName in stepNames:
        if 'LoadHist' in stepName: 
            IsAfterFirstHistoryStep=True
        if 'Homog' in stepName and IsAfterFirstHistoryStep:
            for HistLoad in HistoryLoads:
                HistLoad.deactivate(stepName)




def UniformSurfaceFluxThermalBC(Model, Q1, Q2, Q3, Amp_Q1, Amp_Q2, Amp_Q3):
    """
    UniformSurfaceFluxThermalBC:
    This function defines user specified flux values to non periodic RVE. It converts far field flux values in to 
    surface flux and apply them on RVE surfaces.
    
    Inputs:
    Model-                          A model object from a Model Database

    Q1, Q2, Q3-   Far field flux components
    
    Amp_Q1, Amp_Q2, Amp_Q3-       Amplitudes of flux components
    """
    
    stepNames = Model.steps.keys()

    CreateInitialSets(Model)

    # Value of perturbation
    #perturbation = 0.01
    
    #creating surfaces for applying loads
    surf = GetSurfaceFromNodeSet(Model,'RVE_NegX')
    Model.rootAssembly.Surface( name='RVE_NegX', face1Elements=surf)
    RVE_NegY = GetSurfaceFromNodeSet(Model,'RVE_NegY')
    Model.rootAssembly.Surface( name='RVE_NegY', face1Elements=surf)
    RVE_NegZ = GetSurfaceFromNodeSet(Model,'RVE_NegZ')
    Model.rootAssembly.Surface( name='RVE_NegZ', face1Elements=surf)
    RVE_PosX = GetSurfaceFromNodeSet(Model,'RVE_PosX')
    Model.rootAssembly.Surface( name='RVE_PosX', face1Elements=surf)
    RVE_PosY = GetSurfaceFromNodeSet(Model,'RVE_PosY')
    Model.rootAssembly.Surface( name='RVE_PosY', face1Elements=surf)
    RVE_PosZ = GetSurfaceFromNodeSet(Model,'RVE_PosZ') 
    Model.rootAssembly.Surface( name='RVE_PosZ', face1Elements=surf)
    
    regionNegX = Model.rootAssembly.surfaces['RVE_NegX']
    regionNegY = Model.rootAssembly.surfaces['RVE_NegY']
    regionNegZ = Model.rootAssembly.surfaces['RVE_NegZ']
    regionPosX = Model.rootAssembly.surfaces['RVE_PosX']
    regionPosY = Model.rootAssembly.surfaces['RVE_PosY']
    regionPosZ = Model.rootAssembly.surfaces['RVE_PosZ']
 

    allComponents = [Q1, Q2, Q3]
    allComponentsNames = ['Q1','Q2','Q3']
    allAmplitudes = [Amp_Q1, Amp_Q2, Amp_Q3]
    
    #Create history loads in first history-type step (non-perturbation step)
    HistoryLoads = []
    for stepName in stepNames:
        if 'LoadHist' in stepName:  #This is a history step
            for component,componentName,amplitudeName in zip(allComponents,allComponentsNames,allAmplitudes):
                if not component == None: #If this component is constrained
                    #Determine whether to apply to the normal ref node or the shear ref node
                    if ('Q1' in componentName):
                        regionToApply1 = regionNegX
                        regionToApply2 = regionPosX
                    if ('Q2' in componentName):
                        regionToApply1 = regionNegY
                        regionToApply2 = regionPosY
                    if ('Q3' in componentName):
                        regionToApply1 = regionNegZ
                        regionToApply2 = regionPosZ
                    #Set up the amplitude if it's specified
                    ampArg = UNSET if amplitudeName == None else amplitudeName
                    HistoryLoads.append( Model.SurfaceHeatFlux(name='NegFace_%s_' %stepName+componentName, 
                                                               createStepName=stepName, 
                                                               region=regionToApply1,
                                                               magnitude=-component,
                                                               amplitude=ampArg) )
                    HistoryLoads.append( Model.SurfaceHeatFlux(name='PosFace_%s_' %stepName+componentName, 
                                                               createStepName=stepName, 
                                                               region=regionToApply2,
                                                               magnitude=component,
                                                               amplitude=ampArg) )                                                               
            #We only want to set up loads for the first history step - they will run thru the whole analysis.
            #break
    # Now we need to deactivate the history step in every pertrubation step used for homogenization.
    IsAfterFirstHistoryStep = False
    for stepName in stepNames:
        if 'LoadHist' in stepName: 
            IsAfterFirstHistoryStep=True
        if 'Homog' in stepName and IsAfterFirstHistoryStep:
            for HistLoad in HistoryLoads:
                HistLoad.deactivate(stepName)



def homogLoadCasesFluxDrivenPCB(Model):
    """
    homogLoadCasesFluxDrivenPCB:
    This function loads all perturbation steps found in model. 
    It applies a gradient to RP so that very small temperature difference will be created across opposite faces of RVE. 
    
    Inputs:
    Model-                          A model object from a Model Database
    """
 
    RVEDimensions = CreateInitialSets(Model)
    
    stepNames = Model.steps.keys()
    region1 = Model.rootAssembly.sets['RP-ThermGrad_x']
    region2 = Model.rootAssembly.sets['RP-ThermGrad_y']
    region3 = Model.rootAssembly.sets['RP-ThermGrad_z']    
    
    tempDifferenceAcrossWidth = 1.0
    
    IsAfterFirstHistoryStep = False 
    
    HistoryLoads = []
    for stepName in stepNames: 
        
        if (stepName == 'Initial'):
            continue
        
        if 'X' in stepName:
            regionToApply = region1
            perturbation = (tempDifferenceAcrossWidth/RVEDimensions[0])
        if 'Y' in stepName:
            regionToApply = region2
            perturbation = (tempDifferenceAcrossWidth/RVEDimensions[1])
        if 'Z' in stepName:
            regionToApply = region3
            perturbation = (tempDifferenceAcrossWidth/RVEDimensions[2])
 
        if (IsAfterFirstHistoryStep == True):
            for HistLoad in HistoryLoads:
                HistLoad.deactivate(stepName)

        IsAfterFirstHistoryStep = True
 
        if not ('LoadHist' in stepName):
            HistoryLoads.append( Model.TemperatureBC(name='perturbation_%s' %stepName, 
                                                     createStepName=stepName, 
                                                     region=regionToApply,
                                                     magnitude=perturbation,
                                                     distributionType=UNIFORM,
                                                     fieldName='',
                                                     fixed=OFF,
                                                     amplitude=UNSET) )   




if __name__=="__main__":
    #userFluxDrivenPBC(mdb.models['Model-1'], 0.01,0.2,None,None,None,None,200)  
    homogLoadCasesFluxDrivenPCB(mdb.models['Model-1'])
    
