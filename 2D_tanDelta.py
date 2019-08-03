##############################################
##############################################
## 2D TAN DELTA VARIATION                   ##
## SHANK KULKARNI                           ##
## 08/22/2016                               ##
##############################################
##############################################

    


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
import math
import numpy
import os        # Operating system
import shutil    # copying or moving files

PP = 1000   ### number of random points
#xx=numpy.zeros(PP)
xx1=numpy.zeros(PP)
#yy=numpy.zeros(PP)
yy1=numpy.zeros(PP)
dis=numpy.zeros(PP)

NO_INCL = 30


rad=0.2
Q=2    # Number of iterations
W=19    # change of frequency
V=2    # change of radius or % volume


max_incl = 25      # set number of inclusions required



mdb.models['Model-1'].Material(name='VE')
mdb.models['Model-1'].materials['VE'].Density(table=((1.0e-9, ), ))

mdb.models['Model-1'].materials['VE'].Elastic(table=((
    0.9436270853, 0.49214), ))        
mdb.models['Model-1'].materials['VE'].Viscoelastic(domain=TIME, table=((
    0.001497038998, 0.0, 0.032), (0.005959550149, 0.0, 0.1), (0.033512934, 0.0, 0.316), 
    (0.188457662, 0.0, 1.0), (0.5959550149, 0.0, 3.162), (0.11890852, 0.0, 10.0), (
    0.029868669, 0.0, 31.623), (0.01189109, 0.909091, 100.0), (0.007502477683, 
    0.030303, 316.228), (0.001497038998, 0.0, 1000.0), (0.001011931615, 0.0, 
    3162.278), (0.0007163379279, 0.0, 10000.0), (0.0006841563571, 0.0, 31622.777), (
    0.0006531666964, 0.0, 100000.0)), time=PRONY)

mdb.models['Model-1'].Material(name='Elastic')
mdb.models['Model-1'].materials['Elastic'].Density(table=((2.47*1e-9, ), ))
mdb.models['Model-1'].materials['Elastic'].Elastic(table=
    ((64890, 0.249), ))

    
    

  
for q in range (0,25):
    pow = -7.0 - math.log10(2*math.pi)
    
    
    

    ##-------------------------------------------------------------------------------------------------        
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-2.0, 2.0), 
        point2=(2.0, -2.0))
    mdb.models['Model-1'].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['Part-1'].BaseShell(sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
        sheetSize=2, transform=
        mdb.models['Model-1'].parts['Part-1'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-1'].parts['Part-1'].faces[0], 
        sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
    mdb.models['Model-1'].parts['Part-1'].projectReferencesOntoSketch(filter=
        COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])

    num_incl = 0
    xx1 = []
    yy1 = []

    while (num_incl < max_incl):
        xx=random.uniform(-1.8, 1.8)
        yy=random.uniform(-1.8, 1.8)

        isPointIntersecting = False
        for j in range (0,len(xx1)):
    
    
            dis[j]=sqrt((xx-xx1[j])**2+(yy-yy1[j])**2)

                
            if dis[j] < (2.2*rad):

                isPointIntersecting = True
                break

        if (isPointIntersecting == False):
            xx1.append(xx)
            yy1.append(yy)
            num_incl = num_incl + 1        


    for i in range(max_incl):    

        mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(
            xx1[i], yy1[i]), point1=((xx1[i]-rad), yy1[i]))

        mdb.models['Model-1'].parts['Part-1'].PartitionFaceBySketch(faces=
            mdb.models['Model-1'].parts['Part-1'].faces.getSequenceFromMask(('[#1 ]', 
            ), ), sketch=mdb.models['Model-1'].sketches['__profile__'])
                    

########################################################################
##                   Let's Create Sections 
########################################################################

            
        
    mdb.models['Model-1'].HomogeneousSolidSection(material='VE', name='Matrix', 
        thickness=None)
    mdb.models['Model-1'].HomogeneousSolidSection(material='Elastic', name='Inclusion', 
        thickness=None)



    mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt(((1.9, 
        1.9, 0.0), (0.0, 0.0, 1.0)), )), sectionName='Matrix', 
        thicknessAssignment=FROM_SECTION)


    for i in range (num_incl):    

        mdb.models['Model-1'].parts['Part-1'].SectionAssignment(offset=0.2, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            faces=mdb.models['Model-1'].parts['Part-1'].faces.findAt((((xx1[i]), 
            (yy1[i]), 0.0), (0.0, 0.0, 1.0)), )), sectionName='Inclusion', 
            thicknessAssignment=FROM_SECTION)    

    
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
        part=mdb.models['Model-1'].parts['Part-1'])
    
    
    

    mdb.models['Model-1'].ImplicitDynamicsStep(initialInc=(3*(10**((-pow)-3))), maxNumInc=100000, 
        name='Step-1', noStop=OFF, nohaf=OFF, previous='Initial', 
        timeIncrementationMethod=FIXED, timePeriod=(3*(10**(-pow))))
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((-1.0, -2.0, 0.0), ), )), u1=UNSET, u2=0.0, ur3=UNSET)
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((-2.0, 1.0, 0.0), ), )), u1=0.0, u2=UNSET, ur3=UNSET)

    
    mdb.models['Model-1'].rootAssembly.seedPartInstance(deviationFactor=0.1, 
        minSizeFactor=0.1, regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ), size=0.1)
    
    for i in range (num_incl):
        mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
            elemCode=CPE8RH, elemLibrary=STANDARD), ElemType(elemCode=CPE6M, 
            elemLibrary=STANDARD)),regions=(
            mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((
            (xx1[i]), (yy1[i]), 0.0), )), ))        
    
    #ElemType(elemCode=CPE6M, 
    #    elemLibrary=STANDARD)),
    
    mdb.models['Model-1'].rootAssembly.setElementType(elemTypes=(ElemType(
        elemCode=CPE8RH, elemLibrary=STANDARD), ElemType(elemCode=CPE6M, 
        elemLibrary=STANDARD)),  regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].faces.findAt(((
        1.9, 1.9, 0.0), )), ))
    
    
    mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))    

    mdb.models['Model-1'].rootAssembly.generateMesh(regions=(
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'], ))
    mdb.models['Model-1'].Tie(adjust=ON, master=Region(
        vertices=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].vertices.findAt(
        ((2.0, -2.0, 0.0), ), )), name='Constraint-1', positionToleranceMethod=
        COMPUTED, slave=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(
        ((2.0, -1.0, 0.0), ), )), thickness=ON, tieRotations=ON)
        
        
    #CREATION OF SETS FOR ALL EDGES
    mdb.models['Model-1'].rootAssembly.Set(edges=
        mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(((
        -1.0, -2.0, 0.0), )), name='BottomEdge')


    #HISTORY OUTPUT REQUESTS FOR ALL SETS
    mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', frequency=1
        , name='H-Output-5', rebar=EXCLUDE, region=
        mdb.models['Model-1'].rootAssembly.sets['BottomEdge'], sectionPoints=DEFAULT, 
        variables=('S22', 'E22'))
    
    for w in range (1,W):  

        w0=2*(math.pi)*10**pow ## CHANGING THE FREQUNCY WITH RESPECT TO pow

        mdb.models['Model-1'].steps['Step-1'].setValues(initialInc=(3*(10**((-pow)-3))), timePeriod=
            (3*(10**(-pow)))) 

        mdb.models['Model-1'].PeriodicAmplitude(a_0=0.0, data=((0.0, 1.0), ), 
            frequency=w0, name='Amp-1', start=0.0, timeSpan=STEP)
        mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Step-1'
            , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-3', region=Region(
            edges=mdb.models['Model-1'].rootAssembly.instances['Part-1-1'].edges.findAt(
            ((1.0, 2.0, 0.0), ), )), u1=UNSET, u2=0.01, ur3=UNSET)            
            
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name='Job-%d-%d' %(w,q) , nodalOutputPrecision=SINGLE, 
            numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
            waitHours=0, waitMinutes=0)
        #mdb.jobs['Job-%d-%d' %(w,q)].writeInput()
        mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
        #mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()
        pow = pow + 0.5    ## Next frequency

# Save by skulka17 on Mon Apr 18 10:48:15 2016
