##############################################
##############################################
## 2D ALPHA                                 ##
## SHANK KULKARNI                           ##
## 11/15/2016                               ##
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




dis=numpy.zeros(1000)
rad=0.15
Q=21    # Number of iterations
W=2    # change of frequency
V=2    # change of radius or % volume


max_incl = 12     # set number of inclusions required

pow = 6.0 
Nf=1.0

r = 0.15
def SphericalPartition(i,r):
    a=r*sin(pi/6)
    b=r*cos(pi/6)
    #Datum planes
    dp1=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=xx1[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=yy1[i], principalPlane=XZPLANE)
	
    # Create half circle to sweep
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=2.0, name='__profile__', sheetSize=20.0, transform=
        mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 1.0), ), sketchOrientation=TOP, origin=(xx1[i],yy1[i],zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,r), point2=(r,0))

    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(r,0), point2=(0,-r))

    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0,r), point2=(0,-r))
    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=
        mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01,0.01), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 1.0), ))
    
    # Create circle as path to sweep
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=2.0, name='__profile__', sheetSize=20.0,
        transform=mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 1.0), ), sketchOrientation=TOP, origin=(xx1[i], yy1[i], zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-r, 0))

    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.02,0.02, 
        0.02), )), sketch=mdb.models['Model-1'].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 1.0), ))
		
    ## Create the spherical partition
    m= mdb.models['Model-1'].parts['MATRIX']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.02, 0.02, 0.02),),edges=(m.edges.findAt((xx1[i], yy1[i]-a,zz1[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.02, 0.02, 0.02),),edges=(m.edges.findAt((xx1[i], yy1[i]+a,zz1[i]+b), ),))	

for w in range (1,W):
    
    for q in range (1,Q):
        
        w0=2*(math.pi)*10**pow ## CHANGING THE FREQUNCY WITH RESPECT TO pow
        #w0=10**6
        
        # Visco-elastic matrix material
        mdb.models['Model-1'].Material(name='VE')
        mdb.models['Model-1'].materials['VE'].Density(table=((1.18e-09, ), ))
        mdb.models['Model-1'].materials['VE'].Elastic(moduli=INSTANTANEOUS, table=((
            4061.12, 0.3706), ))        
        mdb.models['Model-1'].materials['VE'].Viscoelastic(domain=TIME,  table=((
            0.0738, 0.0, 463.4), (0.147, 0.0, 0.06407), (0.3134, 0.0, 0.0001163), (
            0.3786, 0.0, 7.321e-07)), time=PRONY)

        # Elastic inclusions material    
        mdb.models['Model-1'].Material(name='Elastic')
        mdb.models['Model-1'].materials['Elastic'].Density(table=((2.47e-09, ), ))
        mdb.models['Model-1'].materials['Elastic'].Elastic(table=
            ((64890.0, 0.249), ))
##--------------------------------------------------------------------------------------------------- 
        ## LETS CREATE MATRIX
        mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)
        mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
            point2=(1.0, 1.0))
        mdb.models['Model-1'].Part(dimensionality=THREE_D, name='MATRIX', type=
            DEFORMABLE_BODY)
        mdb.models['Model-1'].parts['MATRIX'].BaseSolidExtrude(depth=2.0, sketch=
            mdb.models['Model-1'].sketches['__profile__'])    
        del mdb.models['Model-1'].sketches['__profile__']

        #mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.1, name='__profile__', 
        #    sheetSize=5.0, transform=
        #    mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(
        #    sketchPlane=mdb.models['Model-1'].parts['MATRIX'].faces.findAt((0.5, 
        #    0.5, 2.0), ), sketchPlaneSide=SIDE1, 
        #    sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((1.0, 0.5, 
        #    2.0), ), sketchOrientation=RIGHT, origin=(0.5, 0.5, 2.0)))
        #mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=
        #    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
        #mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(-0.12, 0.12), 
        #    point2=(0.12, -0.12))
        #mdb.models['Model-1'].parts['MATRIX'].PartitionFaceBySketch(faces=
        #    mdb.models['Model-1'].parts['MATRIX'].faces.findAt(((0.5, 0.5, 
        #    2.0), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
        #    sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((1.0, 0.5, 
        #    2.0), ))



    
        #################################################################################################
        ##   LET'S DEFINE SECTIONS
        #################################################################################################
        # MATRIX
        mdb.models['Model-1'].HomogeneousSolidSection(material='VE', name='MATRIX', 
            thickness=None)
        # INCLUSION
        mdb.models['Model-1'].HomogeneousSolidSection(material='Elastic', name='INCLUSION', 
            thickness=None)        

        #################################################################################################
        ## LET'S CREATE INCLUSIONS
        #################################################################################################            
        num_incl = 0
        xx1 = []
        yy1 = []
        zz1 = []

        while (num_incl < max_incl):
            xx=random.uniform((0+r), (1.0-r))
            yy=random.uniform((0+r), (1.0-r))
            zz=random.uniform((0.5+r), (2.0-r))
            
            isPointIntersecting = False
            for j in range (0,len(xx1)):
        
        
                dis[j]=sqrt((xx-xx1[j])**2+(yy-yy1[j])**2+(zz-zz1[j])**2)

                    
                if dis[j] < (2.2*rad):

                    isPointIntersecting = True
                    break

            if (isPointIntersecting == False):
                xx1.append(xx)
                yy1.append(yy)
                zz1.append(zz)
                num_incl = num_incl + 1        
        i=0
        while(i<num_incl):
            SphericalPartition(i,r)
            i=i+1
            
            
        #plane1=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=0.5, 
        #    principalPlane=XYPLANE)
        #mdb.models['Model-1'].parts['MATRIX'].PartitionCellByDatumPlane(cells=
        #    mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), )
        #   , datumPlane=mdb.models['Model-1'].parts['MATRIX'].datums[plane1.id])
           
        #mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.1, name='__profile__', 
        #    sheetSize=5.0, transform=
        #    mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(
        #    sketchPlane=mdb.models['Model-1'].parts['MATRIX'].faces.findAt((0.5, 
        #    0.5, 0.5), ), sketchPlaneSide=SIDE1, 
        #    sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((1.0, 0.25, 
        #    0.5), ), sketchOrientation=RIGHT, origin=(0.5, 0.5, 0.5)))
        #mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=
        #    COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
        #mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(
        #    -0.12, 0.12), point2=(0.12, 
        #    -0.12))
        #mdb.models['Model-1'].parts['MATRIX'].PartitionFaceBySketch(faces=
        #    mdb.models['Model-1'].parts['MATRIX'].faces.findAt(((0.5, 0.5, 
        #    0.5), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
        #    sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((1.0, 0.25, 
        #    0.5), ))
            
        #################################################################################################
        ##   LET'S ASSIGN SECTIONS
        #################################################################################################

        mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), )), sectionName='MATRIX', thicknessAssignment=FROM_SECTION)
        mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 1.99), ), )), sectionName='MATRIX', thicknessAssignment=FROM_SECTION)            
        i=0
        while(i<num_incl):
            mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]-0.2*r, zz1[i]+0.2*r), ), )), sectionName='INCLUSION', thicknessAssignment=FROM_SECTION)
            mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]+0.2*r, zz1[i]+0.2*r), ), )), sectionName='INCLUSION', thicknessAssignment=FROM_SECTION)
            #mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
            #    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            #    cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]-1.05*r, zz1[i]), ), )), sectionName='INTERPHASE', thicknessAssignment=FROM_SECTION)
            #mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
            #    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            #    cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]+1.05*r, zz1[i]), ), )), sectionName='INTERPHASE', thicknessAssignment=FROM_SECTION)                
            i=i+1        

        #################################################################################################
        ##   LET'S MESH THE PART 
        ################################################################################################# 
        mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
            mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), )
            , sizeGrowthRate=1.0, technique=FREE)        
        mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(
            elemCode=C3D20, elemLibrary=STANDARD,  
            kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
            distortionControl=DEFAULT), ElemType(elemCode=C3D15, elemLibrary=STANDARD), 
            ElemType(elemCode=C3D10, elemLibrary=STANDARD)), regions=(
            mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), ) 
            , ))
        mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
            mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 1.99), ), )
            , sizeGrowthRate=1.0, technique=FREE)        
        mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(
            elemCode=C3D20, elemLibrary=STANDARD,  
            kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
            distortionControl=DEFAULT), ElemType(elemCode=C3D15, elemLibrary=STANDARD), 
            ElemType(elemCode=C3D10, elemLibrary=STANDARD)), regions=(
            mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 1.99), ), ) 
            , ))
            
        for i in range (num_incl):
            mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
                mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), )
                , sizeGrowthRate=1.0, technique=FREE)        
            mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(
                elemCode=C3D20, elemLibrary=STANDARD,  
                kinematicSplit=AVERAGE_STRAIN, hourglassControl=DEFAULT, 
                distortionControl=DEFAULT), ElemType(elemCode=C3D15, elemLibrary=STANDARD), 
                ElemType(elemCode=C3D10, elemLibrary=STANDARD)), regions=(
                mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), ) 
                , ))             
            
        mdb.models['Model-1'].parts['MATRIX'].seedPart(deviationFactor=0.1, 
            minSizeFactor=0.1, size=0.05)
        mdb.models['Model-1'].parts['MATRIX'].generateMesh()
        
        #################################################################################################
        ##   LET'S CREATE INCTANCE 
        #################################################################################################        
        mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
        mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='MATRIX-1', 
            part=mdb.models['Model-1'].parts['MATRIX'])
        
        #################################################################################################
        ##   LET'S CREATE STEP
        #################################################################################################        
        mdb.models['Model-1'].ImplicitDynamicsStep(initialInc=(1*(10**((-6)-3))), maxNumInc=1000, 
            name='Step-1', noStop=OFF, nohaf=OFF, previous='Initial', 
            timeIncrementationMethod=FIXED, timePeriod=(1*(10**(-6))))
        
        
        
        ########################################################################
        ##                   Let's Create BC's and load
        ########################################################################   
        mdb.models['Model-1'].PeriodicAmplitude(a_0=0.0, data=((0.0, 1.0), ), 
            frequency=Nf*w0, name='Amp-1', start=0.0, timeSpan=STEP)

        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-1', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.5, 0.0, 1.0), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-2', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.0, 0.5, 1.0), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-3', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.5, 1.0, 1.0), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-4', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((1.0, 0.5, 1.0), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-5', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.5, 0.0, 0.25), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-6', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.0, 0.5, 0.25), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-7', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.5, 1.0, 0.25), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
        mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-8', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((1.0, 0.5, 0.25), ), )), u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=
            UNSET)
            
        mdb.models['Model-1'].Pressure(amplitude='Amp-1', createStepName='Step-1', 
            distributionType=UNIFORM, field='', magnitude=10.0, name='Load-1', region=
            Region(
            side1Faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.5, 0.5, 2.0), ), )))
            
        ########################################################################
        ## LET'S CREATE SETS AND HISTORY PUTPUT REQUESETS
        ########################################################################    
        #mdb.models['Model-1'].rootAssembly.Set(faces=
        #    mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(((
        #    0.5, 0.5, 2.0), )), name='Applied wave')
        mdb.models['Model-1'].rootAssembly.Set(faces=
            mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(((
            0.5, 0.5, 0.5), )), name='Received wave')
            
        #mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', frequency=1
        #    , name='H-Output-2', rebar=EXCLUDE, region=
        #    mdb.models['Model-1'].rootAssembly.sets['Applied wave'], sectionPoints=DEFAULT, 
        #    variables=('S33', 'E33'))
        mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', frequency=1
            , name='H-Output-3', rebar=EXCLUDE, region=
            mdb.models['Model-1'].rootAssembly.sets['Received wave'], sectionPoints=DEFAULT, 
            variables=('S33'))
        ########################################################################
        ## LET'S CREATE JOB
        ########################################################################         
        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-1', modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name='Job-%d-%d' %(w,q) , nodalOutputPrecision=SINGLE, 
            numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
            waitHours=0, waitMinutes=0)
        #mdb.jobs['Job-%d-%d' %(w,q)].writeInput()
        mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
        #mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()
    Nf = Nf + 1.0    ## Next frequency
    