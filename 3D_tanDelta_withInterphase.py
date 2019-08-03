##############################################
##############################################
## 3D TAN DELTA VARIATION                   ##
## SHANK KULKARNI                           ##
## 09/16/2016                               ##
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

r = 0.03
def SphericalPartition(i,r):
    a=r*sin(pi/6)
    b=r*cos(pi/6)
    #Datum planes
    dp1=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=xx1[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=yy1[i], principalPlane=XZPLANE)
	
    # Create half circle to sweep
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0, transform=
        mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i],yy1[i],zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,r), point2=(r,0))

    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(r,0), point2=(0,-r))

    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0,r), point2=(0,-r))
    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=
        mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01,0.01), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
    
    # Create circle as path to sweep
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0,
        transform=mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i], yy1[i], zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-r, 0))

    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01,0.01, 
        0.01), )), sketch=mdb.models['Model-1'].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
		
    ## Create the spherical partition
    m= mdb.models['Model-1'].parts['MATRIX']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]-a,zz1[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]+a,zz1[i]+b), ),))	

def subPartition(i,r):
    a=r*sin(pi/6)
    b=r*cos(pi/6)
## Creating Datum Planes 1
    dp1=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=xx1[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-1'].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=yy1[i], principalPlane=XZPLANE)
	
## Creating Partition profile 2
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0, transform=
        mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i],yy1[i],zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,r), point2=(r,0))

    mdb.models['Model-1'].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(r,0), point2=(0,-r))

    mdb.models['Model-1'].sketches['__profile__'].Line(point1=(0,r), point2=(0,-r))
    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=
        mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01,0.01), )), sketch=mdb.models['Model-1'].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp1.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
    
## Creating circle for giving path to sweep
    mdb.models['Model-1'].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0,
        transform=mdb.models['Model-1'].parts['MATRIX'].MakeSketchTransform(sketchPlane=mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i], yy1[i], zz1[i])))
    mdb.models['Model-1'].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-1'].sketches['__profile__'])
    mdb.models['Model-1'].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-r, 0))

    mdb.models['Model-1'].parts['MATRIX'].PartitionCellBySketch(cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01,0.01, 
        0.01), )), sketch=mdb.models['Model-1'].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-1'].parts['MATRIX'].datums[dp2.id], sketchUpEdge=mdb.models['Model-1'].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
		
## Creating the spherical partition
    m= mdb.models['Model-1'].parts['MATRIX']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]-a,zz1[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]+a,zz1[i]+b), ),))	    
    
    
PP = 100    ### number of random points

dis=numpy.zeros(PP)

NO_INCL = 30

Q=51    # Number of iterations
W=11    # change of frequency
V=2    # change of radius or % volume


max_incl = 24     # set number of inclusions required




########################################################################
##                   Let's Create Materials 
########################################################################        
    
mdb.models['Model-1'].Material(name='Interphase')
mdb.models['Model-1'].materials['Interphase'].Density(table=((1.0e-9, ), ))

mdb.models['Model-1'].materials['Interphase'].Elastic(table=((
    0.9436270853, 0.49214), ))        
mdb.models['Model-1'].materials['Interphase'].Viscoelastic(domain=TIME, table=((
    0.001497038998, 0.0, 0.032), (0.005959550149, 0.0, 0.1), (0.033512934, 0.0, 0.316), 
    (0.188457662, 0.0, 1.0), (0.5959550149, 0.0, 3.162), (0.11890852, 0.0, 10.0), (
    0.029868669, 0.0, 31.623), (0.01189109, 0.909091, 100.0), (0.007502477683, 
    0.030303, 316.228), (0.001497038998, 0.0, 1000.0), (0.001011931615, 0.0, 
    3162.278), (0.0007163379279, 0.0, 10000.0), (0.0006841563571, 0.0, 31622.777), (
    0.0006531666964, 0.0, 100000.0)), time=PRONY)

mdb.models['Model-1'].Material(name='Inclusion')
mdb.models['Model-1'].materials['Inclusion'].Density(table=((2.47*1e-9, ), ))
mdb.models['Model-1'].materials['Inclusion'].Elastic(table=
    ((64890, 0.249), ))
    
mdb.models['Model-1'].Material(name='Matrix')
mdb.models['Model-1'].materials['Matrix'].Density(table=((1.0e-9, ), ))

mdb.models['Model-1'].materials['Matrix'].Elastic(table=((
    29.87551867, 0.493776), ))        
mdb.models['Model-1'].materials['Matrix'].Viscoelastic(domain=TIME, table=((
    0.00017617, 0.0, 3.0), (0.00099076, 0.0,10.0), (0.0055714, 0.0, 32.0), 
    (0.01762, 0.0, 100.0), (0.05571, 0.0, 316.0), (0.3133, 0.0, 1000.0), (
    0.55714, 0.0, 3162.0), (0.03133, 0.8333333, 10000.0), (0.0078698, 
    0.0, 31623.0), (0.003133, 0.0, 100000.0), (0.00099076, 0.0, 
    316228.0), (0.00031328, 0.0, 1000000.0), (0.00017618, 0.0, 3162278.0), (
    0.00009906, 0.0, 10000000.0)), time=PRONY)    


    
for q in range (1,Q):

    pow = -7.5 - math.log10(2*math.pi)
    
    
    
    ## LETS CREATE MATRIX
    mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=2.0)
    mdb.models['Model-1'].sketches['__profile__'].sketchOptions.setValues(
        decimalPlaces=4)
    mdb.models['Model-1'].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
        point2=(0.3, 0.3))
    mdb.models['Model-1'].Part(dimensionality=THREE_D, name='MATRIX', type=
        DEFORMABLE_BODY)
    mdb.models['Model-1'].parts['MATRIX'].BaseSolidExtrude(depth=0.3, sketch=
        mdb.models['Model-1'].sketches['__profile__'])
    del mdb.models['Model-1'].sketches['__profile__']        
  


        
#####################################################################################################################
## DEFINE SECTIONS
#####################################################################################################################

    # MATRIX
    mdb.models['Model-1'].HomogeneousSolidSection(material='Interphase', name='MATRIX', 
        thickness=None)
    # INCLUSION
    mdb.models['Model-1'].HomogeneousSolidSection(material='Inclusion', name='INCLUSION', 
        thickness=None)  
    # INTERPHASE    
    #mdb.models['Model-1'].HomogeneousSolidSection(material='Matrix', name='INTERPHASE', 
    #    thickness=None)

    num_incl = 0
    xx1 = []
    yy1 = []
    zz1 = []

    while (num_incl < max_incl):
        xx=random.uniform(0.04, 0.26)
        yy=random.uniform(0.04, 0.26)
        zz=random.uniform(0.04, 0.26)

        isPointIntersecting = False
        for j in range (0,len(xx1)):
    
    
            dis[j]=sqrt((xx-xx1[j])**2+(yy-yy1[j])**2+(zz-zz1[j])**2)

                
            if dis[j] < (2.2*r):

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
        #subPartition(i,1.1*r)
        i=i+1


    ### section assignment ################ 

    mdb.models['Model-1'].parts['MATRIX'].SectionAssignment(offset=0.0, 
        offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
        cells=mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), )), sectionName='MATRIX', thicknessAssignment=FROM_SECTION)
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

    ## Creating Instances
    
    mdb.models['Model-1'].rootAssembly.DatumCsysByDefault(CARTESIAN)
    mdb.models['Model-1'].rootAssembly.Instance(dependent=ON, name='MATRIX-1', 
        part=mdb.models['Model-1'].parts['MATRIX'])


########################################################################
##                   Let's Create Step 
########################################################################        

    mdb.models['Model-1'].ImplicitDynamicsStep(initialInc=(3*(10**((-pow)-3))), maxNumInc=100000, 
        name='Step-1', noStop=OFF, nohaf=OFF, previous='Initial', 
        timeIncrementationMethod=FIXED, timePeriod=(2*(10**(-pow))))
            
########################################################################
##                   Let's Create BC's and load
######################################################################## 
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-1', region=Region(
        faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
        ((0.0, 0.015, 0.075), ), )), u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)
    	
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-2', region=Region(
        faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
        ((0.015, 0.0, 0.075), ), )), u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)				
    
    mdb.models['Model-1'].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
        distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
        'BC-3', region=Region(
        faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
        ((0.075, 0.075, 0.0), ), )), u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, 
        ur2=UNSET, ur3=UNSET)




########################################################################
##                   Let's Mesh
########################################################################
    ## Creating Mesh on the Part
    
    
    mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(
        elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
        regions=(mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01,0.01,0.01), ), ), ))
        
    i=0	
    while(i<num_incl):
        mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
            elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
            regions=(mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), ), ))
        #mdb.models['Model-1'].parts['MATRIX'].setElementType(elemTypes=(ElemType(elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
        #    elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
        #    regions=(mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-1.05*r,zz1[i]), ), ((xx1[i],yy1[i]+1.05*r,zz1[i]), ), ), ))                
        i=i+1

    # Seed the part	 
    mdb.models['Model-1'].parts['MATRIX'].seedPart(deviationFactor=0.1, 
        minSizeFactor=0.1, size=0.015)
    # Increase the seed on edges    
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.075, 0.0, 0.0), 
        )), minSizeFactor=0.1, size=0.03)    
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.0, 0.075, 0.0), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.0, 0.0, 0.075), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.075, 0.3, 0.0), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.075, 0.0, 0.3), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.3, 0.075, 0.0), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.0, 0.075, 0.3), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.3, 0.0, 0.075), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.0, 0.3, 0.075), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.075, 0.3, 0.3), 
        )), minSizeFactor=0.1, size=0.03)
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.3, 0.075, 0.3), 
        )), minSizeFactor=0.1, size=0.03)            
    mdb.models['Model-1'].parts['MATRIX'].seedEdgeBySize(deviationFactor=0.1, 
        edges=mdb.models['Model-1'].parts['MATRIX'].edges.findAt(((0.3, 0.3, 0.075), 
        )), minSizeFactor=0.1, size=0.03)
        
    mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
        mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), ), technique=FREE)
    
    i=0	
    while(i<num_incl):
        mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
    	    mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), ), technique=FREE)
        #mdb.models['Model-1'].parts['MATRIX'].setMeshControls(elemShape=TET, regions=
    	#    mdb.models['Model-1'].parts['MATRIX'].cells.findAt(((xx1[i],yy1[i]-1.05*r,zz1[i]), ), ((xx1[i],yy1[i]+1.05*r,zz1[i]), ), ), technique=FREE)                
        i=i+1
    mdb.models['Model-1'].parts['MATRIX'].generateMesh()

    mdb.models['Model-1'].Tie(adjust=OFF, master=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].edges.findAt(
        ((0.3, 0.15, 0.0), ), )), name='Constraint-1', positionToleranceMethod=
        COMPUTED, slave=Region(
        side1Faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
        ((0.3, 0.15, 0.15), ), )), thickness=ON, tieRotations=OFF)
    mdb.models['Model-1'].Tie(adjust=OFF, master=Region(
        edges=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].edges.findAt(
        ((0.15, 0.3, 0.0), ), )), name='Constraint-2', positionToleranceMethod=
        COMPUTED, slave=Region(
        side1Faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
        ((0.15, 0.3, 0.15), ), )), thickness=ON, tieRotations=ON)

        
        
########################################################################
##                   Let's Create Sets
########################################################################	

    mdb.models['Model-1'].rootAssembly.Set(faces=
        mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(((
        0.075, 0.075, 0.0), )), name='Bottom_Face')
    	
    mdb.models['Model-1'].HistoryOutputRequest(createStepName='Step-1', frequency=1
        , name='H-Output-2', rebar=EXCLUDE, region=
        mdb.models['Model-1'].rootAssembly.sets['Bottom_Face'], sectionPoints=DEFAULT, 
        variables=('S33', 'E33'))

            
    for w in range (1,W):   

        w0=2*(math.pi)*10**pow ## CHANGING THE FREQUNCY WITH RESPECT TO pow
    
        #mdb.models['Model-1'].rootAssembly.regenerate()
        mdb.models['Model-1'].steps['Step-1'].setValues(initialInc=(3*(10**((-pow)-3))), timePeriod=
            (2*(10**(-pow))))            
            
        mdb.models['Model-1'].PeriodicAmplitude(a_0=0.0, data=((0.0, 1.0), ), 
            frequency=w0, name='Amp-1', start=0.0, timeSpan=STEP)
        mdb.models['Model-1'].DisplacementBC(amplitude='Amp-1', createStepName='Step-1', 
            distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
            'BC-4', region=Region(
            faces=mdb.models['Model-1'].rootAssembly.instances['MATRIX-1'].faces.findAt(
            ((0.075, 0.075, 0.3), ), )), u1=UNSET, u2=UNSET, u3=0.003, ur1=UNSET, 
            ur2=UNSET, ur3=UNSET)            
    ########################################################################
    ##                   Let's Create Jobs
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
        
        pow = pow + 1.0    ## Next frequency

# Save by skulka17 on Mon Apr 18 10:48:15 2016
