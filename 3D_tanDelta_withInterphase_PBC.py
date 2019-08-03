##############################################
##############################################
## 3D TAN DELTA VARIATION                   ##
## SHANK KULKARNI                           ##
## 09/16/2016                               ##
##############################################
##############################################

    

import assembly
import mesh
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

#All useful functions
def PinTheNode(Model, NodeToPin):
    """
    PinTheNode:
    This function will fix the node which is given to this function.
    Function will create PinnedBC for initial step and that will not get deactivated ever.
    
    Inputs:
    Model - A model object from a Model Database
    
    NodeToPin -  node set name as a string   
    """
    
    stepNames = Model.steps.keys()

    regionPin = Model.rootAssembly.sets[NodeToPin]

    Model.PinnedBC(createStepName='%s' %stepNames[0], 
         name='The Pin' , region=regionPin)

def GenerateSortedMeshNodeArray(UnsortedMeshNodesArray):
    """
    Creates a MeshNodeArray that is sorted by (x,y,z) of the nodes.
    TODO - This needs to also sort based on the materials or section that a node is connected to
           to properly handle the case of coincident nodes on either side of a cohesive element
    
    Inputs:
    Model -                     Model object from mdb
    UnsortedMeshNodes -         Reference to iterable of unsorted nodes
    
    Output:
    Returns the created MeshNodesArray
    """
    
    # We need to handle random floating point roundoff carefully, since this can be different 
    # for each surfaces and will impact the way that they sort.  We'll create a comparator for
    # the sorting that should take care of the roundoff.  We only have to do this once per pair
    # so create it in this loop and not the next one.
    
    #First, get the maximum range for any coordinate in the nodes so that we can specify a reasonable tolerance automatically
    coords = [ n.coordinates for n in UnsortedMeshNodesArray ]
    TransCoords = zip(*coords)
    #Get ranges for x, y, and Z
    Ranges = [ max(a) - min(a) for a in TransCoords ]
    #Biggest range for any coordinate
    maxrange = max(Ranges)
    #Set a reasonable value for the tolerance relative to the maximum range for any coordinate
    NodeSetSortTolerance = maxrange * 1e-8
    def CompareNodeCoordWithTolerance(a, b):
        """Compares coordinates of nodes a and b using a tolerance"""
        CoordA = tuple(a.coordinates)
        CoordB = tuple(b.coordinates)
        #Check X. If equal, check Y. If equal, check Z.
        # -1 for a is first, 0 for Equal, 1 for b is first
        for (A,B) in zip(CoordA, CoordB):
            if A+NodeSetSortTolerance < B: #a before b beyond tolerance
                return -1
            elif B+NodeSetSortTolerance < A: #b before a beyond tolerance
                return 1
            # otherwise equal within tolerance.
            # go to next coordinate
        #Have gone through all coords and they are equal within tolerance.
        #This is probably bad because you have coincident nodes.
        return 0
    
    SortedNodesList = sorted(list(UnsortedMeshNodesArray), cmp=CompareNodeCoordWithTolerance)
    
    return mesh.MeshNodeArray(SortedNodesList)

def WriteMeshNodeArrayToUnsortedNodeSetDef(Nodes, Name):
    """
    WriteMeshNodeArrayToUnsortedNodeSetDef:
    Writes a string that will define *NSET, UNSORTED using the provided MeshNodesArray.
    This won't work for flattened input.
    
    Inputs:
    Nodes -             Some sequence of MeshNode objects that will be in the NSET
    Name  -             The name to be used for the NSET
    
    Output:
    A string that can be written to an input file defining an NSET.
    """
    
    Result = "*NSET, NSET=" + Name + ", UNSORTED\n"
    Result += '\n'.join([', '.join([n.instanceName + '.' + str(n.label) for n in rownodes]) for rownodes in GetInChunks(Nodes,8)])
    return Result    

# This is where the periodic equation gets defined.  This is the core kernel of this whole exercise.
def DefinePeriodicEquation3D( MasterNodesetName,
                              SlaveNodesetName,
                              RP_vertex ):
    """
    DefinePeriodicEquation3D:
    Creates equation objects in the model that define a 2D periodic relationship between the 
    nodesets with names MasterNodesetName and SlaveNodesetName.
    
    Inputs:
    MasterNodesetName -    The master nodeset (single node)
    SlaveNodesetName -     The slave nodeset (single node)  
    RP_Vertex -            reference point
    
    Outputs:
    String defining all equations to pair nodesets
    """
    EqnDefs = []
    #remember that Python is zero-based, so add 1 whenever using i or j to define DOFs
    
    for i in range(3):
        EqTerms = [ (-1.0,SlaveNodesetName,i+1), 
                    ( 1.0,MasterNodesetName,i+1) ]

        if (RP_vertex == "Vertex_NegXNegYPosZ"):
            EqTerms.append( (-1.0,RP_vertex,i+1) ) 
        else:
            EqTerms.append( (1.0,RP_vertex,i+1) )            
        EqnDefs.append(WriteEquationTerms(EqTerms))
    return '\n'.join(EqnDefs)

def WriteEquationTerms(Terms):
    """
    WriteEquationTerms:
    Creates string defining terms of a *EQUATION.
    
    Inputs:
    Terms -          A list of tuples defining the terms, each tuple being of the form (coefficient, set name, DOF)
    
    Output:
    A string that defines the terms of an equation (can be placed under *EQUATION)
    """
    
    Result = str(len(Terms))+'\n'
    Result += '\n'.join([', '.join([', '.join( [term[1], str(term[2]), '%e'%(term[0])] ) for term in rowterms]) for rowterms in GetInChunks(Terms,4)])
    return Result 
    
def GetInChunks(seq, n):
    """ Provides n-sized chunks from seq """
    
    for i in xrange(0, len(seq), n):
        yield seq[i:i + n]

def AddEquationsAndUnsortedNsetStringToKeywords(Model,KeywordsString):
    """
    AddEquationsAndUnsortedNsetStringToKeywords:
    Adds a string to the model's keywords via the keyword editor to define the equations and unsorted nodesets
    
    Inputs:
    Model -             The model object from the mdb that will be modificed
    KeywordsString -    String containing the equation and unsorted stes keywords
    """
    
    Model.keywordBlock.synchVersions(storeNodesAndElements=False)
    for i,block in enumerate(Model.keywordBlock.sieBlocks):
        KeyWord = block.partition(',')[0].lower().replace(' ','').replace('\t','')
        if KeyWord == "*endassembly" or KeyWord=="*step":
            Model.keywordBlock.insert(position=i-1,text=KeywordsString)
            return 
            

r = 0.03
def SphericalPartition(i,r):
    a=r*sin(pi/6)
    b=r*cos(pi/6)
    #Datum planes
    dp1=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=xx1[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=yy1[i], principalPlane=XZPLANE)
	
    # Create half circle to sweep
    mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0, transform=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i],yy1[i],zz1[i])))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,r), point2=(r,0))

    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(r,0), point2=(0,-r))

    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].Line(point1=(0,r), point2=(0,-r))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].PartitionCellBySketch(cells=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((0.01, 0.01,0.01), )), sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp1.id], sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
    
    # Create circle as path to sweep
    mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0,
        transform=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].MakeSketchTransform(sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i], yy1[i], zz1[i])))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-r, 0))

    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].PartitionCellBySketch(cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((0.01,0.01, 
        0.01), )), sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp2.id], sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
		
    ## Create the spherical partition
    m= mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]-a,zz1[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]+a,zz1[i]+b), ),))	

def subPartition(i,r):
    a=r*sin(pi/6)
    b=r*cos(pi/6)
## Creating Datum Planes 1
    dp1=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=xx1[i], principalPlane=YZPLANE)
    dp2=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].DatumPlaneByPrincipalPlane(offset=yy1[i], principalPlane=XZPLANE)
	
## Creating Partition profile 2
    mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0, transform=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].MakeSketchTransform(
        sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp1.id], sketchPlaneSide=SIDE1, 
        sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i],yy1[i],zz1[i])))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES, sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(0,r), point2=(r,0))

    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].ArcByCenterEnds(center=(0,0), direction=CLOCKWISE, point1=(r,0), point2=(0,-r))

    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].Line(point1=(0,r), point2=(0,-r))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].PartitionCellBySketch(cells=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((0.01, 0.01,0.01), )), sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'], 
        sketchOrientation=TOP, sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp1.id], sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
    
## Creating circle for giving path to sweep
    mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(gridSpacing=0.2, name='__profile__', sheetSize=2.0,
        transform=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].MakeSketchTransform(sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp2.id],
        sketchPlaneSide=SIDE1,sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ), sketchOrientation=TOP, origin=(xx1[i], yy1[i], zz1[i])))
    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].projectReferencesOntoSketch(filter=COPLANAR_EDGES,sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
    mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].CircleByCenterPerimeter(center=(0, 0), point1=(-r, 0))

    mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].PartitionCellBySketch(cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((0.01,0.01, 
        0.01), )), sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'],sketchOrientation=TOP, sketchPlane=
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].datums[dp2.id], sketchUpEdge=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].edges.findAt((0.0, 0.0, 0.075), ))
		
## Creating the spherical partition
    m= mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX']
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]-a,zz1[i]+b), ),))
    m.PartitionCellBySweepEdge(sweepPath=m.edges.findAt((xx1[i]+r, yy1[i], zz1[i]),),cells=m.cells.findAt((0.01, 0.01, 0.01),),edges=(m.edges.findAt((xx1[i], yy1[i]+a,zz1[i]+b), ),))	    
    
    
PP = 100    ### number of random points

dis=numpy.zeros(PP)

NO_INCL = 30

Q=51    # Number of iterations
W=12    # change of frequency
V=2    # change of radius or % volume

applyPBC = True
mixedBC = False

max_incl = 36     # set number of inclusions required

pow = -3.0 - math.log10(2*math.pi)


########################################################################
##                   Let's Create Materials 
########################################################################        
    
  




for w in range (1,W):
    
    for q in range (1,Q):
        
        w0=2*(math.pi)*10**pow ## CHANGING THE FREQUNCY WITH RESPECT TO pow
        mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d-%d' %(w,q))
        mdb.models['Model-%d-%d' %(w,q)].Material(name='Interphase')
        mdb.models['Model-%d-%d' %(w,q)].materials['Interphase'].Density(table=((1.0e-9, ), ))
        
        mdb.models['Model-%d-%d' %(w,q)].materials['Interphase'].Elastic(table=((
            0.9436270853, 0.49214), ))        
        mdb.models['Model-%d-%d' %(w,q)].materials['Interphase'].Viscoelastic(domain=TIME, table=((
            0.001497038998, 0.0, 0.032), (0.005959550149, 0.0, 0.1), (0.033512934, 0.0, 0.316), 
            (0.188457662, 0.0, 1.0), (0.5959550149, 0.0, 3.162), (0.11890852, 0.0, 10.0), (
            0.029868669, 0.0, 31.623), (0.01189109, 0.909091, 100.0), (0.007502477683, 
            0.030303, 316.228), (0.001497038998, 0.0, 1000.0), (0.001011931615, 0.0, 
            3162.278), (0.0007163379279, 0.0, 10000.0), (0.0006841563571, 0.0, 31622.777), (
            0.0006531666964, 0.0, 100000.0)), time=PRONY)
        
        mdb.models['Model-%d-%d' %(w,q)].Material(name='Inclusion')
        mdb.models['Model-%d-%d' %(w,q)].materials['Inclusion'].Density(table=((2.47*e-9, ), ))
        mdb.models['Model-%d-%d' %(w,q)].materials['Inclusion'].Elastic(table=
            ((64890, 0.249), ))
            
        mdb.models['Model-%d-%d' %(w,q)].Material(name='Matrix')
        mdb.models['Model-%d-%d' %(w,q)].materials['Matrix'].Density(table=((1.0e-9, ), ))
        
        mdb.models['Model-%d-%d' %(w,q)].materials['Matrix'].Elastic(table=((
            29.87551867, 0.493776), ))        
        mdb.models['Model-%d-%d' %(w,q)].materials['Matrix'].Viscoelastic(domain=TIME, table=((
            0.00017617, 0.0, 3.0), (0.00099076, 0.0,10.0), (0.0055714, 0.0, 32.0), 
            (0.01762, 0.0, 100.0), (0.05571, 0.0, 316.0), (0.3133, 0.0, 1000.0), (
            0.55714, 0.0, 3162.0), (0.03133, 0.8333333, 10000.0), (0.0078698, 
            0.0, 31623.0), (0.003133, 0.0, 100000.0), (0.00099076, 0.0, 
            316228.0), (0.00031328, 0.0, 1000000.0), (0.00017618, 0.0, 3162278.0), (
            0.00009906, 0.0, 10000000.0)), time=PRONY)          
        ## LETS CREATE MATRIX
        mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(name='__profile__', sheetSize=2.0)
        mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].sketchOptions.setValues(
            decimalPlaces=4)
        mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].rectangle(point1=(0.0, 0.0), 
            point2=(0.3, 0.3))
        mdb.models['Model-%d-%d' %(w,q)].Part(dimensionality=THREE_D, name='MATRIX', type=
            DEFORMABLE_BODY)
        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].BaseSolidExtrude(depth=0.3, sketch=
            mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
        del mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__']        
      


            
#####################################################################################################################
## DEFINE SECTIONS
#####################################################################################################################

        # MATRIX
        mdb.models['Model-%d-%d' %(w,q)].HomogeneousSolidSection(material='Interphase', name='MATRIX', 
            thickness=None)
        # INCLUSION
        mdb.models['Model-%d-%d' %(w,q)].HomogeneousSolidSection(material='Inclusion', name='INCLUSION', 
            thickness=None)  
        # INTERPHASE    
        #mdb.models['Model-%d-%d' %(w,q)].HomogeneousSolidSection(material='Matrix', name='INTERPHASE', 
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

        mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((0.01, 0.01, 0.01), ), )), sectionName='MATRIX', thicknessAssignment=FROM_SECTION)
        i=0
        while(i<num_incl):
            mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]-0.2*r, zz1[i]+0.2*r), ), )), sectionName='INCLUSION', thicknessAssignment=FROM_SECTION)
            mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].SectionAssignment(offset=0.0, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]+0.2*r, zz1[i]+0.2*r), ), )), sectionName='INCLUSION', thicknessAssignment=FROM_SECTION)
            #mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].SectionAssignment(offset=0.0, 
            #    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            #    cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]-1.05*r, zz1[i]), ), )), sectionName='INTERPHASE', thicknessAssignment=FROM_SECTION)
            #mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].SectionAssignment(offset=0.0, 
            #    offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            #    cells=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'].cells.findAt(((xx1[i], yy1[i]+1.05*r, zz1[i]), ), )), sectionName='INTERPHASE', thicknessAssignment=FROM_SECTION)                
            i=i+1
	
        ## Creating Instances
        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.DatumCsysByDefault(CARTESIAN)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Instance(dependent=OFF, name='MATRIX-1', 
            part=mdb.models['Model-%d-%d' %(w,q)].parts['MATRIX'])

 
    ########################################################################
    ##                   Let's Create Step 
    ########################################################################        

        mdb.models['Model-%d-%d' %(w,q)].ImplicitDynamicsStep(initialInc=(3*(10**((-pow)-3))), maxNumInc=100000, 
            name='Step-1', noStop=OFF, nohaf=OFF, previous='Initial', 
            timeIncrementationMethod=FIXED, timePeriod=(2*(10**(-pow))))
                
    ########################################################################
    ##                   Let's Create BC's and load
    ######################################################################## 
        mdb.models['Model-%d-%d' %(w,q)].PeriodicAmplitude(a_0=0.0, data=((0.0, 1.0), ), 
            frequency=w0, name='Amp-1', start=0.0, timeSpan=STEP)    
        if (mixedBC == True):
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-1', region=Region(
                faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.0, 0.015, 0.075), ), )), u1=0.0, u2=UNSET, u3=UNSET, ur1=UNSET, 
                ur2=UNSET, ur3=UNSET)
                
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-2', region=Region(
                faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.015, 0.0, 0.075), ), )), u1=UNSET, u2=0.0, u3=UNSET, ur1=UNSET, 
                ur2=UNSET, ur3=UNSET)				
            
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-3', region=Region(
                faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.075, 0.075, 0.0), ), )), u1=UNSET, u2=UNSET, u3=0.0, ur1=UNSET, 
                ur2=UNSET, ur3=UNSET)

            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude='Amp-1', createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-4', region=Region(
                faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.075, 0.075, 0.3), ), )), u1=UNSET, u2=UNSET, u3=0.01, ur1=UNSET, 
                ur2=UNSET, ur3=UNSET)

    ########################################################################
    ##                   Let's Mesh
    ########################################################################
        ## Creating Mesh on the instance
        
        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setElementType(elemTypes=(ElemType(
            elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
            elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
            regions=(mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((0.01,0.01,0.01), ), ), ))
            
        i=0	
        while(i<num_incl):
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setElementType(elemTypes=(ElemType(elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
                elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
                regions=(mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), ), ))
            #mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setElementType(elemTypes=(ElemType(elemCode=C3D8, elemLibrary=STANDARD), ElemType(elemCode=C3D6, 
            #    elemLibrary=STANDARD), ElemType(elemCode=C3D4, elemLibrary=STANDARD)), 
            #    regions=(mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((xx1[i],yy1[i]-1.05*r,zz1[i]), ), ((xx1[i],yy1[i]+1.05*r,zz1[i]), ), ), ))                
            i=i+1
 
        # Seed the instance	 
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedPartInstance(deviationFactor=0.1, 
            minSizeFactor=0.1, size=0.015, regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'], ))
        # Increase the seed on edges    
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.075, 0.0, 0.0), 
            )), minSizeFactor=0.1, size=0.03)    
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.0, 0.075, 0.0), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.0, 0.0, 0.075), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.075, 0.3, 0.0), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.075, 0.0, 0.3), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.3, 0.075, 0.0), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.0, 0.075, 0.3), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.3, 0.0, 0.075), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.0, 0.3, 0.075), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.075, 0.3, 0.3), 
            )), minSizeFactor=0.1, size=0.03)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.3, 0.075, 0.3), 
            )), minSizeFactor=0.1, size=0.03)            
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedEdgeBySize(deviationFactor=0.1, 
            edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(((0.3, 0.3, 0.075), 
            )), minSizeFactor=0.1, size=0.03)
            
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setMeshControls(elemShape=TET, regions=
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((0.01, 0.01, 0.01), ), ), technique=FREE)
        
        i=0	
        while(i<num_incl):
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setMeshControls(elemShape=TET, regions=
        	    mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((xx1[i],yy1[i]-0.2*r,zz1[i]+0.2*r), ), ((xx1[i],yy1[i]+0.2*r,zz1[i]+0.2*r), ), ), technique=FREE)
            #mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setMeshControls(elemShape=TET, regions=
        	#    mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].cells.findAt(((xx1[i],yy1[i]-1.05*r,zz1[i]), ), ((xx1[i],yy1[i]+1.05*r,zz1[i]), ), ), technique=FREE)                
            i=i+1
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.generateMesh(regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'], ))
   
        if (mixedBC == True):
            mdb.models['Model-%d-%d' %(w,q)].Tie(adjust=OFF, master=Region(
                edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(
                ((0.3, 0.15, 0.0), ), )), name='Constraint-1', positionToleranceMethod=
                COMPUTED, slave=Region(
                side1Faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.3, 0.15, 0.15), ), )), thickness=ON, tieRotations=OFF)
            mdb.models['Model-%d-%d' %(w,q)].Tie(adjust=OFF, master=Region(
                edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].edges.findAt(
                ((0.15, 0.3, 0.0), ), )), name='Constraint-2', positionToleranceMethod=
                COMPUTED, slave=Region(
                side1Faces=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(
                ((0.15, 0.3, 0.15), ), )), thickness=ON, tieRotations=ON)

        if (applyPBC == True):
            # create all nodesets for each face
            NegX=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.0, 0.15, 0.15), )), name='NegX')
            NegY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.15, 0.0, 0.15), )), name='NegY')
            NegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.15, 0.15, 0.0), )), name='NegZ')
            PosX=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.3, 0.15, 0.15), )), name='PosX')
            PosY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.15, 0.3, 0.15), )), name='PosY')
            PosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
                0.15, 0.15, 0.3), )), name='PosZ')                

            NegXNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=NegX.nodes, name='NegXNodes')
            NegYNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=NegY.nodes, name='NegYNodes')
            NegZNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=NegZ.nodes, name='NegZNodes')
            PosXNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=PosX.nodes, name='PosXNodes')
            PosYNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=PosY.nodes, name='PosYNodes')
            PosZNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=PosZ.nodes, name='PosZNodes')
            
            #Final vertex node sets to use in equation
            Vertex_NegXNegYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_NegXNegYNegZ", 
                operation=assembly.INTERSECTION, sets=[ NegXNodes, NegYNodes,NegZNodes])
            Vertex_PosXNegYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_PosXNegYNegZ", 
                operation=assembly.INTERSECTION, sets=[ PosXNodes, NegYNodes,NegZNodes])
            Vertex_NegXPosYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_NegXPosYNegZ", 
                operation=assembly.INTERSECTION, sets=[ NegXNodes, PosYNodes,NegZNodes])
            Vertex_NegXNegYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_NegXNegYPosZ", 
                operation=assembly.INTERSECTION, sets=[ NegXNodes, NegYNodes,PosZNodes])
            Vertex_PosXPosYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_PosXPosYPosZ", 
                operation=assembly.INTERSECTION, sets=[ PosXNodes, PosYNodes,PosZNodes])
            Vertex_NegXPosYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_NegXPosYPosZ", 
                operation=assembly.INTERSECTION, sets=[ NegXNodes, PosYNodes,PosZNodes])
            Vertex_PosXNegYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_PosXNegYPosZ", 
                operation=assembly.INTERSECTION, sets=[ PosXNodes, NegYNodes,PosZNodes])
            Vertex_PosXPosYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_PosXPosYNegZ", 
                operation=assembly.INTERSECTION, sets=[ PosXNodes, PosYNodes,NegZNodes])                

            #Final face node sets to use in equation
            NegX_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="NegX_red", operation=assembly.DIFFERENCE,
                sets=[ NegXNodes, NegYNodes,NegZNodes, PosYNodes, PosZNodes ])
            NegY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="NegY_red", operation=assembly.DIFFERENCE,
                sets=[ NegYNodes, NegZNodes, PosZNodes ])            
            NegZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="NegZ_red", operation=assembly.DIFFERENCE,
                sets=[ NegZNodes, Vertex_NegXNegYNegZ, Vertex_PosXNegYNegZ, Vertex_NegXPosYNegZ])
            PosX_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="PosX_red", operation=assembly.DIFFERENCE,
                sets=[ PosXNodes, NegYNodes,NegZNodes, PosYNodes, PosZNodes ])
            PosY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="PosY_red", operation=assembly.DIFFERENCE,
                sets=[ PosYNodes, NegZNodes, PosZNodes ])
            PosZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="PosZ_red", operation=assembly.DIFFERENCE,
                sets=[ PosZNodes, Vertex_NegXPosYPosZ, Vertex_PosXNegYPosZ, Vertex_NegXNegYPosZ])                
                
            #edge node sets     
            #Edge_NegXNegY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXNegY", operation=assembly.INTERSECTION,
            #    sets=[ NegXNodes, NegYNodes])
            #Edge_NegXNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXNegZ", operation=assembly.INTERSECTION,
            #    sets=[ NegXNodes, NegZNodes])
            #Edge_NegXPosY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXPosY", operation=assembly.INTERSECTION,
            #    sets=[ NegXNodes, PosYNodes])
            #Edge_NegXPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXPosZ", operation=assembly.INTERSECTION,
            #    sets=[ NegXNodes, PosZNodes])
            #Edge_PosXPosY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosY", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, PosYNodes])
            #Edge_PosXPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosZ", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, PosZNodes])
            #Edge_PosXPosY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosY", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, PosYNodes])
            #Edge_PosXPosY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosY", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, PosYNodes])                
            #Edge_PosXNegY=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXNegY", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, NegYNodes])
            #Edge_PosXNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXNegZ", operation=assembly.INTERSECTION,
            #    sets=[ PosXNodes, NegZNodes])
            #Edge_NegYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegYNegZ", operation=assembly.INTERSECTION,
            #    sets=[ NegYNodes, NegZNodes])
            #Edge_NegYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegYPosZ", operation=assembly.INTERSECTION,
            #    sets=[ NegYNodes, PosZNodes])                
            #Edge_PosYNegZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosYNegZ", operation=assembly.INTERSECTION,
            #    sets=[ PosYNodes, NegZNodes])
            #Edge_PosYPosZ=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosYPosZ", operation=assembly.INTERSECTION,
            #    sets=[ PosYNodes, PosZNodes])
                
            #Final edge node sets to use in equation
            #Edge_NegXNegY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXNegY_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegXNegY, Vertex_NegXNegYNegZ, Vertex_NegXNegYPosZ])
            #Edge_NegXNegZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXNegZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegXNegZ, Vertex_NegXNegYNegZ, Vertex_NegXPosYNegZ])
            #Edge_NegXPosY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXPosY_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegXPosY, Vertex_NegXPosYNegZ, Vertex_NegXPosYPosZ])
            #Edge_NegXPosZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegXPosZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegXPosZ, Vertex_NegXNegYPosZ, Vertex_NegXPosYPosZ])
            #    
            #Edge_PosXNegY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXNegY_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosXNegY, Vertex_PosXNegYNegZ, Vertex_PosXNegYPosZ])
            #Edge_PosXNegZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXNegZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosXNegZ, Vertex_PosXNegYNegZ, Vertex_PosXPosYNegZ])
            #Edge_PosXPosY_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosY_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosXPosY, Vertex_PosXPosYNegZ, Vertex_PosXPosYPosZ])
            #Edge_PosXPosZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosXPosZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosXPosZ, Vertex_PosXNegYPosZ, Vertex_PosXPosYPosZ]) 
            #    
            #Edge_NegYNegZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegYNegZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegYNegZ, Vertex_NegXNegYNegZ, Vertex_PosXNegYNegZ])
            #Edge_NegYPosZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_NegYPosZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_NegYPosZ, Vertex_NegXNegYPosZ, Vertex_PosXNegYPosZ])
            #Edge_PosYNegZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosYNegZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosYNegZ, Vertex_NegXPosYNegZ, Vertex_PosXPosYNegZ])
            #Edge_PosYPosZ_red=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Edge_PosYPosZ_red",
            #    operation=assembly.DIFFERENCE, sets=[ Edge_PosYPosZ, Vertex_NegXPosYPosZ, Vertex_PosXPosYPosZ])                
            
            
            PinTheNode(mdb.models['Model-%d-%d' %(w,q)],'Vertex_NegXNegYNegZ')

            # restrict 2 dofs of 3 vertex each
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-Z', 
                region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_NegXNegYPosZ'], u1=SET, u2=
                SET, u3=UNSET)    
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-Y', 
                region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_NegXPosYNegZ'], u1=SET, u2=
                UNSET, u3=SET)                 
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name=
                'BC-X', region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_PosXNegYNegZ'], u1=
                UNSET, u2=SET, u3=SET)
                
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude='Amp-1', createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'loading', region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_NegXNegYPosZ'], u1=UNSET, 
                u2=UNSET, u3=0.003) 

            # The nodesets and equations need to be handled as custom keywords because of limitations in the Python API (no support for *NSET, UNSORTED)
            # We will store them in a string
            KeywordsString = ''

            # Define name to be used for sorted version of each set.
            DataForSorting = [ [NegX_red    , "NegX_red_ordered"],
                               [NegY_red    , "NegY_red_ordered"],
                               [NegZ_red    , "NegZ_red_ordered"],
                               [PosX_red    , "PosX_red_ordered"],
                               [PosY_red    , "PosY_red_ordered"],
                               [PosZ_red    , "PosZ_red_ordered"],
                               [Vertex_NegXNegYNegZ , "Vertex_NegXNegYNegZ_ordered"],
                               [Vertex_PosXNegYNegZ , "Vertex_PosXNegYNegZ_ordered"],
                               [Vertex_NegXPosYNegZ , "Vertex_NegXPosYNegZ_ordered"],
                               [Vertex_NegXNegYPosZ , "Vertex_NegXNegYPosZ_ordered"],
                               [Vertex_PosXPosYPosZ , "Vertex_PosXPosYPosZ_ordered"],
                               [Vertex_NegXPosYPosZ , "Vertex_NegXPosYPosZ_ordered"],
                               [Vertex_PosXNegYPosZ , "Vertex_PosXNegYPosZ_ordered"],
                               [Vertex_PosXPosYNegZ , "Vertex_PosXPosYNegZ_ordered"]]
                               #[Edge_NegXNegY_red , "Edge_NegXNegY_red_ordered"],
                               #[Edge_NegXNegZ_red , "Edge_NegXNegZ_red_ordered"],
                               #[Edge_NegXPosY_red , "Edge_NegXPosY_red_ordered"],
                               #[Edge_NegXPosZ_red , "Edge_NegXPosZ_red_ordered"],
                               #[Edge_PosXNegY_red , "Edge_PosXNegY_red_ordered"],
                               #[Edge_PosXNegZ_red , "Edge_PosXNegZ_red_ordered"],
                               #[Edge_PosXPosY_red , "Edge_PosXPosY_red_ordered"],
                               #[Edge_PosXPosZ_red , "Edge_PosXPosZ_red_ordered"],
                               #[Edge_NegYNegZ_red , "Edge_NegYNegZ_red_ordered"],
                               #[Edge_NegYPosZ_red , "Edge_NegYPosZ_red_ordered"],
                               #[Edge_PosYNegZ_red , "Edge_PosYNegZ_red_ordered"],
                               #[Edge_PosYPosZ_red , "Edge_PosYPosZ_red_ordered"]]
                                          
            
            SortedNodes = {}
            for nodeset, name in DataForSorting:
                SortedNodes[name] = GenerateSortedMeshNodeArray(nodeset.nodes)
                KeywordsString += WriteMeshNodeArrayToUnsortedNodeSetDef(SortedNodes[name],name)+'\n'

            #Do the sorting and write the nodesets
            #equation 1----------------------------------------------------------------------
            NegX_allNodes = []
            sortedNodes_NegX = GenerateSortedMeshNodeArray(NegX_red.nodes)
            for nodeLable in sortedNodes_NegX:
                NegX_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))            

            PosX_allNodes = []
            sortedNodes_PosX = GenerateSortedMeshNodeArray(PosX_red.nodes)
            for nodeLable in sortedNodes_PosX:
                PosX_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))

            Pairs_allNodes_1=[]
            for ii, jj in zip(PosX_allNodes,NegX_allNodes):
                Pairs_allNodes_1.append([ii, jj, "Vertex_PosXNegYNegZ"])

            KeywordsString += "*EQUATION\n"
            for (masterName, slaveName, RP_vertex) in Pairs_allNodes_1:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n'                
            #equation 2------------------------------------------------------------------------
            NegY_allNodes = []
            sortedNodes_NegY = GenerateSortedMeshNodeArray(NegY_red.nodes)
            for nodeLable in sortedNodes_NegY:
                NegY_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))            

            PosY_allNodes = []
            sortedNodes_PosY = GenerateSortedMeshNodeArray(PosY_red.nodes)
            for nodeLable in sortedNodes_PosY:
                PosY_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))

            Pairs_allNodes_2=[]
            for ii, jj in zip(PosY_allNodes,NegY_allNodes):
                Pairs_allNodes_2.append([ii, jj, "Vertex_NegXPosYNegZ"])


            for (masterName, slaveName, RP_vertex) in Pairs_allNodes_2:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n' 
                                                        
            #equation 3 ------------------------------------------------------------------------                                            
            NegZ_allNodes = []
            sortedNodes_NegZ = GenerateSortedMeshNodeArray(NegZ_red.nodes)
            for nodeLable in sortedNodes_NegZ:
                NegZ_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))            

            PosZ_allNodes = []
            sortedNodes_PosZ = GenerateSortedMeshNodeArray(PosZ_red.nodes)
            for nodeLable in sortedNodes_PosZ:
                PosZ_allNodes.append('MATRIX-1'+'.'+str(nodeLable.label))

            Pairs_allNodes_3=[]
            for ii, jj in zip(PosZ_allNodes,NegZ_allNodes):
                Pairs_allNodes_3.append([ii, jj, "Vertex_NegXNegYPosZ"])


            for (masterName, slaveName, RP_vertex) in Pairs_allNodes_3:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n'                                                         
                                                        
                                                        
            #print KeywordsString              
    ########################################################################
    ##                   Let's Create Sets
    ########################################################################	
	
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(faces=
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['MATRIX-1'].faces.findAt(((
            0.075, 0.075, 0.0), )), name='Bottom_Face')
        	
        mdb.models['Model-%d-%d' %(w,q)].HistoryOutputRequest(createStepName='Step-1', frequency=1
            , name='H-Output-2', rebar=EXCLUDE, region=
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Bottom_Face'], sectionPoints=DEFAULT, 
            variables=('S33', 'E33'))

    ########################################################################
    ##                   Let's Create Jobs
    ########################################################################            

        mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
            explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
            memory=90, memoryUnits=PERCENTAGE, model='Model-%d-%d' %(w,q), modelPrint=OFF, 
            multiprocessingMode=DEFAULT, name='JobPBC-%d-%d' %(w,q) , nodalOutputPrecision=SINGLE, 
            numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
            waitHours=0, waitMinutes=0)
        #mdb.jobs['JobPBC-%d-%d' %(w,q)].writeInput()
        #mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
        #mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()
        AddEquationsAndUnsortedNsetStringToKeywords(mdb.models['Model-%d-%d' %(w,q)],KeywordsString)
        mdb.jobs['JobPBC-%d-%d' %(w,q)].writeInput()        
    pow = pow + 0.2    ## Next frequency

# Save by skulka17 on Mon Apr 18 10:48:15 2016
