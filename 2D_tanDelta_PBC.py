##############################################
##############################################
## 2D TAN DELTA VARIATION  _ PBC                ##
## SHANK KULKARNI                           ##
## 10/03/2016                               ##
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
    
    for i in range(2):
        EqTerms = [ (-1.0,SlaveNodesetName,i+1), 
                    ( 1.0,MasterNodesetName,i+1) ]

        if (RP_vertex == "Vertex_3"):
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
        
        
PP = 1000   ### number of random points
dis=numpy.zeros(PP)

NO_INCL = 30

applyPBC = True
mixedBC = False

rad=0.2
Q=51    # Number of iterations
W=12    # change of frequency
V=2    # change of radius or % volume


max_incl = 13      # set number of inclusions required

pow = -3.0 - math.log10(2*math.pi)





for w in range (1,W):
  
    for q in range (1,Q):
        
        w0=2*(math.pi)*10**pow ## CHANGING THE FREQUNCY WITH RESPECT TO pow
        mdb.Model(modelType=STANDARD_EXPLICIT, name='Model-%d-%d' %(w,q))
        mdb.models['Model-%d-%d' %(w,q)].Material(name='VE')
        mdb.models['Model-%d-%d' %(w,q)].materials['VE'].Density(table=((1.0e-9, ), ))
        
        mdb.models['Model-%d-%d' %(w,q)].materials['VE'].Elastic(table=((
            0.9436270853, 0.49214), ))        
        mdb.models['Model-%d-%d' %(w,q)].materials['VE'].Viscoelastic(domain=TIME, table=((
            0.001497038998, 0.0, 0.032), (0.005959550149, 0.0, 0.1), (0.033512934, 0.0, 0.316), 
            (0.188457662, 0.0, 1.0), (0.5959550149, 0.0, 3.162), (0.11890852, 0.0, 10.0), (
            0.029868669, 0.0, 31.623), (0.01189109, 0.909091, 100.0), (0.007502477683, 
            0.030303, 316.228), (0.001497038998, 0.0, 1000.0), (0.001011931615, 0.0, 
            3162.278), (0.0007163379279, 0.0, 10000.0), (0.0006841563571, 0.0, 31622.777), (
            0.0006531666964, 0.0, 100000.0)), time=PRONY)
        
        mdb.models['Model-%d-%d' %(w,q)].Material(name='Elastic')
        mdb.models['Model-%d-%d' %(w,q)].materials['Elastic'].Density(table=((2.47*e-9, ), ))
        mdb.models['Model-%d-%d' %(w,q)].materials['Elastic'].Elastic(table=
            ((64890, 0.249), ))
##---------------------------------------------------------------------------------------------------        
        mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(name='__profile__', sheetSize=2.0)
        mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].rectangle(point1=(-2.0, 2.0), 
            point2=(2.0, -2.0))
        mdb.models['Model-%d-%d' %(w,q)].Part(dimensionality=TWO_D_PLANAR, name='Part-1', type=
            DEFORMABLE_BODY)
        mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].BaseShell(sketch=
            mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
        del mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__']
        mdb.models['Model-%d-%d' %(w,q)].ConstrainedSketch(gridSpacing=1.8, name='__profile__', 
            sheetSize=2, transform=
            mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].MakeSketchTransform(
            sketchPlane=mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].faces[0], 
            sketchPlaneSide=SIDE1, sketchOrientation=RIGHT, origin=(0.0, 0.0, 0.0)))
        mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].projectReferencesOntoSketch(filter=
            COPLANAR_EDGES, sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
    
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

            mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'].CircleByCenterPerimeter(center=(
                xx1[i], yy1[i]), point1=((xx1[i]-rad), yy1[i]))
    
            mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].PartitionFaceBySketch(faces=
                mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].faces.getSequenceFromMask(('[#1 ]', 
                ), ), sketch=mdb.models['Model-%d-%d' %(w,q)].sketches['__profile__'])
                        

    ########################################################################
    ##                   Let's Create Sections 
    ########################################################################

                
            
        mdb.models['Model-%d-%d' %(w,q)].HomogeneousSolidSection(material='VE', name='Matrix', 
            thickness=None)
        mdb.models['Model-%d-%d' %(w,q)].HomogeneousSolidSection(material='Elastic', name='Inclusion', 
            thickness=None)
    

    
        mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].SectionAssignment(offset=0.0, 
            offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
            faces=mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].faces.findAt(((1.9, 
            1.9, 0.0), (0.0, 0.0, 1.0)), )), sectionName='Matrix', 
            thicknessAssignment=FROM_SECTION)

    
        for i in range (num_incl):    

            mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].SectionAssignment(offset=0.2, 
                offsetField='', offsetType=MIDDLE_SURFACE, region=Region(
                faces=mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'].faces.findAt((((xx1[i]), 
                (yy1[i]), 0.0), (0.0, 0.0, 1.0)), )), sectionName='Inclusion', 
                thicknessAssignment=FROM_SECTION)    

        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.DatumCsysByDefault(CARTESIAN)
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Instance(dependent=OFF, name='Part-1-1', 
            part=mdb.models['Model-%d-%d' %(w,q)].parts['Part-1'])
        
        
        

        mdb.models['Model-%d-%d' %(w,q)].ImplicitDynamicsStep(initialInc=(3*(10**((-pow)-3))), maxNumInc=100000, 
            name='Step-1', noStop=OFF, nohaf=OFF, previous='Initial', 
            timeIncrementationMethod=FIXED, timePeriod=(3*(10**(-pow))))

        mdb.models['Model-%d-%d' %(w,q)].PeriodicAmplitude(a_0=0.0, data=((0.0, 1.0), ), 
            frequency=w0, name='Amp-1', start=0.0, timeSpan=STEP)
                
     
        if (mixedBC == True):
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-1', region=Region(
                edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(
                ((-5.0, -10.0, 0.0), ), )), u1=UNSET, u2=0.0, ur3=UNSET)
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-2', region=Region(
                edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(
                ((-10.0, 5.0, 0.0), ), )), u1=0.0, u2=UNSET, ur3=UNSET)

            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude='Amp-1', createStepName='Step-1'
                , distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'BC-3', region=Region(
                edges=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(
                ((5.0, 10.0, 0.0), ), )), u1=UNSET, u2=0.1, ur3=UNSET)
        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.seedPartInstance(deviationFactor=0.1, 
            minSizeFactor=0.1, regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'], ), size=0.1)
        
        for i in range (num_incl):
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setElementType(elemTypes=(ElemType(
                elemCode=CPE8RH, elemLibrary=STANDARD), ElemType(elemCode=CPE6M, 
                elemLibrary=STANDARD)),regions=(
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].faces.findAt(((
                (xx1[i]), (yy1[i]), 0.0), )), ))        
        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.setElementType(elemTypes=(ElemType(
            elemCode=CPE8RH, elemLibrary=STANDARD), ElemType(elemCode=CPE6M, 
            elemLibrary=STANDARD)),  regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].faces.findAt(((
            1.9, 1.9, 0.0), )), ))
        
        
        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.generateMesh(regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'], ))    

        mdb.models['Model-%d-%d' %(w,q)].rootAssembly.generateMesh(regions=(
            mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'], ))
       
        if (applyPBC == True):
            BottomEdge=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(edges=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
                -1.0, -2.0, 0.0), )), name='BottomEdge')        
            TopEdge=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(edges=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
                -1.0, 2.0, 0.0), )), name='TopEdge')
            LeftEdge=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(edges=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
                -2.0, -1.0, 0.0), )), name='LeftEdge')
            RightEdge=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(edges=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.instances['Part-1-1'].edges.findAt(((
                2.0, -1.0, 0.0), )), name='RightEdge')            
            BottomEdgeNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=BottomEdge.nodes, name='BottomEdgeNodes')
            TopEdgeNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=TopEdge.nodes, name='TopEdgeNodes')    
            LeftEdgeNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=LeftEdge.nodes, name='LeftEdgeNodes')   
            RightEdgeNodes=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.Set(nodes=RightEdge.nodes, name='RightEdgeNodes')
            
            #creating 4 vertex
            Vertex_1 = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_1", 
                operation=assembly.INTERSECTION, sets=[ LeftEdgeNodes, BottomEdgeNodes ])
            Vertex_2 = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_2", 
                operation=assembly.INTERSECTION, sets=[ RightEdgeNodes, BottomEdgeNodes ])
            Vertex_3 = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_3", 
                operation=assembly.INTERSECTION, sets=[ LeftEdgeNodes, TopEdgeNodes ])
            Vertex_4 = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="Vertex_4", 
                operation=assembly.INTERSECTION, sets=[ RightEdgeNodes, TopEdgeNodes ])
                
            PinTheNode(mdb.models['Model-%d-%d' %(w,q)],'Vertex_1')  

            # restrict 1 dof of 2 vertex each
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Initial', 
                distributionType=UNIFORM, fieldName='', localCsys=None, name='BC-4', 
                region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_3'], u1=SET, u2=
                UNSET, ur3=UNSET)            
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude=UNSET, createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'V_2', region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_2'], u1=
                UNSET, u2=0.0, ur3=UNSET)
            mdb.models['Model-%d-%d' %(w,q)].DisplacementBC(amplitude='Amp-1', createStepName='Step-1', 
                distributionType=UNIFORM, fieldName='', fixed=OFF, localCsys=None, name=
                'loading', region=mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['Vertex_3'], u1=UNSET, 
                u2=0.01, ur3=UNSET)            
                
            #APPLICATION OF EQUATION CONSTRAINT
            # first reduced edge node sets (without vertex)
            TopEdge_red = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="TopEdge_red", 
                operation=assembly.DIFFERENCE, sets=[ TopEdgeNodes, Vertex_3, Vertex_4 ])            
            BottomEdge_red = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="BottomEdge_red", 
                operation=assembly.DIFFERENCE, sets=[ BottomEdgeNodes, Vertex_1, Vertex_2 ])
            LeftEdge_red = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="LeftEdge_red", 
                operation=assembly.DIFFERENCE, sets=[ LeftEdgeNodes, Vertex_1, Vertex_3 ])
            RightEdge_red = mdb.models['Model-%d-%d' %(w,q)].rootAssembly.SetByBoolean(name="RightEdge_red", 
                operation=assembly.DIFFERENCE, sets=[ RightEdgeNodes, Vertex_2, Vertex_4 ])                
            
            # The nodesets and equations need to be handled as custom keywords because of limitations in the Python API (no support for *NSET, UNSORTED)
            # We will store them in a string
            KeywordsString = ''

            # Define name to be used for sorted version of each set.
            DataForSorting = [ [TopEdge_red    , "TopEdge_red_ordered"],
                               [BottomEdge_red , "BottomEdge_red_ordered"],
                               [LeftEdge_red   , "LeftEdge_red_ordered"],
                               [RightEdge_red  , "RightEdge_red_ordered"],
                               [Vertex_3, "Vertex_3_ordered"],
                               [Vertex_4, "Vertex_4_ordered"]]            

            #Do the sorting and write the nodesets
            Top_allNodes = []
            sortedNodes_Top = GenerateSortedMeshNodeArray(TopEdge_red.nodes)
            for nodeLable in sortedNodes_Top:
                Top_allNodes.append('Part-1-1'+'.'+str(nodeLable.label))            

            Bottom_allNodes = []
            sortedNodes_Bottom = GenerateSortedMeshNodeArray(BottomEdge_red.nodes)
            for nodeLable in sortedNodes_Bottom:
                Bottom_allNodes.append('Part-1-1'+'.'+str(nodeLable.label))

            Pairs_allNodes_1=[]
            for ii, jj in zip(Bottom_allNodes,Top_allNodes):
                Pairs_allNodes_1.append([ii, jj, "Vertex_3"])
 
            Left_allNodes = []
            sortedNodes_Left = GenerateSortedMeshNodeArray(LeftEdge_red.nodes)
            for nodeLable in sortedNodes_Left:
                Left_allNodes.append('Part-1-1'+'.'+str(nodeLable.label))            

            Right_allNodes = []
            sortedNodes_Right = GenerateSortedMeshNodeArray(RightEdge_red.nodes)
            for nodeLable in sortedNodes_Right:
                Right_allNodes.append('Part-1-1'+'.'+str(nodeLable.label))

            Pairs_allNodes_2=[]
            for ii, jj in zip(Right_allNodes,Left_allNodes):
                Pairs_allNodes_1.append([ii, jj, "Vertex_2"])


            
            keyWordEquation = ''
            SortedNodes = {}
            for nodeset, name in DataForSorting:
                SortedNodes[name] = GenerateSortedMeshNodeArray(nodeset.nodes)
                KeywordsString += WriteMeshNodeArrayToUnsortedNodeSetDef(SortedNodes[name],name)+'\n'

                
            # These are the pairs of faces/edges/vertices that need to be constrained.
            Pairs = [["Vertex_3_ordered", "Vertex_4_ordered", "Vertex_2"]]
            # Now loop through and define the equations for each of these nodes.
            KeywordsString += "*EQUATION\n"
            for (masterName, slaveName, RP_vertex) in Pairs_allNodes_1:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n'

            for (masterName, slaveName, RP_vertex) in Pairs_allNodes_2:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n'

            for (masterName, slaveName, RP_vertex) in Pairs:
                KeywordsString += DefinePeriodicEquation3D(masterName,
                                                        slaveName,
                                                        RP_vertex) + '\n'                                                        
                               
            #HISTORY OUTPUT REQUESTS FOR ALL SETS
            mdb.models['Model-%d-%d' %(w,q)].HistoryOutputRequest(createStepName='Step-1', frequency=1
                , name='H-Output-5', rebar=EXCLUDE, region=
                mdb.models['Model-%d-%d' %(w,q)].rootAssembly.sets['BottomEdge'], sectionPoints=DEFAULT, 
                variables=('S22', 'E22'))
                
    
            mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
                explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
                memory=90, memoryUnits=PERCENTAGE, model='Model-%d-%d' %(w,q), modelPrint=OFF, 
                multiprocessingMode=DEFAULT, name='JobPBC-%d-%d' %(w,q) , nodalOutputPrecision=SINGLE, 
                numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
                waitHours=0, waitMinutes=0)
            
            #mdb.jobs['Job-%d-%d' %(w,q) ].submit(consistencyChecking=OFF)    
            #mdb.jobs['Job-%d-%d' %(w,q) ].waitForCompletion()
            AddEquationsAndUnsortedNsetStringToKeywords(mdb.models['Model-%d-%d' %(w,q)],KeywordsString)
            mdb.jobs['JobPBC-%d-%d' %(w,q)].writeInput()
    pow = pow + 0.2    ## Next frequency

# Save by skulka17 on Mon Apr 18 10:48:15 2016
