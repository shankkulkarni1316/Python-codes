##########################################################
##########################################################
## CONTAINS 3 MAIN FUNCTIONS: (FOR THERMAL ANALYSIS)    ##
##  1) APPLICATION OF PBC FOR CONFORMAL MESHING         ##
##  2) APPLICATION ON PBC FOR NON-CONFORMAL MESHING     ##
##  3) UniformSurfaceGradient BC                                        ##
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
#Abaqus API imports
import interaction
import mesh
import assembly
#There is an executive for this script at the bottom of the script.

from ConstraintsUtilities import CreateInitialSets
from ConstraintsUtilities import MapSubsetFromCopiedNodes
from ConstraintsUtilities import GetSurfaceFromNodeSet
from ConstraintsUtilities import AddEquationsAndUnsortedNsetStringToKeywords
from ConstraintsUtilities import WriteEquationTerms
from ConstraintsUtilities import GetInChunks
from ConstraintsUtilities import WriteMeshNodeArrayToUnsortedNodeSetDef
from ConstraintsUtilities import GenerateSortedMeshNodeArray
from ConstraintsUtilities import GetPeriodicityVector
from ConstraintsUtilities import ThermalRP
from ConstraintsUtilities import TheCleanUpFuction


def applyThermalPBC(Model):
    """
    applyThermalPBC:
    this function checks if meshing is conformal or not and then calls appropriate function to apply PBC for thermal analysis.
        
    Inputs:
    Model -                 A model object from a Model Database
    """
 
    CreateInitialSets(Model)
    #Define references to the sets that contain all the nodes on each face
    RVE_NegX = Model.rootAssembly.sets['RVE_NegX']
    RVE_PosX = Model.rootAssembly.sets['RVE_PosX']
    RVE_NegY = Model.rootAssembly.sets['RVE_NegY']
    RVE_PosY = Model.rootAssembly.sets['RVE_PosY']
    RVE_NegZ = Model.rootAssembly.sets['RVE_NegZ']
    RVE_PosZ = Model.rootAssembly.sets['RVE_PosZ']
    
    DataForSortingForCheck = [ [RVE_NegX , "RVE_NegX_ordered"],
                             [RVE_PosX , "RVE_PosX_ordered"],
                             [RVE_NegY , "RVE_NegY_ordered"],
                             [RVE_PosY , "RVE_PosY_ordered"],
                             [RVE_NegZ , "RVE_NegZ_ordered"],
                             [RVE_PosZ , "RVE_PosZ_ordered"] ]    
    

    SortedNodes = {}
    for nodeset, name in DataForSortingForCheck:
        SortedNodes[name] = GenerateSortedMeshNodeArray(nodeset.nodes)
        #KeywordsString += WriteMeshNodeArrayToUnsortedNodeSetDef(SortedNodes[name],name)+'\n'

        

    Pairs = [ [ "RVE_NegX_ordered", "RVE_PosX_ordered" ],                    
              [ "RVE_NegY_ordered", "RVE_PosY_ordered" ],                    
              [ "RVE_NegZ_ordered", "RVE_PosZ_ordered" ] ]

    isMeshConformal = False # to check conformality
    
    for (masterName, slaveName) in Pairs:
        PeriodicityVec = GetPeriodicityVector(SortedNodes[masterName],SortedNodes[slaveName])
        if PeriodicityVec == None:
            isMeshConformal = True

    if (isMeshConformal == False):
        print "Mesh is Conformal."
        return GenerateThermalPBCs(Model)
    if (isMeshConformal == True): 
        print "Mesh is Non Conformal."
        return GenerateThermalPBCsNC(Model)    

    
    
def GenerateThermalPBCs( Model):
    """
    Define thermal periodicity on a CAE model object.
    
    This method generates equations that will enforce periodicity
    on a 3D hexahedral RVE that has a periodic mesh - that is, the
    mesh topology on one face exactly matches the mesh topology on
    the opposite face.  Creates reference node for specifying 
    far-field temperature gradients
    
    Inputs:
    Model -                 A model object from a Model Database

    
    Outputs:
    Returns the string containing the keywords to be written to the input deck.
    
    The Model Object will be modified to include the appropriate constraints.
    """

    #Define references to the sets that contain all the nodes on each face
    RVE_NegX = Model.rootAssembly.sets['RVE_NegX']
    RVE_PosX = Model.rootAssembly.sets['RVE_PosX']
    RVE_NegY = Model.rootAssembly.sets['RVE_NegY']
    RVE_PosY = Model.rootAssembly.sets['RVE_PosY']
    RVE_NegZ = Model.rootAssembly.sets['RVE_NegZ']
    RVE_PosZ = Model.rootAssembly.sets['RVE_PosZ']
        
    #Trim the nodesets to eliminate multiple slaving of nodes.
    RVE_NegY_red = Model.rootAssembly.SetByBoolean(name="RVE_NegY_red", operation=assembly.DIFFERENCE, sets=[ RVE_NegY, RVE_PosX ])
    RVE_PosY_red = Model.rootAssembly.SetByBoolean(name="RVE_PosY_red", operation=assembly.DIFFERENCE, sets=[ RVE_PosY, RVE_PosX ])
    RVE_NegZ_red = Model.rootAssembly.SetByBoolean(name="RVE_NegZ_red", operation=assembly.DIFFERENCE, sets=[ RVE_NegZ, RVE_PosX, RVE_PosY ])
    RVE_PosZ_red = Model.rootAssembly.SetByBoolean(name="RVE_PosZ_red", operation=assembly.DIFFERENCE, sets=[ RVE_PosZ, RVE_PosX, RVE_PosY ])
        
    # After these sets are built, we need to sort each one based on position so that paired sets have matching order.
    # We can't set the sort order in a node set via the Python API, so we're going to make MeshNodeArray objects
    # with the proper sorting and then write the input file code directly.  An alternaitve to this would be to
    # create an individual equation for each pair of nodes, but since the API only supports putting sets into
    # equations, not only would this require lots of equation definitions, it would also require a node set for
    # every single node on the surface of the RVE, so we aren't going to do it that way.
    
    # The nodesets and equations need to be handled as custom keywords because of limitations in the Python API (no support for *NSET, UNSORTED)
    # We will store them in a string
    KeywordsString = ''
    
    # Define name to be used for sorted version of each set.
    DataForSorting = [ [RVE_NegX     , "RVE_NegX_ordered"],
                       [RVE_PosX     , "RVE_PosX_ordered"],
                       [RVE_NegY_red , "RVE_NegY_ordered"],
                       [RVE_PosY_red , "RVE_PosY_ordered"],
                       [RVE_NegZ_red , "RVE_NegZ_ordered"],
                       [RVE_PosZ_red , "RVE_PosZ_ordered"] ]
    
    #Do the sorting and write the nodesets
    SortedNodes = {}
    for nodeset, name in DataForSorting:
        SortedNodes[name] = GenerateSortedMeshNodeArray(nodeset.nodes)
        KeywordsString += WriteMeshNodeArrayToUnsortedNodeSetDef(SortedNodes[name],name)+'\n'
        
    #Create periodic reference nodes and nodesets pointing to the nodes.
    #Each node's DOF 11 will correspond to a different far-field temperature gradient.
    #We have to do it this way because Abaqus doesn't admit using DOFs reserved for
    #mechanical quantities for thermal
    ThermalRP(Model)

        
    # These are the pairs of faces/edges/vertices that need to be constrained.
    Pairs = [ [ "RVE_NegX_ordered", "RVE_PosX_ordered" ],                    
              [ "RVE_NegY_ordered", "RVE_PosY_ordered" ],                    
              [ "RVE_NegZ_ordered", "RVE_PosZ_ordered" ] ]
    # Now loop through and define the equations for each of these nodes.
    KeywordsString += "*EQUATION\n"
    for (masterName, slaveName) in Pairs:
        PeriodicityVec = GetPeriodicityVector(SortedNodes[masterName],SortedNodes[slaveName])
        if PeriodicityVec == None:
            print "ERROR - Could not match nodes on " + masterName + " and " + slaveName
            print "It is possible that the node positions on the faces do not match within the required tolerance"
            raise RuntimeError("Unmatched Faces")
        KeywordsString += DefineScalarFieldPeriodicEquation3D(masterName,
                                                              slaveName,
                                                              11,
                                                              [ ("RP-ThermGrad_x",11), ("RP-ThermGrad_y",11), ("RP-ThermGrad_z",11) ],
                                                              PeriodicityVec) + '\n'
    return KeywordsString                                                          
    #AddEquationsAndUnsortedNsetStringToKeywords(Model,KeywordsString)
    

 
def GenerateThermalPBCsNC( Model):
    """
    Define thermal periodicity on a CAE model object.
    
    This method generates equations that will enforce periodicity
    on a 3D hexahedral RVE that has a non periodic mesh - that is, the
    mesh topology on one face does not matche the mesh topology on
    the opposite face.  Creates reference node for specifying 
    far-field temperature gradients
    
    Inputs:
    Model -                 A model object from a Model Database

    
    Outputs:
    Returns the string containing the keywords to be written to the input deck.
    
    The Model Object will be modified to include the appropriate constraints.
    """

    #Define references to the sets that contain all the nodes on each face
    RVE_NegX = Model.rootAssembly.sets['RVE_NegX']
    RVE_PosX = Model.rootAssembly.sets['RVE_PosX']
    RVE_NegY = Model.rootAssembly.sets['RVE_NegY']
    RVE_PosY = Model.rootAssembly.sets['RVE_PosY']
    RVE_NegZ = Model.rootAssembly.sets['RVE_NegZ']
    RVE_PosZ = Model.rootAssembly.sets['RVE_PosZ'] 
    
    #Trim the nodesets to eliminate multiple slaving of nodes.
    RVE_NegY_red = Model.rootAssembly.SetByBoolean(name="RVE_NegY_red", operation=assembly.DIFFERENCE, sets=[ RVE_NegY, RVE_NegX ])
    RVE_NegZ_red = Model.rootAssembly.SetByBoolean(name="RVE_NegZ_red", operation=assembly.DIFFERENCE, sets=[ RVE_NegZ, RVE_NegX, RVE_NegY ])
	
	####---------------------------------------------------------------------------------------------------------------------------------
	## Creating all FACE, EDGE and VERTEX nodes here for original mesh (befor copying the mesh topology)
	## GUIDELINES (Face node sets don't include edge or vertex nodes and edge node sets don't include vertex nodes)    
    ####--------------------------------------------------------------------------------------------------------------------------------

    Edge_AlignedPosXNegY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXNegY", operation=assembly.INTERSECTION, sets=[ RVE_NegY, RVE_PosX ])
    Edge_AlignedPosXPosY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXPosY", operation=assembly.INTERSECTION, sets=[ RVE_PosY, RVE_PosX ])
    Edge_AlignedPosXNegZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXNegZ", operation=assembly.INTERSECTION, sets=[ RVE_NegZ, RVE_PosX ])
    Edge_AlignedPosXPosZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXPosZ", operation=assembly.INTERSECTION, sets=[ RVE_PosZ, RVE_PosX ])
    Edge_AlignedNegXNegY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXNegY", operation=assembly.INTERSECTION, sets=[ RVE_NegY, RVE_NegX ])
    Edge_AlignedNegXPosY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXPosY", operation=assembly.INTERSECTION, sets=[ RVE_PosY, RVE_NegX ])
    Edge_AlignedNegXNegZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXNegZ", operation=assembly.INTERSECTION, sets=[ RVE_NegZ, RVE_NegX ])
    Edge_AlignedNegXPosZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXPosZ", operation=assembly.INTERSECTION, sets=[ RVE_PosZ, RVE_NegX ])

    Vertex_PosXPosYPosZ = Model.rootAssembly.SetByBoolean(name="Vertex_PosXPosYPosZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedPosXPosY, Edge_AlignedPosXPosZ ])
    Vertex_PosXNegYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_PosXNegYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedPosXNegY, Edge_AlignedPosXNegZ ])
    Vertex_NegXPosYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXPosYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXPosY, Edge_AlignedNegXNegZ ])
    Vertex_NegXNegYPosZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXNegYPosZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXNegY, Edge_AlignedNegXPosZ ])
    Vertex_NegXNegYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXNegYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXNegY, Edge_AlignedNegXNegZ ])	    
	#-------------------------------------------------------------------------------------------------------------------------------------
	# Lets determine the height, width and length of RVE so that we can shift orphan mesh accordingly
	# for that first find out co ordinates of 4 vertices and calculate dimensions accordingly	
	#------------------------------------------------------------------------------------------------------------------------------------
    V0=	Model.rootAssembly.sets['Vertex_NegXNegYNegZ'].nodes[0].coordinates
    V1= Model.rootAssembly.sets['Vertex_PosXNegYNegZ'].nodes[0].coordinates
    V2=	Model.rootAssembly.sets['Vertex_NegXPosYNegZ'].nodes[0].coordinates
    V3= Model.rootAssembly.sets['Vertex_NegXNegYPosZ'].nodes[0].coordinates
	
    Xlength = abs(V0[0]-V1[0])          # length of RVE in X direction
    Ylength = abs(V0[1]-V2[1])			# length of RVE in Y direction
    Zlength = abs(V0[2]-V3[2])			# length of RVE in Z direction    
    ## ###################################################################################################################################
    ## This part of code is written to replace the code from 3D PLM people 
    ## this code is more robust and can handle all type of surfaces (surfaces which are not flat as well)
    ## 
    ## ### FUCTION NAMED "GetSurfaceFromNodeSet"  IS WRITTEN AT THE END OF THE CODE WHICH RETURNS SURFACE WHEN NODE SET IS GIVEN
    ## 
    ## Here using that function surface is created at NegX, NegY and NegZ
    ## then offset shell mesh is created so that it will coincident with opposite surface
    ## Also surface is created at PosX, PosY and PosZ nodes and then both surfaces are tied together with PosX as master and copied as slave
    ## 
    ## ####################################################################################################################################
    Assembly=Model.rootAssembly
    PartNames = Model.parts.keys()
    if len(PartNames)>1:
        print "GeneratePBCsNC is only set up to handle a single part"
        raise NotImplementedError("GeneratePBCsNC is only set up to handle a single part")
    Part = Model.parts[PartNames[0]]
    pp=Part.elements
    Part.Set(elements=pp, name='all_solid_elements')
    InstanceNames = Model.rootAssembly.instances.keys()
    InstanceName = InstanceNames[0]
    if len(InstanceNames)>1:
        print "GeneratePBCs is only set up to handle a single instance"
        raise NotImplementedError("GeneratePBCs is only set up to handle a single instance")
    Instance=Model.rootAssembly.instances[InstanceNames[0]]    
 	#-------------------------------------------------------------------------------------------------------------------------------------
    # This part creats orphan mesh similar to mesh on NegX, NegY and NegZ surface and offset it to PosX, PosY and PosZ surface respectively	
    #--------------------------------------------------------------------------------------------------------------------------------------
    ## New idea of using set 
    Surf1 = GetSurfaceFromNodeSet(Model,'RVE_NegX')

    includedSurface = Model.rootAssembly.Surface( name='NegXoffset', face1Elements=Surf1)
    
    Assembly.generateMeshByOffset(region=Assembly.surfaces['NegXoffset'], meshType=SHELL,distanceBetweenLayers=0.0, numLayers=1, initialOffset=-Xlength)	
    elementsX=Part.elements
    Part.Set(elements=elementsX, name='all_solid_elementsX')    
    Part.SetByBoolean(name='ShellX', operation=DIFFERENCE, sets=(
    Part.sets['all_solid_elementsX'], Part.sets['all_solid_elements'], ))


    Surf2 = GetSurfaceFromNodeSet(Model,'RVE_NegY')
    includedSurface = Model.rootAssembly.Surface( name='NegYoffset', face1Elements=Surf2)    
    Assembly.generateMeshByOffset(region=Assembly.surfaces['NegYoffset'], meshType=SHELL,distanceBetweenLayers=0.0, numLayers=1, initialOffset=-Ylength)	
    elementsY=Part.elements
    Part.Set(elements=elementsY, name='all_solid_elementsY')    
    Part.SetByBoolean(name='ShellY', operation=DIFFERENCE, sets=(
    Part.sets['all_solid_elementsY'], Part.sets['all_solid_elements'], Part.sets['all_solid_elementsX'], ))


    Surf3 = GetSurfaceFromNodeSet(Model,'RVE_NegZ')
    includedSurface = Model.rootAssembly.Surface( name='NegZoffset', face1Elements=Surf3)    
    Assembly.generateMeshByOffset(region=Assembly.surfaces['NegZoffset'], meshType=SHELL,distanceBetweenLayers=0.0, numLayers=1, initialOffset=-Zlength)	
    elementsZ=Part.elements
    Part.Set(elements=elementsZ, name='all_solid_elementsZ')    
    Part.SetByBoolean(name='ShellZ', operation=DIFFERENCE, sets=(
    Part.sets['all_solid_elementsZ'], Part.sets['all_solid_elements'], Part.sets['all_solid_elementsY'],)) 


    
    # First lets calculate number of nodes on each face node set for negative side
    numNegX = Model.rootAssembly.sets['RVE_NegX'].nodes
    numX = len(numNegX)
    numNegY = Model.rootAssembly.sets['RVE_NegY'].nodes
    numY = len(numNegY)
    numNegZ = Model.rootAssembly.sets['RVE_NegZ'].nodes
    numZ = len(numNegZ)    

    
    ################################################################################################################################################
    # Now we have to create node sets again for copied mesh 
    # With same logic we will create all node sets again for copied mesh 
    # here suffix "Proj" will be added to all sets which stands for projected
    ##################################################################################################################################################
    
    nodesInShellX=Instance.sets['ShellX'].nodes    
    RVE_NegX_Proj = Model.rootAssembly.Set(name='RVE_NegX_Proj', nodes=nodesInShellX)
    nodesInShellY=Instance.sets['ShellY'].nodes    
    RVE_NegY_Proj = Model.rootAssembly.Set(name='RVE_NegY_Proj', nodes=nodesInShellY)
    nodesInShellZ=Instance.sets['ShellZ'].nodes    
    RVE_NegZ_Proj = Model.rootAssembly.Set(name='RVE_NegZ_Proj', nodes=nodesInShellZ)
    
   

    # using sorting function
    ##################################
    ###########  X- DIRECTION  ##########
    ##################################

    ##################################
    ###########  Y- DIRECTION  ##########
    ##################################
    #NegY_red_Proj = MapNodesetOnCopiedNodeset(Model,NegY,NegY_Proj,numY,NegY_red,'NegY_red_Proj')

    ##################################
    ###########  Z- DIRECTION  ##########
    ##################################

    #NegZ_red_Proj = MapNodesetOnCopiedNodeset(Model,NegZ,NegZ_Proj,numZ,NegZ_red,'NegZ_red_Proj')
    
    nodesForSet = MapSubsetFromCopiedNodes(Model.rootAssembly.sets['RVE_NegY'].nodes,
                                    Model.rootAssembly.sets['RVE_NegY_Proj'].nodes, 
                                    Model.rootAssembly.sets['RVE_NegY_red'].nodes)
    RVE_NegY_red_Proj = Model.rootAssembly.Set(name='RVE_NegY_red_Proj', nodes=nodesForSet)
    nodesForSet = MapSubsetFromCopiedNodes(Model.rootAssembly.sets['RVE_NegZ'].nodes,
                                    Model.rootAssembly.sets['RVE_NegZ_Proj'].nodes, 
                                    Model.rootAssembly.sets['RVE_NegZ_red'].nodes)
    RVE_NegZ_red_Proj = Model.rootAssembly.Set(name='RVE_NegZ_red_Proj', nodes=nodesForSet)

    # Assign all new elements to be surface elements and assigsection
    #Tuple of surface element types
    SurfElementTypes = ( mesh.ElemType(elemCode=DS3, elemLibrary=STANDARD),
                         mesh.ElemType(elemCode=DS4, elemLibrary=STANDARD),
                         mesh.ElemType(elemCode=DS6, elemLibrary=STANDARD),
                         mesh.ElemType(elemCode=DS8, elemLibrary=STANDARD) )
    AllSurfElements = Part.SetByBoolean(name="AllSurfElements", 
                                        operation=UNION, 
                                        sets=(Part.sets['ShellX'],Part.sets['ShellY'],Part.sets['ShellZ']))
    Part.setElementType(regions=AllSurfElements,
                        elemTypes=SurfElementTypes)
                        
    # Material creation to apply for shell elements in next step
    Model.Material(name='Material_SHELL')
    Model.materials['Material_SHELL'].Conductivity(table=((1e-06, ), ))    
                        
    SurfSec = Model.HomogeneousShellSection(name="ProjectedFaces", material="Material_SHELL", numIntPts=1, thickness=1e-05)
    Part.SectionAssignment(region=AllSurfElements, sectionName=SurfSec.name)
    # get all faces on NegX, NegY and NegZ using node sets we have    
    ###############################################################################################
    # TIE contraint
    # here we have to tie copied mesh surface to original mesh surface. i.e. MASTER = PosX, SLAVE = PosX_Projected
    ###############################################################################################

    
    # for X side
    elementsInShellX=Instance.sets['ShellX'].elements
    Assembly.Surface(side1Elements=elementsInShellX,name='Surf_PosX_Projected')
    surface=Assembly.surfaces['Surf_PosX_Projected']

    Surf4 = GetSurfaceFromNodeSet(Model,'RVE_PosX')
    includedSurface = Model.rootAssembly.Surface( name='PosXoffset', face1Elements=Surf4)
    regionX=Assembly.surfaces['PosXoffset']
    Model.Tie(name='X-Face Tie', master=regionX, 
        slave=surface, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, constraintEnforcement=SURFACE_TO_SURFACE, thickness=ON)


    # for Y side
    elementsInShellY=Instance.sets['ShellY'].elements
    Assembly.Surface(side1Elements=elementsInShellY,name='Surf_PosY_Projected')
    surface=Assembly.surfaces['Surf_PosY_Projected']

    Surf5 = GetSurfaceFromNodeSet(Model,'RVE_PosY')
    includedSurface = Model.rootAssembly.Surface( name='PosYoffset', face1Elements=Surf5)
    regionY=Assembly.surfaces['PosYoffset']
    Model.Tie(name='Y-Face Tie', master=regionY, 
        slave=surface, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, constraintEnforcement=SURFACE_TO_SURFACE, thickness=ON)


    # for Z side
    elementsInShellZ=Instance.sets['ShellZ'].elements
    Assembly.Surface(side1Elements=elementsInShellZ,name='Surf_PosZ_Projected')
    surface=Assembly.surfaces['Surf_PosZ_Projected']

    Surf6 = GetSurfaceFromNodeSet(Model,'RVE_PosZ')
    includedSurface = Model.rootAssembly.Surface( name='PosZoffset', face1Elements=Surf6)
    regionZ=Assembly.surfaces['PosZoffset']
    Model.Tie(name='Z-Face Tie', master=regionZ, 
        slave=surface, positionToleranceMethod=COMPUTED, adjust=ON, 
        tieRotations=ON, constraintEnforcement=SURFACE_TO_SURFACE, thickness=ON)            
    
    #--------------------------------------------------------------------------------------------------------------------------------------------------
    
    # After these sets are built, we need to sort each one based on position so that paired sets have matching order.
    # We can't set the sort order in a node set via the Python API, so we're going to make MeshNodeArray objects
    # with the proper sorting and then write the input file code directly.  An alternaitve to this would be to
    # create an individual equation for each pair of nodes, but since the API only supports putting sets into
    # equations, not only would this require lots of equation definitions, it would also require a node set for
    # every single node on the surface of the RVE, so we aren't going to do it that way.
    
    # The nodesets and equations need to be handled as custom keywords because of limitations in the Python API (no support for *NSET, UNSORTED)
    # We will store them in a string
    KeywordsString = ''
    
    # Define name to be used for sorted version of each set.
    DataForSorting = [ [RVE_NegX,             "RVE_NegX_ordered"],
                       [RVE_NegX_Proj,        "RVE_NegX_Proj_ordered"],
                       [RVE_NegY_red,         "RVE_NegY_red_ordered"],
                       [RVE_NegY_red_Proj,    "RVE_NegY_red_Proj_ordered"],
                       [RVE_NegZ_red,         "RVE_NegZ_red_ordered"],
                       [RVE_NegZ_red_Proj,    "RVE_NegZ_red_Proj_ordered"],]    

    #Do the sorting and write the nodesets
    SortedNodes = {}
    for nodeset, name in DataForSorting:
        SortedNodes[name] = GenerateSortedMeshNodeArray(nodeset.nodes)
        KeywordsString += WriteMeshNodeArrayToUnsortedNodeSetDef(SortedNodes[name],name)+'\n'
        
    #Create periodic reference nodes and nodesets pointing to the nodes.
    #Each node's DOF 11 will correspond to a different far-field temperature gradient.
    #We have to do it this way because Abaqus doesn't admit using DOFs reserved for
    #mechanical quantities for thermal
    ThermalRP(Model)
        
    # These are the pairs of faces/edges/vertices that need to be constrained.
    Pairs = [ [ "RVE_NegX_Proj_ordered",        "RVE_NegX_ordered" ],                    
              [ "RVE_NegY_red_Proj_ordered",    "RVE_NegY_red_ordered" ],                    
              [ "RVE_NegZ_red_Proj_ordered",    "RVE_NegZ_red_ordered" ],]
    # Now loop through and define the equations for each of these nodes.                       
    KeywordsString += "*EQUATION\n"
    for (masterName, slaveName) in Pairs:
        PeriodicityVec = GetPeriodicityVector(SortedNodes[masterName],SortedNodes[slaveName])
        if PeriodicityVec == None:
            print "ERROR - Could not match nodes on " + masterName + " and " + slaveName
            print "It is possible that the node positions on the faces do not match within the required tolerance"
            raise RuntimeError("Unmatched Faces")
        KeywordsString += DefineScalarFieldPeriodicEquation3D(masterName,
                                                              slaveName,
                                                              11,
                                                              [ ("RP-ThermGrad_x",11), ("RP-ThermGrad_y",11), ("RP-ThermGrad_z",11) ],
                                                              PeriodicityVec) + '\n'
    #AddEquationsAndUnsortedNsetStringToKeywords(Model,KeywordsString)                                                              
    return KeywordsString
    #We're done



# This is where the periodic equation gets defined.  This is the core kernel of this whole exercise.
def DefineScalarFieldPeriodicEquation3D( MasterNodesetName,
                                         SlaveNodesetName,
                                         SlaveMasterFieldDOF,
                                         FarFieldGradientDOFs,
                                         PeriodicityVector,
                                         tol=1e-8 ):
    """
    Defines equations to establish periodicity in a saclar field variable for a pair of nodesets.
    
    
    Equation form:    phi_s(x_j+p_j)=phi_m(x_j) + <del phi/del x_j> p_j
    p_j is the vector of periodicity, going from the master node to the slave node
    x_j is a coordinate
    phi_s(x_j+p_j) is the field value at the slave positioned at coordinate x_j+p_j
    phi_m(x_j) is the field value at the master positioned at coordiante x_j
    <del phi/del x_j> is the far-field (macro-scale) gradient of the field
    summation is implied over the index j=1..3
    
    Inputs:
    MasterNodesetName -    The master nodeset
    SlaveNodesetName -     The slave nodeset.  Must be ordered in a manner that corresponds with the master.
                           The degrees of freedom (DOFs) in these nodes will be the slaved DOFs.
    SlaveMasterFieldDOF -  The DOF for the field variable at the slave and master node
    FarFieldGradientDOFs - An iterable of pairs of the nodeset names and DOFs that define the far-field gradients.
                           For example, say that the scalar field is temperature, and you have 3 nodesets
                           "RP-ThermGrad_x", "RP-ThermGrad_y", and "RP-ThermGrad_z", each of which contains
                           one node whose DOF 11 represents dT/dx, dT/dy, and dT/dz, respectively.  Then you
                           would supply the following iterable to the function:
                           [ ("RP-ThermGrad_x",11), ("RP-ThermGrad_y",11), ("RP-ThermGrad_z",11) ]
    PeriodicityVector -    The vector that points from a node on the first nodeset to the second nodeset
    tol -                  Tolerance for determining whether a component of the periodicity vector, as a fraction
                           of the total periodicity vector length, is zero (which will result in fewer terms in 
                           the equation)
    
    Outputs:
    String defining all equations to pair nodesets
    """
    
    EqnDefs = []
    PeriodVecZeroThresh = math.sqrt(sum([a*a for a in PeriodicityVector]))*tol
    EqTerms = [ (-1.0,SlaveNodesetName,SlaveMasterFieldDOF), 
                ( 1.0,MasterNodesetName,SlaveMasterFieldDOF) ]
    for PVecComp,FarFieldGradVal in zip(PeriodicityVector, FarFieldGradientDOFs):
        if abs(PVecComp)>PeriodVecZeroThresh:
            EqTerms.append( (PVecComp,FarFieldGradVal[0],FarFieldGradVal[1]) )
    EqnDefs.append(WriteEquationTerms(EqTerms))
    return '\n'.join(EqnDefs)




def UniformSurfaceGradientThermalBC(Model):
    """
    UniformSurfaceGradientThermalBC:
    This method generates equations that will enforce UniformSurfaceGradient BC on a RVE.
    Creates reference nodes for specifying far-field temperature gradients.
    
    Inputs:
    Model -                 A model object from a Model Database

    
    Outputs:
    Returns the string containing the keywords to be written to the input deck.
    """

    CreateInitialSets(Model)
    TempCenter = 10
    #Define references to the sets that contain all the nodes on each face
    RVE_NegX = Model.rootAssembly.sets['RVE_NegX']
    RVE_PosX = Model.rootAssembly.sets['RVE_PosX']
    RVE_NegY = Model.rootAssembly.sets['RVE_NegY']
    RVE_PosY = Model.rootAssembly.sets['RVE_PosY']
    RVE_NegZ = Model.rootAssembly.sets['RVE_NegZ']
    RVE_PosZ = Model.rootAssembly.sets['RVE_PosZ']
    
    AllSurfaceNodes = Model.rootAssembly.SetByBoolean(name="AllSurfaceNodes", operation=assembly.UNION, sets=[ RVE_PosX, RVE_PosY, RVE_PosZ, RVE_NegX, RVE_NegY, RVE_NegZ ])    

    # We have to get coordinates of centre of RVE to apply equations 
    Edge_AlignedPosXNegY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXNegY", operation=assembly.INTERSECTION, sets=[ RVE_NegY, RVE_PosX ])
    Edge_AlignedPosXNegZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedPosXNegZ", operation=assembly.INTERSECTION, sets=[ RVE_NegZ, RVE_PosX ])
    Edge_AlignedNegXNegY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXNegY", operation=assembly.INTERSECTION, sets=[ RVE_NegY, RVE_NegX ])
    Edge_AlignedNegXPosY = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXPosY", operation=assembly.INTERSECTION, sets=[ RVE_PosY, RVE_NegX ])
    Edge_AlignedNegXNegZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXNegZ", operation=assembly.INTERSECTION, sets=[ RVE_NegZ, RVE_NegX ])
    Edge_AlignedNegXPosZ = Model.rootAssembly.SetByBoolean(name="Edge_AlignedNegXPosZ", operation=assembly.INTERSECTION, sets=[ RVE_PosZ, RVE_NegX ])    

    Vertex_PosXNegYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_PosXNegYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedPosXNegY, Edge_AlignedPosXNegZ ])
    Vertex_NegXPosYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXPosYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXPosY, Edge_AlignedNegXNegZ ])
    Vertex_NegXNegYPosZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXNegYPosZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXNegY, Edge_AlignedNegXPosZ ])
    Vertex_NegXNegYNegZ = Model.rootAssembly.SetByBoolean(name="Vertex_NegXNegYNegZ", operation=assembly.INTERSECTION, sets=[ Edge_AlignedNegXNegY, Edge_AlignedNegXNegZ ])    

    V0=	Model.rootAssembly.sets['Vertex_NegXNegYNegZ'].nodes[0].coordinates
    V1= Model.rootAssembly.sets['Vertex_PosXNegYNegZ'].nodes[0].coordinates
    V2=	Model.rootAssembly.sets['Vertex_NegXPosYNegZ'].nodes[0].coordinates
    V3= Model.rootAssembly.sets['Vertex_NegXNegYPosZ'].nodes[0].coordinates
	
    ## here lengths are kept with sigh so that center will be correct
    Xlength = (V0[0]-V1[0])          # length of RVE in X direction
    Ylength = (V0[1]-V2[1])			# length of RVE in Y direction
    Zlength = (V0[2]-V3[2])			# length of RVE in Z direction  

    RPcenterX = V0[0]+(Xlength/2)
    RPcenterY = V0[1]+(Ylength/2)
    RPcenterZ = V0[2]+(Zlength/2)
    
    RVECenter = (RPcenterX, RPcenterY, RPcenterZ)

    ThermalRP(Model)
    #One extra RP for just UniformSurfaceGradient 
    ThermalCenter_UniformSurfaceGradientRefNode = Model.rootAssembly.Node( (0,0,0) )
    Model.rootAssembly.Set("RP-ThermCenter", mesh.MeshNodeArray([ThermalCenter_UniformSurfaceGradientRefNode]))

    
    for i in AllSurfaceNodes.nodes:
        x1 = i.coordinates[0]
        x2 = i.coordinates[1]
        x3 = i.coordinates[2]

    KeywordsString = ''
    KeywordsString += "*EQUATION\n"
    for slaveNode in AllSurfaceNodes.nodes:
        KeywordsString += DefineUniformSurfaceGradientEquation3DThermal(slaveNode,
                                                        "RP-ThermGrad_x",
                                                        "RP-ThermGrad_y",
                                                        "RP-ThermGrad_z",
                                                        "RP-ThermCenter",
                                                        TempCenter,
                                                        RVECenter) + '\n'

    #AddEquationsAndUnsortedNsetStringToKeywords(Model,KeywordsString) 
 
    return KeywordsString    

def DefineUniformSurfaceGradientEquation3DThermal( SlaveNode,
                                   ThermGrad_x,
                                   ThermGrad_y,
                                   ThermGrad_z,
                                   ThermCenter,
                                   TempCenter,
                                   RVECenter,
                                   tol=1e-8 ):
    """
    DefineUniformSurfaceGradientEquation3DThermal:
    Creates equation objects in the model that cause this node's temperature to conform to a uniform temperature gradient
    field specified through the degrees of freedom of three reference nodes.
    
    Inputs:
    SlaveNode -            The slave node - an Abaqus MeshNode object
                           The DOFs in this node will be the slaved to the global temperature gradient field
    ThermGrad_x -  The nodeset name of the node whose DOF 11 are the far-field temperature gradient in x direction
    ThermGrad_y -  The nodeset name of the node whose DOF 11 are the far-field temperature gradient in y direction
    ThermGrad_z -   The nodeset name of the node whose DOF 11 are the far-field temperature gradient in y direction                      

    tol -                  Tolerance for determining whether a coordinate is zero.  If the coordinate's absolute
                           value is less than this fraction of the RVEs dimension, it is considered zero and no
                           entry is written to the equation for this term.
                           
    Outputs:
    String defining all equations to constrain the node to the global temperature field.
    """
    
    EqnDefs = []
    slaveName = ''
    if not SlaveNode.instanceName == None:
        slaveName = SlaveNode.instanceName + '.' + str(SlaveNode.label)
    else:
        slaveName = str(SlaveNode.label)    
    EqTerms = [ (-1.0,slaveName,11) ]
    if (TempCenter != None):
        EqTerms.append( (1.0,ThermCenter,11) )  
    
    EqTerms.append( ((SlaveNode.coordinates[0]-RVECenter[0]),ThermGrad_x,11) )
    EqTerms.append( ((SlaveNode.coordinates[1]-RVECenter[1]),ThermGrad_y,11) )
    EqTerms.append( ((SlaveNode.coordinates[2]-RVECenter[2]),ThermGrad_z,11) )
    EqnDefs.append(WriteEquationTerms(EqTerms))
    
    
    # Maps strain components to NormNode and ShearNode DOFs
    #ijmap = [ [ (NormNodeNodesetName,1),  (ShearNodeNodesetName,1), (ShearNodeNodesetName,2) ],
    #          [ (ShearNodeNodesetName,1), (NormNodeNodesetName,2),  (ShearNodeNodesetName,3) ],
    #          [ (ShearNodeNodesetName,2), (ShearNodeNodesetName,3), (NormNodeNodesetName,3)  ] ]
    # A little code to handle a node that may or may not be in an instance...
    #slaveName = ''
    #if not SlaveNode.instanceName == None:
    #    slaveName = SlaveNode.instanceName + '.' + str(SlaveNode.label)
    #else:
    #    slaveName = str(SlaveNode.label)
    #remember that Python is zero-based, so add 1 whenever using i or j to define DOFs
    #for i in range(3):
    #    EqTerms = [ (-1.0,slaveName,i+1) ]
    #    for j in range(3):
    #        if abs(SlaveNode.coordinates[j]/ModelDimensions[j])>tol:
    #            (nset,dof) = ijmap[i][j]
    #            EqTerms.append( (SlaveNode.coordinates[j],nset,dof) )
    #    EqnDefs.append(WriteEquationTerms(EqTerms))
    return '\n'.join(EqnDefs)
    
    

if __name__=="__main__":
    # Executive for running function on a CAE model.
    applyThermalPBC( mdb.models['Model-1'])

