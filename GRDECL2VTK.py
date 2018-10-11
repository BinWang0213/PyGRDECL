#########################################################################
#       (C) 2017 Department of Petroleum Engineering,                   # 
#       Univeristy of Louisiana at Lafayette, Lafayette, US.            #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# PyGRDECL Code                                                         #
# Author: Bin Wang                                                      # 
# Email: binwang.0213@gmail.com                                         # 
#########################################################################

import numpy as np
import os
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning) #Supress some Future Warnning due to outdated library

try:
    import vtk
    import vtk.util.numpy_support as ns
except ImportError:
    warnings.warn("No vtk module loaded.")



from GRDECL_Parser import *
from GRDECL_FaultProcess import *
#from GRDECL_CADExporter import *


#############################################
#
#  Core Eclipse File Read/View/Analysis Class
#
#############################################

class GeologyModel:
    def __init__(self,filename=''):
        """Eclipse Input file(GRDECL) Visulazation and Analysis  
        Keywords Reference: file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ         -- Grid dimension.
        Trans(i01,j01,k01) -- Transmisability in i,j,k direction
        fault(i01,j01)     -- Fault indicator in i,j direction(0-sealing, 0.5-partially connecting, 1-fully connecting)
        GRID_type           - 0-Cartesian 1-Corner point

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        self.fname=filename
        self.GRDECL_Data=GRDECL_Parser()
        self.FaultProcessor=None

        self.VTK_Grids=vtk.vtkUnstructuredGrid()


        if(self.fname!=''):
            self.GRDECL_Data=GRDECL_Parser(self.fname)

    def readGRDECL(self,filename):
        #* Using GRDECL parser to load the raw model data, DX,DY,COORD,etc
        self.GRDECL_Data.fname=filename
        self.GRDECL_Data.read_GRDECL()

    def GRDECL2VTK(self):
        #* Convert corner point grid/cartesian grid into VTK unstructure grid
        print('[Geometry] Converting GRDECL to Paraview Hexahedron mesh data....')
        NX,NY,NZ=self.GRDECL_Data.NX,self.GRDECL_Data.NY,self.GRDECL_Data.NZ
        if(self.GRDECL_Data.GRID_type=='CornerPoint'):
            
            #1.Collect Points from the raw CornerPoint data [ZCORN]&[COORD]
            #X,Y has to be interpolated from [ZCORN]
            Points = vtk.vtkPoints()
            Points.SetNumberOfPoints(len(self.GRDECL_Data.ZCORN)) #=2*NX*2*NY*2*NZ
            
            ptsid=0
            for k in range(NZ):
                for j in range(NY):
                    for i in range(NX):
                        CellCoords=self.GRDECL_Data.getCellCoords(i,j,k)
                        for pi in range(8): # Loop 8 point for each cell, see getCellCoords(i,j,k) for node ordering
                            Points.SetPoint(ptsid,CellCoords[pi])
                            ptsid+=1
            self.VTK_Grids.SetPoints(Points)
            
            #2. Recover Cells which follows the convention of [ZCORN]
            cellArray = vtk.vtkCellArray()
            Cell=vtk.vtkHexahedron()

            cellid=0
            for k in range(NZ):
                for j in range(NY):
                    for i in range(NX):
                        for pi in range(8):
                            #Convert GRDECL node index convention to VTK convention
                            #https://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
                            #0,1,2,3(GRDECL)->0,1,3,2(VTK,anti-clockwise)
                            if(pi==2 or pi==6): VTKid=pi+1
                            elif(pi==3 or pi==7): VTKid=pi-1 
                            else: VTKid=pi
                            Cell.GetPointIds().SetId(pi,cellid*8+VTKid)
                        cellArray.InsertNextCell(Cell)
                        cellid+=1
            
            self.VTK_Grids.SetCells(Cell.GetCellType(), cellArray)

            print("     NumOfPoints",self.VTK_Grids.GetNumberOfPoints())
            print("     NumOfCells",self.VTK_Grids.GetNumberOfCells())

            #3. Load grid properties data if applicable
            for keyword,data in self.GRDECL_Data.SpatialDatas.items():
                self.AppendScalarData(keyword,data)

            print('.....Done!')
    
    def decomposeModel(self):
        '''#* Identify and extract boundary/falut faces
        
        Fault-based model decomposition,subdividing the geology model along fault face
        **Currently, fault only happens on X,Y plane. No fault in Z direction

          6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1         

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        '''

        if(self.FaultProcessor is None):
            #1. Find the fault line and split the domain
            self.FaultProcessor=FaultProcess(self.GRDECL_Data)
            self.FaultProcessor.findFaultLines()
            self.FaultProcessor.findBoundaryLines()
            self.FaultProcessor.SplitDomainByFault()

        #2. Convert Splitted Domain (2D) into 3D sub volume
        #Convert 2D sub domain into 3D subvolume
        CellCenter=[(i+0.5,j+0.5) for j in range(self.GRDECL_Data.NY) for i in range(self.GRDECL_Data.NX)]

        DomainMarker2D=np.zeros(self.GRDECL_Data.NY*self.GRDECL_Data.NX)
        for i,poly in enumerate(self.FaultProcessor.SplitPolygons):
            #print(len(poly),poly)
            #mark the (i,j) as value of [flag] if cell center is in polygon
            flag=points_in_polygon(CellCenter,poly,flag=i+1) 
            #https://stackoverflow.com/questions/31862704/numpy-combine-all-nonzero-elements-of-one-array-in-to-another
            DomainMarker2D = np.where(flag == 0, DomainMarker2D, flag)
        
        DomainMarker3D=np.tile(DomainMarker2D,self.GRDECL_Data.NZ)
        self.AppendScalarData('SubVolumeIDs',DomainMarker3D)

        DomainMarker2D=np.zeros(self.GRDECL_Data.NY*self.GRDECL_Data.NX)
        RandomColor=10*np.random.rand(len(self.FaultProcessor.SplitPolygons))
        for i,poly in enumerate(self.FaultProcessor.SplitPolygons):
            #print(len(poly),poly)
            #mark the (i,j) as value of [flag] if cell center is in polygon
            flag=points_in_polygon(CellCenter,poly,flag=RandomColor[i]) 
            #https://stackoverflow.com/questions/31862704/numpy-combine-all-nonzero-elements-of-one-array-in-to-another
            DomainMarker2D = np.where(flag == 0, DomainMarker2D, flag)
        
        DomainMarker3D=np.tile(DomainMarker2D,self.GRDECL_Data.NZ)
        self.AppendScalarData('SubVolumes',DomainMarker3D)
    
    def AppendScalarData(self,name,numpy_array):
        #* Append scalar cell data (numpy array) into vtk object 
        data = ns.numpy_to_vtk(numpy_array.ravel(order='F'),deep=True, array_type=vtk.VTK_FLOAT)
        data.SetName(str(name))
        data.SetNumberOfComponents(1)
        self.VTK_Grids.GetCellData().AddArray(data)

    def Write2VTU(self):
        basename=os.path.splitext(os.path.basename(self.fname))[0]
        if not os.path.exists("Results"):
            os.makedirs('Results')
        path=os.path.join('Results',basename + '.vtu')
        print('[Output] Writing "%s" Paraview file....'%(path),end='')

        xmlWriter = vtk.vtkXMLUnstructuredGridWriter()
        xmlWriter.SetFileName(path)
        xmlWriter.SetInputData(self.VTK_Grids)
        xmlWriter.Write()
        print('Done!')


                        



    