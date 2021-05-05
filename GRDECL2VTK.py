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
# Contributing Author: Mustapha Zakari (MZ)                             #
# Email: mustapha.zakari@univ-lorraine.fr                               #
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
from utils import *

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

        self.VTK_Grids=None

        if(self.fname!=''):
            self.GRDECL_Data=GRDECL_Parser(self.fname)

    def readGRDECL(self,filename):
        #* Using GRDECL parser to load the raw model data, DX,DY,COORD,etc
        self.GRDECL_Data.fname=filename
        self.GRDECL_Data.read_GRDECL()

    def buildCartGrid(self,physDims=[100.0,100.0,10.0],gridDims=[10,10,1]):
        #* Create simple cartesian grid
        self.GRDECL_Data.buildCartGrid(physDims,gridDims)

    def GRDECL2VTK(self):
        #* Convert corner point grid/cartesian grid into VTK unstructure grid
        print('[Geometry] Converting GRDECL to Paraview Hexahedron mesh data....')
        NX,NY,NZ=self.GRDECL_Data.NX,self.GRDECL_Data.NY,self.GRDECL_Data.NZ
        self.VTK_Grids=vtk.vtkUnstructuredGrid()

        if(self.GRDECL_Data.GRID_type=='Cartesian'):
            #1. Collect points from Cartesian data
            debug=0
            DX,DY,DZ,TOPS=self.GRDECL_Data.DX,self.GRDECL_Data.DY,self.GRDECL_Data.DZ,self.GRDECL_Data.TOPS
            coordX,coordY,coordZ=Cartesian2UnstructGrid(DX,DY,DZ,TOPS,NX,NY,NZ)

            Points = vtk.vtkPoints()
            Points.SetNumberOfPoints(2*NX*2*NY*2*NZ) #=2*NX*2*NY*2*NZ

            ptsid=0
            for k in range(2*NZ):
                for j in range(2*NY):
                    for i in range(2*NX):
                        x,y,z=coordX[i][j][k],coordY[i][j][k],coordZ[i][j][k]
                        Points.SetPoint(ptsid,[x,y,z])
                        ptsid+=1
            self.VTK_Grids.SetPoints(Points)

            #2. Recover Cells which follows the cartesian data
            cellArray = vtk.vtkCellArray()
            Cell=vtk.vtkHexahedron()

            cellid=0
            for k in range(NZ):
                for j in range(NY):
                    for i in range(NX):
                        idx_GB=[getIJK(2*i,2*j,2*k,2*NX,2*NY,2*NZ),  #index-0
                            getIJK(2*i+1,2*j,2*k,2*NX,2*NY,2*NZ),    #index-1
                            getIJK(2*i+1,2*j+1,2*k,2*NX,2*NY,2*NZ),  #index-2
                            getIJK(2*i,2*j+1,2*k,2*NX,2*NY,2*NZ),    #index-3
                            getIJK(2*i,2*j,2*k+1,2*NX,2*NY,2*NZ),    #index-4
                            getIJK(2*i+1,2*j,2*k+1,2*NX,2*NY,2*NZ),  #index-5
                            getIJK(2*i+1,2*j+1,2*k+1,2*NX,2*NY,2*NZ),#index-6
                            getIJK(2*i,2*j+1,2*k+1,2*NX,2*NY,2*NZ)]  #index-7
                        for pi in range(8):
                            Cell.GetPointIds().SetId(pi,idx_GB[pi])
                        cellArray.InsertNextCell(Cell)
                        cellid+=1

            self.VTK_Grids.SetCells(Cell.GetCellType(), cellArray)


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

        self.Update()

        print('     .....Done!')

    def Update(self):
        #Load all available keywords/cellarrays into VTK container
        for keyword,data in self.GRDECL_Data.SpatialDatas.items():
            self.AppendScalarData2VTK(keyword,data) #VTK will automatically overwrite the data with the same keyword


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
        self.AppendScalarData2VTK('SubVolumeIDs',DomainMarker3D)

        DomainMarker2D=np.zeros(self.GRDECL_Data.NY*self.GRDECL_Data.NX)
        RandomColor=10*np.random.rand(len(self.FaultProcessor.SplitPolygons))
        for i,poly in enumerate(self.FaultProcessor.SplitPolygons):
            #print(len(poly),poly)
            #mark the (i,j) as value of [flag] if cell center is in polygon
            flag=points_in_polygon(CellCenter,poly,flag=RandomColor[i])
            #https://stackoverflow.com/questions/31862704/numpy-combine-all-nonzero-elements-of-one-array-in-to-another
            DomainMarker2D = np.where(flag == 0, DomainMarker2D, flag)

        DomainMarker3D=np.tile(DomainMarker2D,self.GRDECL_Data.NZ)
        self.AppendScalarData2VTK('SubVolumes',DomainMarker3D)


    def CreateCellData(self,varname="SW",val=0.0,val_array=[]):
        """Create a new data field

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """
        assert varname not in self.GRDECL_Data.SpatialDatas, "     [Error] Variable [%s] is existed! Please use UpdateCellData function" %(varname)

        if(len(val_array)==0):
            self.GRDECL_Data.SpatialDatas[varname]=np.ones(self.GRDECL_Data.N)*val
            print('     New variable [%s] created with a value of %lf!'%(varname,val))
        else:
            assert len(val_array)==self.GRDECL_Data.N, print('     [Error] Input array is not compatible with number of cells!')
            self.GRDECL_Data.SpatialDatas[varname]=np.array(val_array)
            print('     New variable [%s] created with a given array!'%(varname,val))

    def LoadCellData(self,varname="SW",filename="123.txt"):
        """Create a new data field and load from a file

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """
        Data=KeyWordReader(filename,varname)
        if(Data is not None):
            assert len(Data)==self.GRDECL_Data.N, print('     [Error] Input array is not compatible with number of cells!')
            self.GRDECL_Data.SpatialDatas[varname]=Data
            print('     New variable [%s] loaded from file!'%(varname))

        return Data


    def UpdateCellData(self,varname="PERMX",val=100.0,nx_range=(1,-1),ny_range=(1,-1),nz_range=(1,-1),array=[]):
        """Update/modify Cell data field (Permeability/Porosity) with given grid block range

        Arguments
        ---------
        var         -- The varable name you want to update, e.g PERMX, PERMY or PERMZ
        val         -- The varable value you want to update
        nx_range    -- The specifc grid range in x for updating, 1-based index
        nx_range    -- The specifc grid range in y for updating, 1-based index
        nz_range    -- The specifc grid range in z for updating, 1-based index

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """

        if (varname not in self.GRDECL_Data.SpatialDatas): #New cell data
            print("[Warnning] Variable [%s] is not existed!")

        nx_range=np.array(nx_range)
        ny_range=np.array(ny_range)
        nz_range=np.array(nz_range)

        #If no nx,ny,nz range are defined, all perm will be updated
        if(nx_range[1] == -1):
            nx_range[1] = self.GRDECL_Data.NX
        if(ny_range[1] == -1):
            ny_range[1] = self.GRDECL_Data.NY
        if(nz_range[1] == -1):
            nz_range[1] = self.GRDECL_Data.NZ

        #Convert 1-based index to 0-based index
        nx_range=nx_range-1
        ny_range=ny_range-1
        nz_range=nz_range-1

        #Update perm field with specific grid range
        for k in range(self.GRDECL_Data.NZ):
            for j in range(self.GRDECL_Data.NY):
                for i in range(self.GRDECL_Data.NX):
                    if(i>=nx_range[0] and i<=nx_range[1]):
                        if(j>=ny_range[0] and j<=ny_range[1]):
                            if(k>=nz_range[0] and k<=nz_range[1]):
                                ijk = getIJK(i, j, k, self.GRDECL_Data.NX, self.GRDECL_Data.NY, self.GRDECL_Data.NZ)

                                if(len(array)>0):#We load up value from a 1D array which storaged in the following loop order
                                    val=array[ijk]
                                self.GRDECL_Data.SpatialDatas[varname][ijk]=val


    def WriteNPSL(self):
        """Write the permeability/porosity field for NPSL

            filename_permx.txt
            filename_permy.txt
            filename_permz.txt
            filename_poro.txt

            Arguments
            ---------
            filename    -- The surname of permeability and porosity data files

            Programmer: Bin Wang (yin.feng@louisiana.edu)
            Creation:   Feb, 2018
        """
        basename=os.path.splitext(os.path.basename(self.fname))[0]
        if not os.path.exists("Results"):
            os.makedirs('Results')
        fnames=[os.path.join('Results',basename + '_permx.txt'),
               os.path.join('Results',basename + '_permy.txt'),
               os.path.join('Results',basename + '_permz.txt'),
               os.path.join('Results',basename + '_poro.txt')]

        np.savetxt(fnames[0], self.GRDECL_Data.SpatialDatas['PERMX'], delimiter="\n",fmt='%1.4f')
        np.savetxt(fnames[1], self.GRDECL_Data.SpatialDatas['PERMY'], delimiter="\n",fmt='%1.4f')
        np.savetxt(fnames[2], self.GRDECL_Data.SpatialDatas['PERMZ'], delimiter="\n",fmt='%1.4f')
        np.savetxt(fnames[3],  self.GRDECL_Data.SpatialDatas['PORO'], delimiter="\n",fmt='%1.4f')

        for name in fnames:
            print('NPSL file [%s] successfully genetrated, pelase use NPSL to load it!' % (name))

        #Output special NPSL init sw field
        if('SW_NPSL' in self.GRDECL_Data.SpatialDatas):
            fnames.append(os.path.join('Results',basename + '_sw.txt'))
            header='#CheckPoint_Data\nTIMESTEP 1\nCELL_DATA 3600\nNUMBER_OF_SPECIES 1\n\nSCALARS C_PseudoOil float'
            np.savetxt(fnames[-1],  self.GRDECL_Data.SpatialDatas['SW_NPSL'], delimiter="\n",fmt='%1.4f',header=header,comments="")
            print('NPSL file [%s] successfully genetrated, pelase use NPSL to load it!' % (fnames[-1]))



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

    def Write2VTP(self):
        self.Write2VTU()
        basename=os.path.splitext(os.path.basename(self.fname))[0]
        inputFile=os.path.join('Results',basename + '.vtu')
        outFile=os.path.join('Results',basename + '.vtp')
        print('[Output] Writing "%s" VTP file..'%(outFile),end='')
        reader = vtk.vtkXMLUnstructuredGridReader()
        reader.SetFileName(inputFile)
        reader.Update()
        ugrid = reader.GetOutput()
        geometryFilter = vtk.vtkGeometryFilter()
        geometryFilter.SetInputData(ugrid)
        geometryFilter.Update()
        polydata = geometryFilter.GetOutput()
        writer =vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(outFile)
        writer.SetInputData(polydata)
        writer.Write()
        print("vtp file created.")




    def AppendScalarData2VTK(self,name,numpy_array):
        #* Append scalar cell data (numpy array) into vtk object, should not directly called by user
        data = ns.numpy_to_vtk(numpy_array.ravel(order='F'),deep=True, array_type=vtk.VTK_FLOAT)
        data.SetName(str(name))
        data.SetNumberOfComponents(1)
        self.VTK_Grids.GetCellData().AddArray(data)

    #####################################################################
    # MZ::GRDECL cartesian writer functions
    # MZ::Write spatial data at a specified output format:
    # MZ::Writes rows "linestring" at given (j,k) and varying i indices
    #####################################################################
    def writeCartGrid_forGRDECL(self,filename):
        # MZ:: Dictionary keyword+data storing Grid attributes to be written
        Grid=self.GRDECL_Data;
        self.grid_Data = {
            "TOPS":Grid.TOPS,
            "DX": Grid.DX,"DY": Grid.DY,"DZ": Grid.DZ
        }

        # Output format for floats
        floatformat = "%5.2f"

        # Open file  filename path
        self.fname=filename;
        f = open(self.fname, 'w')
        print('[Output] Writing "%s" GRDECL file..'%filename,end='')

        # Write DIMENSIONS
        f.write("GRID\n\n")
        f.write("DIMENS\n")
        f.write(str(Grid.NX)+" "+str(Grid.NY)+" "+str(Grid.NZ)+" /"+"\n")

        # Write Grid&Spatial information using merged dictionaries,
        for key,value in {**self.grid_Data, **self.GRDECL_Data.SpatialDatas}.items():
            # print(" Writing %s"%key)
            f.write("\n"+key)
            write_grid_data_forGRDECL(f, Grid, value, floatformat)
            f.write("\n")

        # Close file
        f.close()
        print("...done")

#####################################################################
# MZ::GRDECL cartesian writer functions
# MZ::Write spatial data at a specified output format:
# MZ::Writes rows "linestring" at given (j,k) and varying i indices
#####################################################################
def write_grid_data_forGRDECL(f,Grid,data,data_format):
    # Add space between scalar values
    data_format=data_format+" "
    for k in range(Grid.NZ):
        for j in range(Grid.NY):
            # Set local flat indices range for i in (0,NX-1)
            i0 = 0 + (Grid.NX) * (j + k * (Grid.NY))
            iend = Grid.NX - 1 + (Grid.NX) * (j + k * (Grid.NY))

            # Write First Elt of linestring at new line
            linestring="\n"+custom_format(data_format,data[i0])
            # Append next ELts: from i0+1 to iend
            for i in range(i0+1, iend+1):
                linestring=linestring+ custom_format(data_format,data[i])

            # Find duplicates and merge them in string "nbelt*Elt"
            linestring=merge_duplicates_forGRDECL(linestring)

            f.write("\n"+linestring)
    f.write("/")

# MZ::Change output format to "%d " for integer values
def custom_format(data_format,data):
    if (int(data)==data):
        data_format= "%d "
    return data_format % data

# MZ::Find duplicates and merge them in string "nbelt*Elt"
from itertools import groupby
def merge_duplicates_forGRDECL(linestring):
    new_str = ""
    for k, v in groupby(linestring.split()):
        group = list(v)
        if (len(group) > 1): new_str += str(len(group)) + "*"
        new_str += group[0] + " "
    return new_str






