#########################################################################
#       (C) 2017 Department of Petroleum Engineering,                   # 
#       Univeristy of Louisiana at Lafayette, Lafayette, US.            #
#                                                                       #
# This code is released under the terms of the BSD license, and thus    #
# free for commercial and research use. Feel free to use the code into  #
# your own project with a PROPER REFERENCE.                             #
#                                                                       #
# PyGRDECL Code                                                          #
# Author: Bin Wang                                                      # 
# Email: binwang.0213@gmail.com                                         # 
#########################################################################

import numpy as np
from pyvtk import *

#############################################
#
#  Core Eclipse File Read/View/Analysis Class
#
#############################################

class GRDECL_Viewer:
    def __init__(self,filename,nx,ny,nz):
        """Eclipse Input file(GRDECL) Visulazation and Analysis  
        Keywords Reference: file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ         -- Grid dimension.
        Trans(i01,j01,k01) -- Transmisability in i,j,k direction
        fault(i01,j01)     -- Fault indicator in i,j direction(0-sealing, 0.5-partially connecting, 1-fully connecting)
        GRID_type           - 0-Cartesian 1-Corner point

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        self.fname=filename
        self.NX=nx
        self.NY=ny
        self.NZ=nz
        self.N=nx*ny*nz
        self.GRID_type=0

        #Cartesian gridblock data KeyWords
        self.TOPS=[]
        self.DX=[]
        self.DY=[]
        self.DZ=[]

        #Corner point gridblock data KeyWrods (not support now)
        self.COORD=[]
        self.ZCORN=[]

        #Petrophysics data Keywords
        self.PERMX=[]
        self.PERMY=[]
        self.PERMZ=[]
        self.PORO=[]

        #Read GRDECL file when initializing the class
        if(len(filename)>0):
            self.read_GRDECL()

        #Derived variabls 
        self.Trans_i0,self.Trans_i1=np.zeros(self.N),np.zeros(self.N)
        self.Trans_j0,self.Trans_j1=np.zeros(self.N),np.zeros(self.N)
        self.Trans_k0,self.Trans_k1=np.zeros(self.N),np.zeros(self.N)

        self.fault_i0,self.fault_i1=np.zeros(self.N),np.zeros(self.N)
        self.fault_j0,self.fault_j1=np.zeros(self.N),np.zeros(self.N)

    def print_self(self):
        print('Dimension(NX,NY,NZ): (%s X %s X %s)'%(self.NX,self.NY,self.NZ))
        print('NumofGBs: %s'%(self.N))
    
    def read_GRDECL(self):
        """Read input file(GRDECL) of Reservoir Simulator- Petrel (Eclipse)  
        file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ -- Grid dimension.
        blockData_raw -- [0] Keywords [1] values
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2017
        """
        debug=0

        f=open(self.fname)
        contents=f.read()
        contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
        
        for i,block in enumerate(contents_in_block):#Block-wise
            blockData_raw=block.strip().split()
            block_id=0
            
            block_dataset=[]
            block_id=0
            if(len(blockData_raw)>0):
                if(".dat" not in blockData_raw[1]):
                    block_id,block_dataset=getBlkdata(blockData_raw,self.N)
                    if(debug):
                        print('block',i)
                        print(block_id,block_dataset)
                
                    block_dataset=np.array(block_dataset,dtype=float)
                else:
                    block_id=block_indicator(blockData_raw[0])
                    block_dataset=self.read_IncludeFile('./Example/'+blockData_raw[1],self.N)
                        #print(block_dataset.shape,block_dataset)

            #Assign read dataset to corrsponding Variables(DX,DY...etc)
                if(block_id==1):
                    self.TOPS=block_dataset
                if(block_id==2):
                    self.DX=block_dataset
                if(block_id==3):
                    self.DY=block_dataset
                if(block_id==4):
                    self.DZ=block_dataset
                if(block_id==5):
                    self.PORO=block_dataset
                if(block_id==6):
                    self.PERMX=block_dataset
                if(block_id==7):
                    self.PERMY=block_dataset
                if(block_id==8):
                    self.PERMZ=block_dataset
    
    def read_IncludeFile(self,filename_include,NumData):
        """Read Include data file
        this data file just a series of values
        e.g. 0.2 0.3 12.23 ....
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Aug. 2018
        """

        f=open(filename_include)
        contents=f.read()
        block_dataset=contents.strip().split() #Sepeart input file by slash /
        block_dataset=np.array(block_dataset,dtype=float)
        if(len(block_dataset)!=NumData):
            print('Data size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(len(block_dataset),NumData))
        return block_dataset

    def field_cutter(self,nx_range=(0,-1),ny_range=(0,-1),nz_range=(0,-1)):
        """Extract the subset of a domain
        
        Arguments
        ---------
        nx_range    -- The specifc grid range in x for the subset 
        nx_range    -- The specifc grid range in y for the subset 
        nz_range    -- The specifc grid range in z for the subset 
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """

        #If no nx,ny,nz range are defined, all perm will be updated       
        if(nx_range[1]==-1):
            nx_range[1] = self.NX
        if(ny_range[1] == -1):
            ny_range[1] = self.NY
        if(nz_range[1]==-1):
            nz_range[1] = self.NZ

        NX_new=nx_range[1]
        NY_new=ny_range[1]
        NZ_new=nz_range[1]
        N_new=NX_new*NY_new*NZ_new

        TOPS_new=np.zeros(N_new)
        DX_new=np.zeros(N_new)
        DY_new=np.zeros(N_new)
        DZ_new=np.zeros(N_new)

        PERMX_new=np.zeros(N_new)
        PERMY_new=np.zeros(N_new)
        PERMZ_new=np.zeros(N_new)
        PORO_new=np.zeros(N_new)

        ijk_new=0
        for k in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    if(i>=nx_range[0] and i<=nx_range[1]):
                        if(j>=ny_range[0] and j<=ny_range[1]):
                            if(k>=nz_range[0] and k<=nz_range[1]):
                                ijk = getIJK(i, j, k, self.NX, self.NY, self.NZ)
                                DX_new[ijk_new]=self.DX[ijk]
                                DY_new[ijk_new]=self.DY[ijk]
                                DZ_new[ijk_new]=self.DZ[ijk]
                                TOPS_new[ijk_new]=self.TOPS[ijk]

                                PERMX_new[ijk_new]=self.PERMX[ijk]
                                PERMY_new[ijk_new]=self.PERMY[ijk]
                                PERMZ_new[ijk_new]=self.PERMZ[ijk]
                                PORO_new[ijk_new]=self.PORO[ijk]
                                ijk_new=ijk_new+1
        
        NewGrid=GRDECL_Viewer('',NX_new,NY_new,NZ_new)
        NewGrid.DX=DX_new
        NewGrid.DY=DY_new
        NewGrid.DZ=DZ_new
        NewGrid.TOPS=self.TOPS
        NewGrid.PERMX=PERMX_new
        NewGrid.PERMY=PERMY_new
        NewGrid.PERMZ=PERMZ_new
        NewGrid.PORO=PORO_new
        NewGrid.print_self()
        return NewGrid


    def field_update(self,var="PERMX",val=100.0,nx_range=(0,-1),ny_range=(0,-1),nz_range=(0,-1)):
        """Update/modify Permeability/Porosity field with given grid block range
        
        Arguments
        ---------
        var         -- The varable name you want to update, e.g PERMX, PERMY or PERMZ
        val         -- The varable value you want to update
        nx_range    -- The specifc grid range in x for updating 
        nx_range    -- The specifc grid range in y for updating
        nz_range    -- The specifc grid range in z for updating
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2018
        """
        nx_range=np.array(nx_range)
        ny_range=np.array(ny_range)
        nz_range=np.array(nz_range)
        #If no nx,ny,nz range are defined, all perm will be updated       
        if(nx_range[1]==-1):
            nx_range[1] = self.NX
        if(ny_range[1] == -1):
            ny_range[1] = self.NY
        if(nz_range[1]==-1):
            nz_range[1] = self.NZ

        #Update perm field with specific grid range
        for k in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    if(i>=nx_range[0] and i<=nx_range[1]):
                        if(j>=ny_range[0] and j<=ny_range[1]):
                            if(k>=nz_range[0] and k<=nz_range[1]):
                                ijk = getIJK(i, j, k, self.NX, self.NY, self.NZ)
                                if(var=="PERMX"):
                                    self.PERMX[ijk]=val
                                if(var=="PERMY"):
                                    self.PERMY[ijk]=val
                                if(var=="PERMZ"):
                                    self.PERMZ[ijk]=val
                                if(var=="PORO"):
                                    self.PORO[ijk]=val
        
    def write_Ascii(self,filehead):
        """Output Perm Field Into several files
        permx.dat
        permy.dat
        permz.dat
        poro.dat

        first line include NX NY NZ info

        Programmer: Bin Wang (yin.feng@louisiana.edu)
        Creation:   Sep, 2017
        """
        headline='Dimension(NX,NY,NZ): (%s X %s X %s) N=%s'%(self.NX,self.NY,self.NZ,self.N)
        np.savetxt(filehead + "_permx.txt", self.PERMX, delimiter="\n",fmt='%1.4f',header=headline)
        np.savetxt(filehead + "_permy.txt", self.PERMY, delimiter="\n",fmt='%1.4f',header=headline)
        np.savetxt(filehead + "_permz.txt", self.PERMZ, delimiter="\n",fmt='%1.4f',header=headline)
        np.savetxt(filehead + "_poro.txt",  self.PORO, delimiter="\n",fmt='%1.4f',header=headline)
        print(headline)
        print('%s successfully genetrated!' % (filehead + "_permx.txt"))
        print('%s successfully genetrated!' % (filehead + "_permy.txt"))
        print('%s successfully genetrated!' % (filehead + "_permz.txt"))
        print('%s successfully genetrated!' % (filehead + "_poro.txt"))

    def write_VTU(self,filename,mode=0):
        """Convert GRDECL format to unstructured VTU for visulization using Paraview

        Arguments
        ---------
        mode    --  0-Perm_Poro  1-transmisability 2-fault condition

        Programmer: Bin Wang (yin.feng@louisiana.edu)
        Creation:   Sep, 2017
        """
        
        debug=0
        fname=filename

        NX,NY,NZ=self.NX,self.NY,self.NZ
        print('\n--------Output----------\nDimension(NX,NY,NZ): (%s X %s X %s)'%(NX,NY,NZ))
        
        #Convert Grid paramters size(NX,NY,NZ) to nodes parameters size(2NX,2NY,2NZ)
        coordX=DX2CoordX(self.DX,NX,NY,NZ)
        coordY=DY2CoordY(self.DY,NX,NY,NZ)
        coordZ=TOPSDZ2CoordZ(self.TOPS,self.DZ,NX,NY,NZ)
        
        ######Grid Geometry Information######
        points=[]
        for k in range(2*NZ):
            for j in range(2*NY):
                for i in range(2*NX):
                    x,y,z=coordX[i][j][k],coordY[i][j][k],coordZ[i][j][k]
                    points.append([x,y,z])
        
        #print(coordX)
        Cell_hex=[]
        for k in range(NZ):  
            for j in range(NY):
                for i in range(NX):
                    idx_GB=[getIJK(2*i,2*j,2*k,2*NX,2*NY,2*NZ),      #index-0
                            getIJK(2*i+1,2*j,2*k,2*NX,2*NY,2*NZ),    #index-1
                            getIJK(2*i+1,2*j+1,2*k,2*NX,2*NY,2*NZ),  #index-2
                            getIJK(2*i,2*j+1,2*k,2*NX,2*NY,2*NZ),    #index-3
                            getIJK(2*i,2*j,2*k+1,2*NX,2*NY,2*NZ),    #index-4
                            getIJK(2*i+1,2*j,2*k+1,2*NX,2*NY,2*NZ),  #index-5
                            getIJK(2*i+1,2*j+1,2*k+1,2*NX,2*NY,2*NZ),#index-6
                            getIJK(2*i,2*j+1,2*k+1,2*NX,2*NY,2*NZ)]  #index-7
                    #debug
                    if(debug):
                        pts_000=points[getIJK(2*i,2*j,2*k,2*NX,2*NY,2*NZ)]
                        pts_100=points[getIJK(2*i+1,2*j,2*k,2*NX,2*NY,2*NZ)]
                        pts_110=points[getIJK(2*i+1,2*j+1,2*k,2*NX,2*NY,2*NZ)]
                        pts_010=points[getIJK(2*i,2*j+1,2*k,2*NX,2*NY,2*NZ)]
                        pts_001=points[getIJK(2*i,2*j,2*k+1,2*NX,2*NY,2*NZ)]
                        pts_101=points[getIJK(2*i+1,2*j,2*k+1,2*NX,2*NY,2*NZ)]
                        pts_111=points[getIJK(2*i+1,2*j+1,2*k+1,2*NX,2*NY,2*NZ)]
                        pts_011=points[getIJK(2*i,2*j+1,2*k+1,2*NX,2*NY,2*NZ)]
                    
                        #print(idx_GB)
                        #print(pts_000,pts_100, pts_110,pts_010)
                        #print(pts_001,pts_101, pts_111,pts_011)
                    Cell_hex.append(idx_GB)
        
        ######Field Data Information######
        if(mode==0):
            celldata=CellData(Scalars(self.PORO,name='phi'),
                        Scalars(self.PERMX,name='kx'),
                        Scalars(self.PERMY,name='ky'),
                        Scalars(self.PERMZ,name='kz'))
        if(mode==1):
            celldata=CellData(Scalars(self.Trans_i0,name='Tx0'),
                        Scalars(self.Trans_i1,name='Tx1'),
                        Scalars(self.Trans_j0,name='Ty0'),
                        Scalars(self.Trans_j1,name='Ty1'),
                        Scalars(self.Trans_k0,name='Tz0'),
                        Scalars(self.Trans_k1,name='Tz1'))
        if(mode==2):
            celldata=CellData(Scalars(self.fault_i0,name='Fault_x0'),
                        Scalars(self.fault_i1,name='Fault_x1'),
                        Scalars(self.fault_j0,name='Fault_y0'),
                        Scalars(self.fault_j1,name='Fault_y1'))
        
        vtk = VtkData(UnstructuredGrid(points,hexahedron=Cell_hex),celldata,'Unstructured Grid File')
        vtk.tofile(fname)
        
        print('Number of nodes:',len(points))
        print('Number of GB:',len(Cell_hex))
        print('%s.vtk successfully genetrated, pelase use Paraview to open it!'%(fname))


    def write_NPSL(self,filehead):
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
        np.savetxt(filehead + "_permx.txt", self.PERMX, delimiter="\n",fmt='%1.4f')
        np.savetxt(filehead + "_permy.txt", self.PERMY, delimiter="\n",fmt='%1.4f')
        np.savetxt(filehead + "_permz.txt", self.PERMZ, delimiter="\n",fmt='%1.4f')
        np.savetxt(filehead + "_poro.txt",  self.PORO, delimiter="\n",fmt='%1.4f')
        print('\n--------Output----------\nDimension(NX,NY,NZ): (%s X %s X %s)'%(self.NX,self.NY,self.NZ))
        print('Number of GB:', self.N)
        print('%s successfully genetrated, pelase use NPSL to load it!' % (filehead + "_permx.txt"))
        print('%s successfully genetrated, pelase use NPSL to load it!' % (filehead + "_permy.txt"))
        print('%s successfully genetrated, pelase use NPSL to load it!' % (filehead + "_permz.txt"))
        print('%s successfully genetrated, pelase use NPSL to load it!' % (filehead + "_poro.txt"))






#############################################
#
#  Auxilary function
#
#############################################

######[read_GRDECL]######
def is_number(s):
    #Determine a string is a number or not
    #Used in [read_GRDECL] [getBlkdata]
    try:
        float(s)
        return True
    except ValueError:
        pass
 
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    
    try:
        num, val = s.split('*')
        return 2
    except ValueError:
        pass
 
    return False

def getBlkdata(blockData,N):
    #Determine the block type and extract all numbers
    #Used in [read_GRDECL]

    blk_id=0
    blk_dataset=[]
    #Find the delete the header
    for i,data in enumerate(blockData):
        if(is_number(data)==False):#Header
            blk_id=block_indicator(data)
        elif(is_number(data)==2):#num*Vals 
            num,val=data.split('*')
            for it in range(int(num)):
                blk_dataset.append(val)
        else:#Normal Data
            blk_dataset.append(data)
    
    if(len(blk_dataset)!=N and len(blk_dataset)!=2):
        print('Data%s size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(blk_id,len(blk_dataset),N))
    #print(blockData)
    #print(blk_id,blk_dataset)    
    return blk_id,blk_dataset

def block_indicator(name):
    #Assign a specific id for different Keywords
    #Used in [read_GRDECL]
    block_indicator=0
    if(name=='TOPS'):
        block_indicator=1
    if(name=='DX'):
        block_indicator=2
    if(name=='DY'):
        block_indicator=3
    if(name=='DZ'):
        block_indicator=4
    if(name=='PORO'):
        block_indicator=5
    if(name=='PERMX'):
        block_indicator=6
    if(name=='PERMY'):
        block_indicator=7
    if(name=='PERMZ'):
        block_indicator=8
    return block_indicator
#/

######[write_VTU]######
def getI_J_K(ijk,NX,NY,NZ):
    #Find index [i,j,k] from a flat 3D matrix index [i,j,k]
    i,j,k=0,0,0
    
    i=ijk%NX
    j=((int)(ijk / NX)) % NY
    k=(int)(ijk / (NX*NY))
    
    return i,j,k

def getIJK(i,j,k,NX,NY,NZ):
    #Convert index [i,j,k] to a flat 3D matrix index [ijk]
    return i + (NX)*(j + k*(NY))


def DX2CoordX(DX,NX,NY,NZ):
    #Convert self.DX to coordiantes-X of unstructure gridblocks
    #Used in [write_VTU]
    debug=0
    coordX=np.zeros((2*NX,2*NY,2*NZ))
    for k in range(2*NZ):
        for j in range(2*NY):
            for i in range(2*NX):
                #print(DX[i][j])
                i_GB,j_GB,k_GB=(int)(i/2),(int)(j/2),(int)(k/2)
                if(i==0):
                    coordX[0][j][k]=0
                    if(debug):print('First Node',(0),'-')
                if(i>0 and i%2==0):
                    ijk_GB=getIJK(i_GB-1,j_GB,k_GB,NX,NY,NZ)
                    coordX[i-1][j][k]=DX[ijk_GB] #odd
                    coordX[i][j][k]=DX[ijk_GB] #even
                    if(debug):print('pair',(i-1,i),(i_GB-1,j_GB))
                    if(i_GB>1):
                        coordX[i-1][j][k]+=coordX[i-1-1][j][k]
                        coordX[i][j][k]=coordX[i-1][j][k]
                if(i==2*NX-1):
                    ijk_GB=getIJK(NX-1,j_GB,k_GB,NX,NY,NZ)
                    coordX[2*NX-1][j][k]=DX[ijk_GB]+coordX[2*NX-1-1][j][k]
                    if(debug):print('Last Node',(2*NX-1),(NX-1,j_GB))
    return coordX

def DY2CoordY(DY,NX,NY,NZ):
    #Convert self.DY to coordiantes-Y of unstructure gridblocks
    #Used in [write_VTU]
    debug=0
    coordY=np.zeros((2*NX,2*NY,2*NZ))
    for k in range(2*NZ):
        for j in range(2*NY):
            for i in range(2*NX):
                #print(DX[i][j])
                i_GB,j_GB,k_GB=(int)(i/2),(int)(j/2),(int)(k/2)
                if(j==0):
                    coordY[i][0][k]=0
                    if(debug and i==0 and k==0):print('First Node',(0),'-')
                if(j>0 and j%2==0):
                    ijk_GB=getIJK(i_GB,j_GB-1,k_GB,NX,NY,NZ)
                    coordY[i][j-1][k]=DY[ijk_GB] #odd
                    coordY[i][j][k]=DY[ijk_GB] #even
                    if(debug and i==0 and k==0):print('pair',(j-1,j),(i_GB,j_GB-1))
                    if(j_GB>1):
                        coordY[i][j-1][k]+=coordY[i][j-1-1][k]
                        coordY[i][j][k]=coordY[i][j-1][k]
                if(j==2*NY-1):
                    ijk_GB=getIJK(i_GB,NY-1,k_GB,NX,NY,NZ)
                    coordY[i][2*NY-1][k]=DY[ijk_GB]+coordY[i][2*NY-1-1][k]
                    if(debug and i==0 and k==0):print('Last Node',(2*NY-1),(i_GB,NY-1))
    return coordY

def TOPSDZ2CoordZ(TOPS,DZ,NX,NY,NZ):
    #Convert self.DZ and Self.TOPS to coordiantes-z of unstructure gridblocks
    #Used in [write_VTU]
    debug=0
    coordZ=np.zeros((2*NX,2*NY,2*NZ))
    for k in range(2*NZ):
        for j in range(2*NY):
            for i in range(2*NX):
                #print(DX[i][j])
                i_GB,j_GB,k_GB=(int)(i/2),(int)(j/2),(int)(k/2)
                if(k==0):
                    coordZ[i][j][0]=0
                    if(debug and i==0 and k==0):print('First Node',(0),'-')
                if(k>0 and k%2==0):
                    ijk_GB=getIJK(i_GB,j_GB,k_GB-1,NX,NY,NZ)
                    coordZ[i][j][k-1]=DZ[ijk_GB] #odd
                    coordZ[i][j][k]=DZ[ijk_GB] #even
                    if(debug and i==0 and k==0):print('pair',(k-1,k),(i_GB,k_GB-1))
                    if(k_GB>1):
                        coordZ[i][j][k-1]+=coordZ[i][j][k-1-1]
                        coordZ[i][j][k]=coordZ[i][j][k-1]
                if(k==2*NZ-1):
                    ijk_GB=getIJK(i_GB,NY-1,k_GB,NX,NY,NZ)
                    coordZ[i][j][2*NZ-1]=DZ[ijk_GB]+coordZ[i][j][2*NZ-1-1]
                    if(debug and i==0 and k==0):print('Last Node',(2*NY-1),(i_GB,NY-1))
                
    return coordZ
#/

######[calc_Trans]######
def check_fault(ijk1,ijk2,TOPS,DZ):
    # Check the fault situation between two neighboring blocks (ijk1, ijk2)
    # Gridblock1    TOPS1------TOPS1+DZ1
    # Gridblock2        TOPS2------TOPS2+DZ2

    TOPS1,DZ1=TOPS[ijk1],DZ[ijk1]
    TOPS2,DZ2=TOPS[ijk2],DZ[ijk2]

    Length_overlap=overlap(TOPS1,TOPS1+DZ1,TOPS2,TOPS2+DZ2)

    if(Length_overlap==DZ1 and Length_overlap==DZ2): #Fully connected fault, not fault any more
        return 1
    elif(Length_overlap==0):#Sealing fault
        return 0
    else:# Partially connected
        return 0.5

def overlap(min1, max1, min2, max2):
    #Math: overlap of 1D segments
    #Reference: https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments?rq=1
    return max(0, min(max1, max2) - max(min1, min2))


#/
