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
        self.read_GRDECL()

        #Derived variabls 
        self.Trans_i0,self.Trans_i1=np.zeros(self.N),np.zeros(self.N)
        self.Trans_j0,self.Trans_j1=np.zeros(self.N),np.zeros(self.N)
        self.Trans_k0,self.Trans_k1=np.zeros(self.N),np.zeros(self.N)

        self.fault_i0,self.fault_i1=np.zeros(self.N),np.zeros(self.N)
        self.fault_j0,self.fault_j1=np.zeros(self.N),np.zeros(self.N)

    
    def read_GRDECL(self):
        """Read input file(GRDECL) of Reservoir Simulator- Petrel (Eclipse)  
        file format:http://petrofaq.org/wiki/Eclipse_Input_Data
        
        Arguments
        ---------
        NX, NY, NZ -- Grid dimension.
        
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
            blk_id=0
            if(len(blockData_raw)>0):
                block_id,block_dataset=getBlkdata(blockData_raw,self.N)
                if(debug):
                    print('block',i)
                    print(block_id,block_dataset)
            
            block_dataset=np.array(block_dataset,dtype=float)
            #print(block_dataset.shape,block_dataset)
            
            #Assign read dataset to corrsponding Variables(DX,DY...etc)
            if(block_id==1):
                self.TOPS=block_dataset
            if(block_id==2):
                self.DX=block_dataset
            if(block_id==3):
                self.DY=block_dataset
            if(block_id==4):
                self.DZ=block_dataset*1
            if(block_id==5):
                self.PORO=block_dataset
            if(block_id==6):
                self.PERMX=block_dataset
            if(block_id==7):
                self.PERMY=block_dataset
            if(block_id==8):
                self.PERMZ=block_dataset

    def calc_Trans(self):
        """Calculate the transmisability for each gridblock
            Fault has to be detected, which affects the flow cross-section

        Arguments
        ---------
        GB index:  m-previous  p-next
        i-direction    ...|imjk| ijk |ipjk|...
        j-direction    ...|ijmk| ijk |ijpk|...
        k-direction    ...|ijkm| ijk |ijkp|...
        Khalfijk(p,m)	--harmonic average permeability in i,j,k direction
        TransIJKP,M		--Transmisability in i,j,k direction

        Programmer: Yin Feng (yin.feng@louisiana.edu)
        Creation:   May, 2016
        """
        debug=0
        NX,NY,NZ=self.NX,self.NY,self.NZ

        #Adding virtual gridblocks around the boundary
        kx=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        ky=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        kz=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        A_i0,A_i1=np.zeros((NX + 2)*(NY + 2)*(NZ + 2)),np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        A_j0,A_j1=np.zeros((NX + 2)*(NY + 2)*(NZ + 2)),np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        A_k0,A_k1=np.zeros((NX + 2)*(NY + 2)*(NZ + 2)),np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        dx=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        dy=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))
        dz=np.zeros((NX + 2)*(NY + 2)*(NZ + 2))

        for k in range(1,NZ+1): #[1,NZ]
            for j in range(1,NY+1): #[1,NY]
                for i in range(1,NX+1):#[1,NX]
                    ijk_flow=getIJK(i,j,k,NX+2,NY+2,NZ+2)
                    ijk=getIJK(i-1,j-1,k-1,NX,NY,NZ)
                    
                    kx[ijk_flow]=self.PERMX[ijk]
                    ky[ijk_flow]=self.PERMY[ijk]
                    kz[ijk_flow]=self.PERMZ[ijk]
                    dx[ijk_flow]=self.DX[ijk]
                    dy[ijk_flow]=self.DY[ijk]
                    dz[ijk_flow]=self.DZ[ijk]

                    if(self.GRID_type==0):#Cartesian GB
                        A_i0[ijk_flow]=self.DY[ijk]*self.DZ[ijk]
                        A_i1[ijk_flow]=A_i0[ijk_flow]
                        A_j0[ijk_flow]=self.DX[ijk]*self.DZ[ijk]
                        A_j1[ijk_flow]=A_j0[ijk_flow]
                        A_k0[ijk_flow]=self.DX[ijk]*self.DY[ijk]
                        A_k1[ijk_flow]=A_k0[ijk_flow]

        #Calcualte transmisability in six direction(i01,j01,k01)
        for k in range(NZ):  
            for j in range(NY):
                for i in range(NX):
                    #7 point scheme
                    ijk,ijk_flow=getIJK(i,j,k,NX,NY,NZ),getIJK(i+1,j+1,k+1,NX+2,NY+2,NZ+2)
                    ijk_i0,ijk_i0_flow=getIJK(i-1,j,k,NX,NY,NZ),getIJK(i,j+1,k+1,NX+2,NY+2,NZ+2)
                    ijk_i1,ijk_i1_flow=getIJK(i+1,j,k,NX,NY,NZ),getIJK(i+2,j+1,k+1,NX+2,NY+2,NZ+2)
                    ijk_j0,ijk_j0_flow=getIJK(i,j-1,k,NX,NY,NZ),getIJK(i+1,j,k+1,NX+2,NY+2,NZ+2)
                    ijk_j1,ijk_j1_flow=getIJK(i,j+1,k,NX,NY,NZ),getIJK(i+1,j+2,k+1,NX+2,NY+2,NZ+2)
                    ijk_k0,ijk_k0_flow=getIJK(i,j,k-1,NX,NY,NZ),getIJK(i+1,j+1,k,NX+2,NY+2,NZ+2)
                    ijk_k1,ijk_k1_flow=getIJK(i,j,k+1,NX,NY,NZ),getIJK(i+1,j+1,k+2,NX+2,NY+2,NZ+2)

                    #debug
                    if(debug==1 and i==5 and j==3 and k==0):
                        #print(ijk)
                        #print(ijk_flow,ijk_i0_flow,ijk_i1_flow,ijk_j0_flow,ijk_j1_flow,ijk_k0_flow,ijk_k1_flow)
                        print(self.PERMX[ijk],kx[ijk_flow],self.PERMX[ijk_i0],kx[ijk_i0_flow],self.PERMX[ijk_j0],kx[ijk_j0_flow])
                        print(self.PERMY[ijk],ky[ijk_flow],self.PERMY[ijk_i0],ky[ijk_i0_flow],self.PERMY[ijk_j0],ky[ijk_j0_flow])
                        print(A_i0[ijk_flow],kx[ijk_flow],dx[ijk_flow],A_i1[ijk_i0_flow],kx[ijk_i0_flow],dx[ijk_i0_flow],ijk,ijk_i0) #Tx0

                    #T_avg=A1k1*A2k2/(A1K1DX2+A2K2DX1)   A1,K1 always ijk  A2k2 always six neighbors
                    self.fault_i0[ijk],self.Trans_i0[ijk]=self.calc_FaceT(A_i0[ijk_flow],kx[ijk_flow],dx[ijk_flow],
                                                                          A_i1[ijk_i0_flow],kx[ijk_i0_flow],dx[ijk_i0_flow],
                                                                          ijk,ijk_i0)
                    self.fault_i1[ijk],self.Trans_i1[ijk]=self.calc_FaceT(A_i1[ijk_flow],kx[ijk_flow],dx[ijk_flow],
                                                                          A_i0[ijk_i1_flow],kx[ijk_i1_flow],dx[ijk_i1_flow],
                                                                          ijk,ijk_i1)
                    self.fault_j0[ijk],self.Trans_j0[ijk]=self.calc_FaceT(A_j0[ijk_flow],ky[ijk_flow],dy[ijk_flow],
                                                                          A_j1[ijk_j0_flow],ky[ijk_j0_flow],dy[ijk_j0_flow],
                                                                          ijk,ijk_j0)
                    self.fault_j1[ijk],self.Trans_j1[ijk]=self.calc_FaceT(A_j1[ijk_flow],ky[ijk_flow],dy[ijk_flow],
                                                                          A_j0[ijk_j1_flow],ky[ijk_j1_flow],dy[ijk_j1_flow],
                                                                          ijk,ijk_j1)
                    temp,self.Trans_k0[ijk]=self.calc_FaceT(A_k0[ijk_flow],kz[ijk_flow],dz[ijk_flow],
                                                       A_k1[ijk_k0_flow],kz[ijk_k0_flow],dz[ijk_k0_flow])
                    temp,self.Trans_k1[ijk]=self.calc_FaceT(A_j1[ijk_flow],ky[ijk_flow],dz[ijk_flow],
                                                       A_j0[ijk_k1_flow],ky[ijk_k1_flow],dz[ijk_k1_flow])
                    
    def calc_FaceT(self,A1,k1,dx1,A2,k2,dx2,ijk=-1,ijk_neigh=-1):
        """Calculate the transmisability for each face
           T_avg=A1k1*A2k2/(A1K1DX2+A2K2DX1)
           

        Arguments
        ---------
        ijk         -- The index of center GB
        ijk_neigh   -- The index of connected GB
        A1,K1       -- Area, Perm of center GB  
        A2,k2       -- Area, Perm of connected GB
        fault_face  -- fault condition of center GB at interested face, 0-sealing 1-fully connect 0.5-partially connect

        Programmer: Bin Wang (binwang.0213@gmail.com)
        Creation:   Sep, 2017
        """
        Trans_face=0
        fault_face=0
        
        if(A2*k2==0): #Hit boundary
            Trans_face=0
            fault_face=0
        else:
            if(ijk!=-1 and ijk_neigh!=-1):#X and Y direction for fault check
                fault_face=check_fault(ijk,ijk_neigh,self.TOPS,self.DZ)
            Trans_face=2*1.127e-3*A1*k1*A2*k2/(A1*k1*dx2+A2*k2*dx1)
            if(fault_face==0): #Sealing fault in this face 
                Trans_face=0
        
        return fault_face,Trans_face
        

    def write_VTU(self,filename,mode=0):
        """#Convert GRDECL format to unstructured VTU for visulization using Paraview

        Arguments
        ---------
        mode    --  0-Perm_Poro  1-transmisability 2-fault condition

        Programmer: Yin Feng (yin.feng@louisiana.edu)
        Creation:   May, 2016
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
            celldata=CellData(Scalars(self.PORO*0.01,name='phi'),
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
    
    if(len(blk_dataset)!=N):
        print('Data size %s is not equal to defined block dimension (NX*NY*NZ) %s'%(len(blk_dataset),N))
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
    debug=1
    coordZ=np.zeros((2*NX,2*NY,2*NZ))
    for k in range(2*NZ):
        for j in range(2*NY):
            for i in range(2*NX):
                #print(DX[i][j])
                i_GB,j_GB,k_GB=(int)(i/2),(int)(j/2),(int)(k/2)
                ijk_GB=getIJK(i_GB,j_GB,k_GB,NX,NY,NZ)
                
                if(k==0):#First layer is always start from top
                    coordZ[i][j][0]=TOPS[ijk_GB]
                else:
                    coordZ[i][j][k]=TOPS[ijk_GB]+DZ[ijk_GB]
    return coordZ
#/

######[calc_Trans]######
def check_fault(ijk1,ijk2,TOPS,DZ):
    #Check the fault situation between two neighboring blocks (ijk1, ijk2)
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