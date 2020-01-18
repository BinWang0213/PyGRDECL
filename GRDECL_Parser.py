#########################################################################
#       (C) 2017-2018 Department of Petroleum Engineering,              # 
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

import os
import numpy as np

SupportKeyWords=[
    'SPECGRID', #Dimenion of the corner point grid
    'DIMENS',   #Define the dimension of the cartesian grid
    'TOPS','DX','DY','DZ',
    'COORD','ZCORN',
    'PORO',
    'PERMX' , 'PERMXY', 'PERMXZ', 
    'PERMYX', 'PERMY' , 'PERMYZ', 
    'PERMZX', 'PERMZY', 'PERMZ',
    'SW_NPSL'
]

KeyWordsDatatypes=[#Corrsponding data types
    int,
    float,
    int,int,int,int,
    float,float,
    float,
    float,float,float,
    float,float,float,
    float,float,float
]


class GRDECL_Parser:
    def __init__(self,filename='',nx=0,ny=0,nz=0):
        """Eclipse Input file(GRDECL) Parser 
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
        self.GRID_type='NaN'

        #Cartesian gridblock data KeyWords
        self.TOPS=[]
        self.DX=[]
        self.DY=[]
        self.DZ=[]

        #Corner point gridblock data KeyWrods (not support now)
        self.COORD=[]
        self.ZCORN=[] #http://maoxp9.blog.163.com/blog/static/122653420093894133671/

        #Petrophysics data Keywords
        self.SpatialDatas={}
        self.SkipedKeywords=0

        #Read GRDECL file when initializing the class
        if(len(filename)>0):
            self.read_GRDECL()

        #Derived variabls 
        self.CELL_FAULT=[]
    
    ######[read_GRDECL]######
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

        print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(self.fname))

        #Read whole file into list
        f=open(self.fname)
        contents=f.read()
        contents=RemoveCommentLines(contents,commenter='--')
        contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
        contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
        NumKeywords=len(contents_in_block)

        GoodFlag=0
        for i,block in enumerate(contents_in_block):#Keyword, Block-wise
            #Clean the data where no spliter \ provided
            block=scanKeyword(block)

            blockData_raw=block.strip().split()
            Keyword=''
            DataArray=[]
            if(len(blockData_raw)>1):
                if(blockData_raw[0]=='ECHO'): #This keyword may next to real keyword
                    Keyword,DataArray=blockData_raw[1],blockData_raw[2:]                    
                else:
                    Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

            #Read Grid Dimension [SPECGRID] or [DIMENS] 
            if(Keyword=='DIMENS'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.GRID_type='Cartesian'
                self.NX,self.NY,self.NZ=DataArray[0],DataArray[1],DataArray[2]
                self.N=self.NX*self.NY*self.NZ
                print("     Grid Type=%s Grid" %(self.GRID_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
                print("     NumOfGrids=%s"%(self.N))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keyword %d [%s] " %(i+1,Keyword),end='')
                GoodFlag=1
                continue
            elif(Keyword=='SPECGRID'):
                DataArray=np.array(DataArray[:3],dtype=int)
                self.GRID_type='CornerPoint'
                self.NX,self.NY,self.NZ=DataArray[0],DataArray[1],DataArray[2]
                self.N=self.NX*self.NY*self.NZ
                print("     Grid Type=%s" %(self.GRID_type))
                print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
                print("     NumOfGrids=%s"%(self.N))
                print('     NumOfKeywords=%s'%(NumKeywords))
                print("     Reading Keywords [%s] " %(Keyword),end='')
                GoodFlag=1
                continue
            
            if(self.GRID_type=='NaN'):#Skip unnecessary keywords
                continue

            if(Keyword in SupportKeyWords): #We need parse the special format in 
                if(len(DataArray)==1 and '.' in DataArray[0]):
                    folder_name=os.path.dirname(self.fname)
                    DataArray=self.read_IncludeFile(os.path.join(folder_name,DataArray[0]),self.N)
                #print(Keyword,DataArray)
                DataArray=parseDataArray(DataArray)
            

            #Read Grid spatial information, x,y,z ordering
            #Corner point cell
            if(Keyword=='COORD'):# Pillar coords
                assert len(DataArray)==6*(self.NX+1)*(self.NY+1),'[Error] Incompatible COORD data size!'
                self.COORD=np.array(DataArray,dtype=float)       
            elif(Keyword=='ZCORN'):# Depth coords
                assert len(DataArray)==8*self.N, '[Error] Incompatible ZCORN data size!'
                self.ZCORN=np.array(DataArray,dtype=float)
            
            #Cartesian cell
            elif(Keyword=='DX'):# Grid size in X dir
                assert len(DataArray)==self.N, '[Error] Incompatible DX data size!'
                self.DX=np.array(DataArray,dtype=float)
            elif(Keyword=='DY'):# Grid size in Y dir
                assert len(DataArray)==self.N, '[Error] Incompatible DY data size!'
                self.DY=np.array(DataArray,dtype=float)
            elif(Keyword=='DZ'):# Grid size in Z dir
                assert len(DataArray)==self.N, '[Error] Incompatible DZ data size!'
                self.DZ=np.array(DataArray,dtype=float)
            elif(Keyword=='TOPS'):# TOP position
                assert len(DataArray)==self.N, '[Error] Incompatible TOPS data size!'
                self.TOPS=np.array(DataArray,dtype=float)

            #Read Grid Properties information
            else:
                self.LoadVar(Keyword,DataArray,DataSize=self.N)

        f.close()
        assert GoodFlag==1,'Can not find grid dimension info, [SPECGRID] or [DIMENS]!'
        print('.....Done!')


        #Genetrate TOPS for cartesian grid if TOPS if not given
        if(self.GRID_type=='Cartesian' and len(self.TOPS)==0):
            self.TOPS=np.zeros(self.N)
            for k in range(self.NZ-1):
                for j in range(self.NY):
                    for i in range(self.NX):
                        ijk=getIJK(i,j,k,self.NX,self.NY,self.NZ)
                        ijk_next=getIJK(i,j,k+1,self.NX,self.NY,self.NZ)
                        self.TOPS[ijk_next] = self.TOPS[ijk] + self.DZ[ijk]
    

    def buildCartGrid(self,physDim=[100.0,100.0,10.0],gridDims=[10,10,1]):
        """Build simple cartesian grid 
        
        Arguments
        ---------
        physDim  -- physical dimensions of system
        gridDims -- grid dimension of system

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Feb. 2019
        """
        self.NX,self.NY,self.NZ=gridDims
        self.N=self.NX*self.NY*self.NZ
        self.GRID_type='Cartesian'


        print("     Grid Type=%s Grid" %(self.GRID_type))
        print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
        print("     NumOfGrids=%s"%(self.N))

        #Assign value to cart grid
        self.DX=np.ones(self.N)*physDim[0]/self.NX
        self.DY=np.ones(self.N)*physDim[1]/self.NY
        self.DZ=np.ones(self.N)*physDim[2]/self.NZ
        self.TOPS=np.zeros(self.N)
        for k in range(self.NZ-1):
                for j in range(self.NY):
                    for i in range(self.NX):
                        ijk=getIJK(i,j,k,self.NX,self.NY,self.NZ)
                        ijk_next=getIJK(i,j,k+1,self.NX,self.NY,self.NZ)
                        self.TOPS[ijk_next] = self.TOPS[ijk] + self.DZ[ijk]
        
        #Build up basic spatial propertis
        self.SpatialDatas["PERMX"]=np.ones(self.N)*10.0
        self.SpatialDatas["PERMY"]=np.ones(self.N)*10.0
        self.SpatialDatas["PERMZ"]=np.ones(self.N)*10.0
        self.SpatialDatas["PORO"]=np.ones(self.N)*0.3


    def LoadVar(self,Keyword,DataArray,DataSize):
        """Load varables into class
        example:
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(Keyword in SupportKeyWords):#KeyWords Check
            assert len(DataArray)==DataSize,'\n     [Error-%s] Incompatible data size! %d-%d' %(Keyword,len(DataArray),DataSize)
            KeywordID=SupportKeyWords.index(Keyword)
            print('     [%s] '%(Keyword),end='')
            self.SpatialDatas[Keyword]=np.array(DataArray,dtype=KeyWordsDatatypes[KeywordID])
        else:
            if(self.SkipedKeywords==0):print()
            print('     [Warnning] Unsupport keywords[%s]' % (Keyword))
            self.SkipedKeywords+=1


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

                                PERMX_new[ijk_new]=self.SpatialDatas["PERMX"][ijk]
                                PERMY_new[ijk_new]=self.SpatialDatas["PERMY"][ijk]
                                PERMZ_new[ijk_new]=self.SpatialDatas["PERMZ"][ijk]
                                PORO_new[ijk_new]=self.SpatialDatas["PORO"][ijk]
                                ijk_new=ijk_new+1
        
        NewGrid=GRDECL_Parser('',NX_new,NY_new,NZ_new)
        NewGrid.GRID_type='Cartesian' #Currently only support CartGrid
        NewGrid.DX=DX_new
        NewGrid.DY=DY_new
        NewGrid.DZ=DZ_new
        NewGrid.TOPS=self.TOPS
        NewGrid.SpatialDatas["PERMX"]=PERMX_new
        NewGrid.SpatialDatas["PERMY"]=PERMY_new
        NewGrid.SpatialDatas["PERMZ"]=PERMZ_new
        NewGrid.SpatialDatas["PORO"]=PORO_new
        NewGrid.print_info()
        return NewGrid


    def print_info(self):
        print("     Grid Type=%s Grid" %(self.GRID_type))
        print("     Grid Dimension(NX,NY,NZ): (%s x %s x %s)"%(self.NX,self.NY,self.NZ))
        print("     NumOfGrids=%s"%(self.N))

    ######[DataInterperator]######
    def getPillar(self,Pid):
        """Get a pillar line from COORD
        Pillar is the vertical cell edge line (top point-bottm point)
        
        IndexMap of COORD
        [Row1] xtop ytop ztop xbottom ybottom zbottom
        [Row2] xtop ytop ztop xbottom ybottom zbottom
        ....
        Row follows an order of X->Y->Z

        Arguments
        ---------
        Pid -- Pillar index in [COORD]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        id_top=[6*Pid+0,6*Pid+1,6*Pid+2]
        id_bottom=[6*Pid+3,6*Pid+4,6*Pid+5]
        TopPoint=np.array([self.COORD[i] for i in id_top])
        BottomPoint=np.array([self.COORD[i] for i in id_bottom])
        return [TopPoint,BottomPoint]

    def getCellPillars(self,i,j):
        """Obtain the four pillars (p0,p1,p2,p3) of a corner point cell
        The index of pillar
        
        3x3x1 system (2D X-Y plane)
        12--- 13  --- 14  ---15
        |      |       |      |  <- Cell 6,7,8
        8 ---  9  --- 10  ---11
        |      |       |      |  <- Cell 3,4,5
        4 ---  5  ---  6  --- 7
        |      |       |      |  <- Cell 0,1,2
        0 ---  1 ---   2 ---  3
        
        
        The pillars index for a grid follows below ordering (XY Plane)
        pil2   pil3
        *------*
        |      |
        |      |
        *------*
        pil0   pil1

        0   12  3
        1. neighboring cell share one common edge index
        2. ZCORN follows the same order for a cell
        3. Try a 3x3x1 grid system using mrst

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        nx,ny=self.NX+1,self.NY+1
        pil0_id,pil1_id=getIJK(i,j,0,nx,ny,0),getIJK(i+1,j,0,nx,ny,0)
        pil2_id,pil3_id=getIJK(i,j+1,0,nx,ny,0),getIJK(i+1,j+1,0,nx,ny,0)

        return [self.getPillar(pil0_id),self.getPillar(pil1_id),self.getPillar(pil2_id),self.getPillar(pil3_id)]     

    def getCellZ(self,i,j,k):
        """Get the Z coords for a cell
        
        Follow getCornerPointCellIdx convention:
        Z, [0,1,2,3,4,5,6,7]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        CellIds=self.getCornerPointCellIdx(i,j,k)
        return [self.ZCORN[i] for i in CellIds]

    def getCellFaceZ(self,i,j,k,Face='X-,X+,Y-,Y+'):
        """Get the Z coords for a cell
        
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   

        Follow getCornerPointCellIdx convention:
        X-, [0,2,4,6]
        X+, [1,3,5,7]
        Y-, [0,1,4,5]
        Y+, [2,3,6,7]

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        CellIds=self.getCornerPointCellIdx(i,j,k)
        if(Face=="X-"): FaceIds=[CellIds[0],CellIds[2],CellIds[4],CellIds[6]]
        if(Face=="X+"): FaceIds=[CellIds[1],CellIds[3],CellIds[5],CellIds[7]]
        if(Face=="Y-"): FaceIds=[CellIds[0],CellIds[1],CellIds[4],CellIds[5]]
        if(Face=="Y+"): FaceIds=[CellIds[2],CellIds[3],CellIds[6],CellIds[7]]
        
        return [self.ZCORN[i] for i in FaceIds]

    def getCellCoords(self,i,j,k):
        """Get XYZ coords for eight node of a cell

        6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1   

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        XYZ=[]
        Pillars=self.getCellPillars(i,j)
        Zs=self.getCellZ(i,j,k)
        for pi in range(8): # Loop 8 point for each cell
            Pillar_i=pi%4
            XYZ.append(self.interpPtsOnPillar(Zs[pi],Pillars[Pillar_i]))
        return XYZ
    
    def getCell1Coord(self,i,j,k,nodeID):
        """Get XYZ coords for one node of 8 nodes for a cell

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        Pillars=self.getCellPillars(i,j)
        Zs=self.getCellZ(i,j,k)
        Pillar_i=nodeID%4
        return self.interpPtsOnPillar(Zs[nodeID],Pillars[Pillar_i])


    def getCornerPointCellIdx(self,i,j,k):
        """Obtain the eight coords index for a cell

        3x3x1 system (2D X-Y plane)
        30---31,32---33,34---35
        |      |       |      |  <- Cell 6,7,8
        24---25,26---27,28---29
        18---19,20---21,22---23
        |      |       |      |  <- Cell 3,4,5
        12---13,14---15,16---17
        6 --- 7,8 --- 9,10---11
        |      |       |      |  <- Cell 0,1,2
        0 --- 1,2 --- 3,4 --- 5

        Node order convention for a 3D cell
         6----7
         -   -   <-Bottom Face
        4----5
          2----3
         -    -  <-Top Face
        0----1

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        nx,ny,nz=2*self.NX,2*self.NY,2*self.NZ

        p1_id,p2_id=getIJK(2*i,2*j,2*k,nx,ny,nz),getIJK(2*i+1,2*j,2*k,nx,ny,nz)
        p3_id,p4_id=getIJK(2*i,2*j+1,2*k,nx,ny,nz),getIJK(2*i+1,2*j+1,2*k,nx,ny,nz)

        p5_id,p6_id=getIJK(2*i,2*j,2*k+1,nx,ny,nz),getIJK(2*i+1,2*j,2*k+1,nx,ny,nz)
        p7_id,p8_id=getIJK(2*i,2*j+1,2*k+1,nx,ny,nz),getIJK(2*i+1,2*j+1,2*k+1,nx,ny,nz)

        #print(p1_id,p2_id,p3_id,p4_id)#Top Layer
        #print(p5_id,p6_id,p7_id,p8_id)#Bottom Layer

        return p1_id,p2_id,p3_id,p4_id,p5_id,p6_id,p7_id,p8_id

    def interpPtsOnPillar(self,z,Pillar):
        """Obtain the eight coords for a cell
           X,Y coords has to be interpolated from Z
        xy1=xy0+k*z

        Pillar=[(x0 y0 z0),(x1 y1 z1)]
        (x,y,z) is somewhere between (x0 y0 z0) and (x1 y1 z1)

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """
        if(abs(Pillar[1][2]-Pillar[0][2])>1e-8):
            k=(z-Pillar[0][2])/(Pillar[1][2]-Pillar[0][2])
        else:#Degenrated cell
            k=0.0

        x=Pillar[0][0]+k*(Pillar[1][0]-Pillar[0][0])
        y=Pillar[0][1]+k*(Pillar[1][1]-Pillar[0][1])
        
        return np.array([x,y,z])
    
    def detectFaceFault(self,Z_ijk,Z_neigh):
        """#* Check Fault type for a face (X-,X+,Y-,Y+)
        TypeID
        -1  --Fully connected
        0  --Sealing
        >0--Partially Fault, fault gap value, e.g. 1.5..

        Z_ijk face node order
        p2     p3
        *------*
        |      |
        |      |
        *------*
        p0     p1

        Arguments
        ---------
        Z_ijk   -- z value of a cell face,e.g. [0.5,0.2,0.6,0.7]
        Z_neigh -- z value of a neighbor cell face

        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        '''
        overlap_p02=overlap(Z_ijk[0],Z_ijk[2],Z_neigh[0],Z_neigh[2])
        overlap_p13=overlap(Z_ijk[1],Z_ijk[3],Z_neigh[1],Z_neigh[3])

        length_p02=Z_ijk[2]-Z_ijk[0]
        length_p13=Z_ijk[3]-Z_ijk[1]

        gap_p02=length_p02-overlap_p02
        gap_p13=length_p13-overlap_p13

        print('Overlap',overlap_p02,overlap_p13)
        print('Gap',gap_p02,gap_p13)
       

        if(abs(gap_p02)+abs(gap_p13)<1e-10): #Fully connected
           
            return -1.0
        elif(abs(overlap_p02)+abs(overlap_p13)<1e-10): #Sealing fault
            return 0.0
        else:#Partially connected
            return 0.5
        '''
        #Simple method
        diffVec=np.array(Z_ijk)-np.array(Z_neigh)
        if(sum(abs(diffVec))>0): 
            return 0.5
        else:
            return -1


    def isBoundaryCell(self,Cell=[0,0,0],Dim='3D'):
        ''' Check the a given cell is boundary cell or not
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        '''
        count=0
        face=[]
        #Boundary Point
        if(Cell[0]==0): 
            count+=1
            face.append('X-')
        if(Cell[0]==self.NX-1): 
            count+=1
            face.append('X+')
        if(Cell[1]==0): 
            count+=1
            face.append('Y-')
        if(Cell[1]==self.NY-1): 
            count+=1
            face.append('Y+')
        
        if(Dim=="3D"):
            if(Cell[2]==0): 
                count+=1
                face.append('Z-')
            if(Cell[2]==self.NZ-1): 
                count+=1
                face.append('Z-')

        return count,face

    def findCellFault(self,Cell=[0,0,0]):
        ''' Check the fault for 4 faces of a cell [X-,X+,Y-,Y+] 2D
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        '''
        i,j,k=Cell
        Faces=['X-','X+','Y-','Y+']
        Fault=[False,False,False,False]
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'X-')
        if(i!=0): 
            Z_neigh=self.getCellFaceZ(i-1,j,k,'X+')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
            #print(Z_ijk,Z_neigh,FaultMarker)
        if(FaultMarker!=-1.0 or i==0):#This is a fault here
            Fault[0]=True
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'X+')
        if(i!=self.NX-1): 
            Z_neigh=self.getCellFaceZ(i+1,j,k,'X-')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
        if(FaultMarker!=-1.0 or i==self.NX-1):#This is a fault here
            Fault[1]=True

        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'Y-')
        if(j!=0):
            Z_neigh=self.getCellFaceZ(i,j-1,k,'Y+')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
        if(FaultMarker!=-1.0 or j==0):#This is a fault here
            Fault[2]=True
        
        FaultMarker=-1
        Z_ijk=self.getCellFaceZ(i,j,k,'Y+')
        if(j!=self.NY-1): 
            Z_neigh=self.getCellFaceZ(i,j+1,k,'Y-')
            FaultMarker=self.detectFaceFault(Z_ijk,Z_neigh)
            #print(Z_ijk,Z_neigh,FaultMarker)
        if(FaultMarker!=-1.0 or j==self.NY-1):#This is a fault here
            Fault[3]=True

        return Fault

#############################################
#
#  Auxilary function
#
#############################################

def parseDataArray(DataArray):
        """Parse special dataArray format in GRDECL 
        example:
            5*3.0=[3.0 3.0 3.0 3.0 3.0]
            1.0 2*3.0 5.0=[1.0 3.0 3.0 5.0]
        
        Author:Bin Wang(binwang.0213@gmail.com)
        Date: Sep. 2018
        """

        data=[]
        error_count=0
        for value in DataArray:
            if(is_number(value)==2):
                num,val=value.split('*')
                for i in range(int(num)): data.append(val)
            elif(is_number(value)==1):
                data.append(value)
            else:
                error_count+=1
        
        if(error_count>0):
            print(DataArray)
        
        assert error_count==0, '[Error] Can not find any numeric value!'
        
        

        return data

def KeyWordReader(fname,varname, datatype=float):
    #Simple reader to read a file with a input keyword name

    f=open(fname)
    print('[Input] Reading ECLIPSE/PETREL file \"%s\" ....'%(fname))
    contents=f.read()
    f.close()

    contents=RemoveCommentLines(contents,commenter='--')
    contents_in_block=contents.strip().split('/') #Sepeart input file by slash /
    contents_in_block = [x for x in contents_in_block if x]#Remove empty block at the end
    NumKeywords=len(contents_in_block)

    for i,block in enumerate(contents_in_block):#Keyword, Block-wise
        #Clean the data where no spliter \ provided
        block=scanKeyword(block)

        blockData_raw=block.strip().split()
        Keyword,DataArray=blockData_raw[0],blockData_raw[1:]

        if(Keyword==varname):
            print("     Reading Keywords [%s] " %(Keyword))
            DataArray=parseDataArray(DataArray)
            return np.array(DataArray,dtype=datatype)

    print('     [Warnning] Can not find keywords on file[%s]',varname)

    return None


def RemoveCommentLines(data,commenter='--'):
    #Remove comment and empty lines
    data_lines=data.strip().split('\n')
    newdata=[]
    for line in data_lines:
        if line.startswith(commenter) or not line.strip():
            # skip comments and blank lines
            continue   
        newdata.append(line)
    return '\n'.join(newdata)

def scanKeyword(data):
    #scan and find the keyword
    #e.g. ['INIT','DX','2500*30.0'] -> ['DX','2500*30.0']
    for key in SupportKeyWords:
        if (key in data) and (data.find(key)!=0):
            return data[data.find(key):-1]
    return data

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
    
    try: #Special format N*val= [val val val ....]
        num, val = s.split('*')
        return 2
    except ValueError:
        pass
 
    return False

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

def overlap(min1, max1, min2, max2):
    #Math: overlap of 1D segments
    #Reference: https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments?rq=1
    return max(0.0, min(max1, max2) - max(min1, min2))
