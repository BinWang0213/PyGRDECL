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
import math
SupportKeyWords=[
    'SPECGRID', #Dimenion of the corner point grid
    'DIMENS',   #Define the dimension of the cartesian grid
    'TOPS','DX','DY','DZ',
    'COORD','ZCORN','ACTNUM',
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
    float,float,float,
    float,
    float,float,float,
    float,float,float,
    float,float,float,
    float
]

class buildCPGGrid_opt:
    def __init__(self,disturbed=True,flat= False,fault_drop = 0.,fault_nx = 2):
        self.disturbed= disturbed
        self.flat= flat
        self.fault_drop = fault_drop
        self.fault_nx = fault_nx

class nodes:
    def __init__(self):
        self.coords=[]
        self.num=0

class cells:
    def __init__(self):
        self.num=0
        self.indexMap=[]

class faces:
    def __init__(self):
        self.nodes    = np.empty((0,1), int)#next appending rows requires 0 in first dim
        self.nodePos  = np.empty((0,1), int)
        self.neighbors= np.empty((0,2), int)
        self.tag=np.empty((0,1), int)
        self.cellTags=np.empty((0,2), int)

class GRDECL_Parser:
    def __init__(self,filename='',nx=0,ny=0,nz=0):
        """Eclipse Input file(GRDECL) Parser 
        Keywords Reference: file format:http://petrofaq.org/wiki/Eclipse_Input_Data

        Arguments
        ---------
        NX, NY, NZ         -- Grid dimension.
        Trans(i01,j01,k01) -- Transmissibility in i,j,k direction
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
        self.ACTNUM=[]

        #MRST Format grid attributes
        self.gridDim=3
        self.cartDims=[]
        self.nodes=nodes()
        self.cells=cells()
        self.faces=faces()

        # Upscaling and block centered fluid simulation fields
        self.fault_nx=0  #localize simlplegrdecl fault
        self.VOL=[]
        self.coarse2fine_ratio=[1,1,1]
        self.Rx=[]
        self.Ry=[]
        self.Rz=[]

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

    def buildCPGGrid(self, physDims=[1.0, 1.0, 0.5], gridDims=[3, 3, 3], opt=None):
        """Build simple corner point grid
        Arguments
        ---------
        physDim   -- physical dimensions of system
        gridDims  -- grid dimension of system
        opt['disturbed'] -- disturb from cartesian grid: skew pillars
        opt['flat'] -- sinusoidal z
        faultDrop -- fault drop at half the CPG in x direction


        Author:Mustapha Zakari(mustapha.zakari@univ-lorraine.fr)
        Date: May. 2021
        Python version based on a translation of MRST simpleGrdecl
        """
        self.NX, self.NY, self.NZ = gridDims
        self.N = self.NX * self.NY * self.NZ
        self.GRID_type = 'CornerPoint'
        if opt is None:   opt = buildCPGGrid_opt()
        # Correct fault_nx to align with coarse grid
        self.fault_nx = (opt.fault_nx // self.coarse2fine_ratio[0]) * self.coarse2fine_ratio[0]

        print("     Creating Grid:")
        print("       Type: %s Grid" % (self.GRID_type))
        print("       Grid Dimensions (NX,NY,NZ): (%s x %s x %s)" % (self.NX, self.NY, self.NZ))
        print("       Number Of Grid Cells: %s" % (self.N))

        # Create x,y,z vectors. Disturb x,z
        x_ = np.linspace(0, physDims[0], self.NX + 1)
        y_ = np.linspace(0, physDims[1], self.NY + 1)
        z_ = np.linspace(0, physDims[2], self.NZ + 1)
        x, y, z = np.meshgrid(x_, y_, z_, indexing='ij')

        if (opt.disturbed):
            # skew pillars
            xmid = x[self.NX // 2, 0, 0];
            zmid = z[0, 0, self.NZ // 2];
            x = x + 3 * 0.025 * (xmid - abs(x - xmid)) * (z - physDims[2]) / physDims[
                2];  # disturb more nead xmid and zmax
            x = physDims[0] * x / x.max()  # rescale
        if (not opt.flat):
            # add sinusoidal variations for z
            def fz(x, y):
                z = np.sin(math.pi * x) + 0.5 * np.sin(1.5 * math.pi * (y + x));
                return z

            for k in range(gridDims[2] + 1):
                xi = x[:, :, k] / physDims[0];
                eta = y[:, :, k] / physDims[1];
                z[:, :, k] = z[:, :, k] + 1.5 * 0.05 * physDims[0] * fz(xi, eta);
            z = np.sort(z, 2);
        z = physDims[2] * (z - z.min()) / (z.max() - z.min())

        # Create pillars
        npillar = (self.NX + 1) * (self.NY + 1)
        ## Each pillar bottom & top xyz COORD (6 coordinates)
        lines = np.zeros((npillar, 6));
        ## Set columns Bottom/Top (z[0]/z[end]) coordinates: xB yB zB xT yT zT
        for i, V in enumerate([x, y, z]):
            lines[:, i + 0] = V[:, :, 0].reshape((npillar), order='F');  # Bottom:iz=0
            lines[:, i + 3] = V[:, :, -1].reshape((npillar), order='F');  # Top:iz=lastindex=-1
        ## flatten lines to create COORD
        self.COORD = lines.reshape((np.prod(lines.size)))

        # Assign z-coordinates
        ## Add repeated indices to z in each direction:
        ## 0,    1,1,       2,2,   ... ,gridDims(d)-1,gridDims(d)-1, gridDims(d)
        ## 0, 2//2,3//2, 4//2,5//2,... ,                           ,2*gridDims(d)//2
        def repInd(d):
            return np.linspace(1, 2 * gridDims[d], 2 * gridDims[d], dtype=int) // 2

        IX, IY, IZ = repInd(0), repInd(1), repInd(2)
        self.z = z
        self.ZCORN = np.zeros((2 * self.NX, 2 * self.NY, 2 * self.NZ), dtype=float)
        for k in range(2 * self.NZ):
            for j in range(2 * self.NY):
                for i in range(2 * self.NX):
                    self.ZCORN[i, j, k] = z[IX[i], IY[j], IZ[k]]

        # Add a fault drop after half the x range: (odd index)
        ## [0, 1],[2, 3],4, ... ,2*(gridDims(d)-1)],[2*gridDims(d)-1, 2*gridDims(d)]
        iFault = self.fault_nx  # self.NX//2
        self.ZCORN[2 * iFault:, :, :] += opt.fault_drop

        # flatten ZCORN
        self.ZCORN = self.ZCORN.reshape((np.prod(self.ZCORN.size)), order='F')
        # Build up basic spatial propertis
        self.ACTNUM = np.ones((self.N))
        self.SpatialDatas["PERMX"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PERMY"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PERMZ"] = np.ones(self.N) * 10.0
        self.SpatialDatas["PORO"] = np.ones(self.N) * 0.3
        created = "***"
        for keyword, data in self.SpatialDatas.items():
            created += keyword + "***"
        print("       Created: Fields: %s" % created)

    def fill_coarse_grid(self, CoarseMod):
        CoarseGrid = CoarseMod.GRDECL_Data
        CoarseGrid.GRID_type = self.GRID_type
        [rx, ry, rz] = self.coarse2fine_ratio
        Nx = self.NX
        Ny = self.NY
        Nz = self.NZ
        ix = self.fault_nx

        # Define ratios
        # For Rx cut at fault index ix
        self.Rx = coarse_sizes_vect(ix, rx)
        CoarseGrid.fault_nx = len(self.Rx)
        self.Rx += coarse_sizes_vect(Nx - ix, rx, inv=True)

        self.Ry = coarse_sizes_vect(Ny, ry)
        self.Rz = coarse_sizes_vect(Nz, rz)

        # Assign Coarse grid dimensions
        CoarseGrid.NX = len(self.Rx)
        CoarseGrid.NY = len(self.Ry)
        CoarseGrid.NZ = len(self.Rz)
        CoarseGrid.N = CoarseGrid.NX * CoarseGrid.NY * CoarseGrid.NZ

        # Create partitioning indices
        P0 = self.create_partition_indices(CoarseGrid)

        # Assign GRDECL cartesian attributes DX,DY,DZ,TOPS,VOL to coarse grid
        if self.GRID_type == 'CornerPoint':
            self.buildDXDYDZTOPS()
        self.fill_VOL()
        self.fill_coarse_grdecl_DXDYDZTOPS(CoarseGrid)

        # Assign coarse COORD and ZCORN attributes
        if self.GRID_type == 'CornerPoint':
            self.fill_coarse_grdecl_COORD_ZCORN(CoarseGrid)

        # Build up basic spatial propertis
        CoarseGrid.SpatialDatas["PERMX"] = np.ones(CoarseGrid.N) * 10.0
        CoarseGrid.SpatialDatas["PERMY"] = np.ones(CoarseGrid.N) * 10.0
        CoarseGrid.SpatialDatas["PERMZ"] = np.ones(CoarseGrid.N) * 10.0
        CoarseGrid.SpatialDatas["PORO"] = np.ones(CoarseGrid.N) * 0.3
        return P0

    def buildCornerPointNodes(self, addZlayers=True):
        # Construct nodal coordinates for Corner Point Grid
        nx, ny, nz = self.NX, self.NY, self.NZ
        ncell = nx * ny * nz
        npillar = (nx + 1) * (ny + 1)

        # Enumerate pillars in grid associate one pillar per node
        # Id from  reshaping of [1 2 3 ... npillar] (column major Fortran ordering 11,21,12,22)
        pillarId = (np.linspace(1, npillar, npillar, dtype=int)).reshape((nx + 1, ny + 1), order='F')
        pillarId -= 1

        # Assign separately 4 pillar Ids per cell 00,10,01,11
        p1 = pillarId[0:nx, 0:ny];
        p2 = pillarId[1:nx + 1, 0:ny];
        p3 = pillarId[0:nx, 1:ny + 1];
        p4 = pillarId[1:nx + 1, 1:ny + 1];
        del pillarId;
        # Assign 4 pillar Ids per cell (with repetitions) and vector reshape it
        lineID = np.zeros((2 * nx, 2 * ny, 2 * nz), dtype=int)
        for k in range(2 * nz):
            lineID[0:2 * nx:2, 0:2 * ny:2, k] = p1;
            lineID[1:2 * nx:2, 0:2 * ny:2, k] = p2;
            lineID[0:2 * nx:2, 1:2 * ny:2, k] = p3;
            lineID[1:2 * nx:2, 1:2 * ny:2, k] = p4;
            lineID2 = lineID.reshape((8 * ncell), order='F')
        del lineID, p1, p2, p3, p4;

        # Use COORD and 4 pillar Ids per cell to get
        # xmymzmxMyMzM pillar coord for their 8 cell nodes
        COORD2 = self.COORD.reshape((npillar, 6))
        node_pillar = np.zeros((8 * ncell, 6), dtype=float)
        for k in range(8 * ncell):
            for j in range(6):
                node_pillar[k, j] = COORD2[lineID2[k], j]
        del COORD2, lineID2;

        # Reconstruct nodal coordinates from pillars and ZCORN
        linFactor = (self.ZCORN[:] - node_pillar[:, 2]) / (node_pillar[:, 5] - node_pillar[:, 2]);
        self.X = node_pillar[:, 0] + linFactor * (node_pillar[:, 3] - node_pillar[:, 0])
        self.Y = node_pillar[:, 1] + linFactor * (node_pillar[:, 4] - node_pillar[:, 1])

        self.X = (self.X).reshape((2 * nx, 2 * ny, 2 * nz), order='F')
        self.Y = (self.Y).reshape((2 * nx, 2 * ny, 2 * nz), order='F')
        self.Z = (self.ZCORN).reshape((2 * nx, 2 * ny, 2 * nz), order='F')

        ## Reverse z if dz<0
        # Expand actnum by 2 in each direction
        # from numpy import tile
        # self.ACTNUM = (self.ACTNUM).reshape((self.NX, self.NY, self.NZ), order='F')
        # a = tile(self.ACTNUM, (2, 2, 2))
        # z=self.Z; z[a==0]=float('NaN')
        # dz=np.diff(z)
        # self.Z=dz;
        if addZlayers:
            # Add top+bottom layers to ensure correct processing of outer bdry at faults
            minz = (self.Z).min();
            maxz = (self.Z).max();

            e = np.zeros((2 * self.NX, 2 * self.NY, 1))
            self.Z = np.concatenate((minz - 2 + e, minz - 1 + e, self.Z, maxz + 1 + e, maxz + 2 + e), axis=2)

            e1 = (self.X[:, :, 0]).reshape((2 * self.NX, 2 * self.NY, 1));
            e2 = (self.X[:, :, -1]).reshape((2 * self.NX, 2 * self.NY, 1));
            self.X = np.concatenate((e1, e1, self.X, e2, e2), axis=2)

            e1 = (self.Y[:, :, 0]).reshape((2 * self.NX, 2 * self.NY, 1));
            e2 = (self.Y[:, :, -1]).reshape((2 * self.NX, 2 * self.NY, 1));
            self.Y = np.concatenate((e1, e1, self.Y, e2, e2), axis=2)

            # Mark active new layers
            actnum = (self.ACTNUM).reshape(self.NX, self.NY, self.NZ);
            e = np.ones((self.NX, self.NY, 1), dtype=bool)
            self.actnum = np.concatenate((e, actnum, e), axis=2)
            self.numAuxiliaryCells = 2 * np.prod(e.shape)

        # Replace nan coordinates in X,Y,Z by inf to avoid nan
        from numpy import isnan
        for V in [self.X, self.Y, self.Z]:
            V[isnan(V)] = float('Inf')

    def buildDXDYDZTOPS(self):
        # Construct DX,DY,DZ,TOPS from CPG X,Y,Z
        self.DX  =np.zeros((self.NX,self.NY,self.NZ))
        self.DY  =np.zeros((self.NX,self.NY,self.NZ))
        self.DZ  =np.zeros((self.NX,self.NY,self.NZ))
        self.TOPS=np.zeros((self.NX,self.NY,self.NZ))
        self.buildCornerPointNodes(addZlayers=False)
        for k in range(self.NZ):
            for j in range(self.NY):
                for i in range(self.NX):
                    I,J,K=2*i,2*j,2*k
                    self.DX[i,j,k]     = self.X[I+1,J,K]-self.X[I,J,K]
                    self.DY[i,j,k]     = self.Y[I,J+1,K]-self.Y[I,J,K]
                    self.DZ[i,j,k]     = self.Z[I,J,K+1]-self.Z[I,J,K]

        # INitialize TOPS with coordinale COORD lines minimum

        for j in range(self.NY):
            for i in range(self.NX):
                iglob=i+j*(self.NX+1)
                self.TOPS[i, j, 0]=self.COORD.reshape(((self.NX+1)*(self.NY+1),6))[iglob,2]

        for k in range(self.NZ-1):
            for j in range(self.NY):
                for i in range(self.NX):
                    self.TOPS[i, j, k+1] = self.TOPS[i,j,k] + self.DZ[i,j,k]
        self.DX=self.DX.reshape((self.N),order='F')
        self.DY=self.DY.reshape((self.N),order='F')
        self.DZ=self.DZ.reshape((self.N),order='F')
        self.TOPS=self.TOPS.reshape((self.N),order='F')

    def createInitialGrid(self):
        # Find unique points
        self.GRID_type='INVALID'
        self.dim=3;
        N=np.prod(self.Z.size);
        # Return the sorted list/array of unique elements of the coord array
        # sorted following i,j,k implies Z must be the first lexicographic axis
        self.nodes.coords=np.zeros((N,3))
        self.nodes.coords[:,0]=self.Z.reshape((N), order='F')
        self.nodes.coords[:,1]=self.Y.reshape((N), order='F')
        self.nodes.coords[:,2]=self.X.reshape((N), order='F')

        # ia: indices of kept unique values in initial array: a,ia=unique(ar0) <=> a=ar0[ia]
        # ic: indices to reconstruct initial array: a,ia,ic=unique(ar0) <=> ar0=a[ic]
        self.nodes.coords, ia,ic =np.unique(self.nodes.coords, return_index=True, return_inverse=True, axis=0)
        self.nodes.coords=np.fliplr( self.nodes.coords)
        self.nodes.num=self.nodes.coords[:,0].size
        self.cartDims= [number // 2 for number in self.X.shape] #nb of cells in each direction

        self.cells.num=np.prod(self.cartDims)
        # Cell indices growing i,j,k
        self.cells.indexMap=np.array(range(self.cells.num))

        # Store indices in P to reconstruct and repeat all point coords
        self.P=ic.reshape(self.X.shape,order='F')
        # Assign each cart cell each unique Id in array B : growing i,j,k
        self.B=np.array(range(self.cells.num)).reshape(self.cartDims,order='F')

    def index(self,i,j,k,sz):
        I,J,K=  np.meshgrid(i,j,k, indexing='ij')
        #find linear index from 3D indices equiv to sub2ind in matlab
        return np.ravel_multi_index([I,J,K], sz, order='F').reshape(-1,order='F')

    def findFaces(self):
        # Find regular faces including their nodes, cell neighbors and direction
        #         PARAMETERS:
        #    G      - Grid struct that faces will be filled into.
        #
        #    P      - 3d-array of point numbers for each cell in the grid, .i.e.,
        #             the grid cell (i,j,k) has corner points numbers
        #                 n1 = P(2*i-1,2*j-1,2*k-1)
        #                 n2 = P(2*i  ,2*j-1,2*k-1)
        #                 n3 = P(2*i-1,2*j  ,2*k-1)
        #                 n4 = P(2*i  ,2*j  ,2*k-1)
        #                 n5 = P(2*i-1,2*j-1,2*k  )
        #                 n6 = P(2*i  ,2*j-1,2*k  )
        #                 n6 = P(2*i-1,2*j  ,2*k  )
        #                 n8 = P(2*i  ,2*j  ,2*k  )
        #
        #    B      - Block numbers in the grid.  These are included for convenience
        #             as they could easily have been computed from (i,j,k).
        #
        #    actnum - Array of 0/1, 0 indicate inactive cell, 1 active cell.
        #
        #    tags   - Face tag (numeric) that should be inserted into G.
        #
        #    opt    - options struct from which opt.verbose is used.
        #
        sz  = self.P.shape
        szb = self.B.shape

        # Find face corners
        di=1;#i step next index
        dj=sz[0]#j step next index
        dk=dj*sz[1]#k step next index
        # first internal faces (West)
        k=self.index(np.array(range(1,sz[0]-1,2)),np.array(range(0,sz[1],2)),np.array(range(0,sz[2]-1,2)),sz)
        # Internal faces
        # converting unique indices k in first global indices from P
        # f = [P(k), P(k+dj), P(k+dj+dk), P(k+dk)];
        #     [c1    c2       c3          c4     ]
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        self.f=np.ones((c1.size,4))
        C=[c1,c2,c3,c4]
        for i in range(4):
            self.f[:,i]=C[i]
        # first internal unique indices of faces (East)
        k=k+di
        self.g=np.ones((c1.size,4))
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        C=[c1,c2,c3,c4]
        for i in range(4):
            self.g[:,i]=C[i]

        # Neighbor cells
        k=self.index(np.array(range(szb[0]-1)),np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.c1 = self.B[np.unravel_index(k, self.B.shape, 'F')]
        k=self.index(np.array(range(1,szb[0])),np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.c2 = self.B[np.unravel_index(k, self.B.shape, 'F')]

        # Keep only perfectly matching cell pairs (c1,c2)
        self.h=np.all(self.f==self.g,1)
        shape=[i//2 for i in sz];shape[0]-=1
        self.h=np.all(self.h.reshape((shape),order='F'),2)
        self.h = np.repeat(self.h[:, :, np.newaxis], sz[2]//2, axis=2)
        self.c1=self.c1[np.where(self.h.reshape((self.h.size),order='F'))]
        self.c2=self.c2[np.where(self.h.reshape((self.h.size),order='F'))]
        self.h=np.where(self.h.reshape((self.h.size), order='F'))
        self.f = self.f[self.h[0][:],:]
        del self.g,self.h

        # Regular Boundary faces
        k=self.index(np.array([0,sz[0]-1]),np.array(range(0,sz[1],2)),np.array(range(0,sz[2]-1,2)),sz)
        # fB = [self.P(k), P(k + dj), P(k + dj + dk), P(k + dk)];
        c1=self.P[np.unravel_index(k, self.P.shape, 'F')]
        c2=self.P[np.unravel_index(k+dj, self.P.shape, 'F')]
        c3=self.P[np.unravel_index(k+dj+dk, self.P.shape, 'F')]
        c4=self.P[np.unravel_index(k+dk, self.P.shape, 'F')]
        C=[c1,c2,c3,c4]
        self.fB=np.ones((c1.size,4))
        for i in range(4):
            self.fB[:,i]=C[i]
        k=self.index([0,szb[0]-1],np.array(range(szb[1])),np.array(range(szb[2])),szb)
        self.cB = self.B[np.unravel_index(k, self.B.shape, 'F')]
        # # tagB = repmat(tags(:), [numel(cB) / 2, 1]);
        self.tagB= np.tile(self.tag, self.cB.size//2)

        # Append boundary faces
        self.f=np.vstack((self.f,self.fB))
        del self.fB
        self.cB[self.tagB==self.tag[0]]=-1
        self.c1=np.hstack((self.c1,self.cB))#
        self.cB = self.B[np.unravel_index(k, self.B.shape, 'F')]
        self.cB[self.tagB==self.tag[1]]=-1
        self.c2=np.hstack((self.c2,self.cB))#

        # Filter out inactive and degenerate cell pairs
        # Remove inactive cells
        self.actnum_lin=np.reshape(self.actnum,-1)
        self.c1[self.c1 != -1] *= (self.actnum_lin[self.c1 != -1]).astype('int32')
        self.c2[self.c2 != -1] *= (self.actnum_lin[self.c2 != -1]).astype('int32')

        # Remove faces with no neighbors and pinched faces
        ind = ((self.c1 != -1) | (self.c2 != -1)) & ~((self.f[:, 0] == self.f[:, 3]) & (self.f[:, 1] == self.f[:, 2]))
        self.f = self.f[ind, :]
        self.c1 = self.c1[ind]
        self.c2 = self.c2[ind]

        # Remove zero-area faces
        ind =(self.f[:,0]==self.f[:,3]) & (self.f[:,1]==self.f[:,2]);
        self.f = self.f[~ind, :]
        self.c1 = self.c1[~ind]
        self.c2 = self.c2[~ind]


        # Remove repeated nodes from pinch
        self.f[self.f[:,0]==self.f[:,3], 0]=float("nan");
        self.f[self.f[:,1]==self.f[:,2], 0]=float("nan");
        self.F=np.reshape(np.hstack((self.f,float("inf")*np.ones((self.f.shape[0],1),dtype=float))),(-1,1))
        self.nF=np.where(self.F==float("inf"))[0]
        self.nF=np.diff(np.hstack(([-1],self.nF)))-1
        self.F=self.F[self.F[:,0]!=float("inf"),:]

        # Write results to grid srtucture
        n=self.c1.size
        print('Found ',n,' new regular faces')
        # self.faces.nodes =np.empty((0,1),dtype=int)
        self.faces.nodes     = np.vstack((self.faces.nodes,self.F))
        self.faces.neighbors = np.vstack((self.faces.neighbors,np.column_stack((self.c1,self.c2))))
        # self.arr1=np.diff(self.faces.nodePos.astype(float),axis=0),
        self.arr=np.vstack((0,self.nF[:,np.newaxis]))
        self.faces.nodePos   = np.cumsum(self.arr);
        self.faces.tag        = np.vstack((self.faces.tag, np.zeros((n, 1))));
        # self.faces.cellTags   = [self.faces.cellTags,   np.tile(self.tag, (self.c1.size, 1))][1];
        self.faces.cellTags   = np.vstack((self.faces.cellTags, np.repeat(self.tag[np.newaxis,:], self.c1.size, axis=0)));

    def findConnections(self,za,zb):
        #           (1)                       (2)
        #            |                         |
        #  za(i+1,1) o-------------------------o za(i+1,2)
        #            |                         |
        #  zb(j+1,1) * * * *                   |
        #            |       * * * *           |            ^
        #            |               * * * *   |            |
        #            |                       * * zb(j+1,2)  |z positive
        #            |                         |            |
        #            |                         |            |
        #            |                         |
        #    za(i,1) o-------------------------o za(i,2)
        #            |                         |
        #    zb(j,1) * * * * * * * *           |
        #            |               * * * * * * zb(j,2)
        #            |                         |
        #
        # Given:
        # vectors of point numbers a and b where a(:,1), b(:,1) refer to pillar (1)
        # and a(:,2), b(:,2) refer to pillar (2).
        # Each cell a_i is assumed to lie between lines a(i,:) and a(i+1,:) and
        # each cell b_j is assumed to lie between the lines b(j,:) and b(j+1,:).
        # Walk along the stack of cells a_i  and find all connections to cell b_j.
        C = np.zeros((0,2),dtype=int);
        j1 = 0;
        j2 = 0;
        for i in range(za.shape[0]-1):
            j = min(j1, j2);  # Largest j where both
                              # zb(j,1) < za(i,1) and
                              # zb(j,2) < za(i,2)
            # While  b(j,:) (bottom of cell b_j) is below a(i+1,:) (top of cella_i)
            while np.any((zb[j,:] < za[i+1,:])):
                # Precise check to avoid adding pinched layers
                if doIntersect(za[i,:], za[i+1,:], zb[j,:], zb[j+1,:]):
                    C = np.vstack((C, np.array([i, j])));

                # Update candidates for next start of j iteration
                if zb[j,0] < za[i + 1, 0]: j1 = j
                if zb[j,1] < za[i + 1, 1]: j2 = j
                j = j + 1;
                # Precise check to avoid adding pinched layers
                # if  doIntersect(za(i,:), za(i+1,:), zb(j,:), zb(j+1,:)),
        return C

    def intersection(self,La,Lb):
        # Find coordinates of intersection between lines.
        # Each row in La and Lb has two point numbers that are start and
        # endpoints of line segments.  Point coordinates are specified in PTS.
        # This function computes [x,y,z] of intersection between each line pair in
        # La,Lb.
        z  = self.nodes.coords[:,2];
        za = np.reshape(z[La],La.shape,'F'); zb = np.reshape(z[Lb],Lb.shape,'F');
        del z

        # Find parameter t
        t  = np.divide( (zb[:,0]-za[:,0])[:,np.newaxis] , np.diff(za)-np.diff(zb) );
        del zb;

        x = self.nodes.coords[:,0];
        xa = np.reshape(x[La],La.shape,'F'); xb = np.reshape(x[Lb],Lb.shape,'F');

        y = self.nodes.coords[:,1];
        ya = np.reshape(y[La],La.shape,'F'); yb = np.reshape(y[Lb],Lb.shape,'F');

        # compute coordinates
        pts=np.zeros((La.shape[0],3))
        pts[:,0][:,np.newaxis]=xa[:,0][:,np.newaxis]+np.multiply(t,np.diff(xa))
        pts[:,1][:,np.newaxis]=ya[:,0][:,np.newaxis]+np.multiply(t,np.diff(ya))
        pts[:,2][:,np.newaxis]=za[:,0][:,np.newaxis]+np.multiply(t,np.diff(za))
        print('f')
        return pts

    def computeFaceGeometry(self,C):
        pa=np.column_stack((self.a[C[:,0],:] , self.a[C[:,0]+1,:]))
        pb=np.column_stack((self.b[C[:,1],:] , self.b[C[:,1]+1,:]))
        n=pa.shape[0]
        if (n<1):
            numnodes=[]
            Corners=[]
        if np.all(np.all (pa==pb)):
            # All faces match exactly --> no faults here.
            numnodes = 4*np.ones((n,1),dtype=int)
            pa_tr=np.transpose(pa)
            Corners = np.reshape(pa_tr,(-1,1),'F')
        else:
            # We only use z-coordinate of corners to determine intersections
            z=self.nodes.coords[:,2]
            az= np.reshape(z[pa],pa.shape,'F')
            bz= np.reshape(z[pb],pb.shape,'F')
            # Find possible points along each pillar: For each pair of cells, these
            # are the min of upper cornes and the max of lower corners along each
            # pillar
            i = np.column_stack((az[:,:2] < bz[:,:2], az[:,2:] > bz[:,2:]))
            I=np.copy(pa)
            I[np.where(i)[:4]] = pb[np.where(i)[:4]];

            # Are there intersections between upper and lower lines?
            # PP = np.vstack((pa[:, :2],pa[:, 2:],pa[:, :2],pa[:, 2:]))
            # QQ = np.vstack((pb[:, :2],pb[:, 2:],pb[:, 2:],pb[:, :2]))
            PP = np.vstack((np.column_stack((pb[:,0],pa[:,1])),pa[:, 2:],pa[:, :2],pa[:, 2:]))
            QQ = np.vstack((np.column_stack((pa[:,0],pb[:,1])),pb[:, 2:],pb[:, 2:],pb[:, :2]))

            i = (z[PP[:, 0]]-z[QQ[:, 0]] )*(z[PP[:, 1]]-z[QQ[:, 1]] )< 0;
            # Qtemp=self.intersection(np.array([[0, 40]]), np.array([[119, 47]]));
            if (np.sum(i*1) > -10):#If nb of intersection >0
                # Compute intersection coordinates
                Q = self.intersection(PP[np.where(i)[0],:], QQ[np.where(i)[0],:]);
            else:
                Q=np.empty((0,3))
            Q, a, b = np.unique(Q, return_index=True, return_inverse=True, axis=0)

            # Store pt nbs (IDs) of intersections in f.
            # Each column of f represent a type of intersection:
            # f(:,1)  -- A12 x B12  = p2
            # f(:,2)  -- A34 x B34  = p6
            # f(:,3)  -- A12 x B34  = p4 or p8
            # f(:,4)  -- A34 x B12  = p4 or p8
            f=np.ones((n,4))
            f[:,:]=np.nan
            # i[1]=True
            size=np.where(i)[0].size
            for i,ind in enumerate(np.where(i)[0]):
                f[np.unravel_index(ind,f.shape,'F')]=b[i]+self.nodes.coords.shape[0]
            del i,PP,QQ
            # Add intersection pts to list
            self.nodes.coords=np.vstack((self.nodes.coords,Q))
            del Q

            # Point numbers for each face are stored in rows of J.
            # NaN is unassigned or deleted point
            J=np.ones((n,8))
            J[:,:]=np.nan

            # Points on pillars that may be part of face (p1, p3, p5, p7)
            # Just remove duplicate points first...
            I=I.astype(float)
            I[I[:,0]==I[:,2],0]=np.nan;
            I[I[:,1]==I[:,3],1]=np.nan;
            J[:,[0,2,4,6]]=I[:,[0,1,3,2]]
            del I
            # Assign bottom-bottom and top-top intersections. In essence, the upper
            # and lower envelopes of each face are now stored in J.
            J[:,1] = f[:,0]; # p2
            J[:,5] = f[:,1]; # p6


            intersect_botA_topB = np.isnan(f[:,2],where=False);     #  A12 x B34
            intersect_topA_botB = np.isnan(f[:,3],where=False);  #  A34 x B12

            intersect_left = az[:,0] > bz[:,2];

            # Case 1
            ind = intersect_botA_topB & intersect_left;
            J[ind, 0] = J[ind, 6]  = np.nan;
            J[ind, 7]      = f[ind,2]; # p8

            # Case 2
            ind = intersect_botA_topB & ~intersect_left;
            J[ind, 2] = J[ind, 4]  = np.nan;
            J[ind, 3]      = f[ind,2]; # p4

            intersect_left = bz[:, 0] > az[:, 2];  # B12 > A34 on pillar 1

            # Case 4
            ind = intersect_topA_botB & intersect_left;
            J[ind, 0] = J[ind, 6]  = np.nan;
            J[ind, 7]      = f[ind,3]; # p8

            # Case 3
            ind = intersect_topA_botB & ~intersect_left;
            J[ind, 2] = J[ind, 4]  = np.nan;
            J[ind, 3]      = f[ind,3]; # p4

            # Remove repeated points arising from pinch:
            Corners  = np.ones((J.shape[0],1)) * float("Inf") ;
            Corners = (np.hstack((J,Corners))).transpose()
            print('f')
        pa=np.column_stack((self.a[C[:,0],:] , self.a[C[:,0]+1,:]))

    def findFaults(self):
        print('findFaults')
        # Move z-coordinate to first index
        self.P=self.P.transpose((2,0, 1))
        self.B=self.B.transpose((2,0, 1))

        szP = self.P.shape
        szB = self.B.shape
        di=szP[0];  # i step next index
        dj=di * szP[1]  # j step next index
        dk = 1  # k step next index
        # first internal faces (West)
        k = self.index(np.array(range(szP[0])), np.array(range(1, szP[1]-1, 2)),np.array(range(0,szP[2],2)),szP)
        # Point numbers for faces on side A of pillar pairs
        self.a=linear_ind_val_array(self.P,np.column_stack((k,k+dj)))
        # Point numbers for faces on side B of pillar pairs
        self.b=linear_ind_val_array(self.P,np.column_stack((k+di,k+di+dj)))

         # Cells associated with each face on side A and ond side B
        k=self.index(np.array(range(szB[0])),np.array(range(szB[1]-1)),np.array(range(szB[2])),szB)
        self.cA = self.B[np.unravel_index(k, self.B.shape, 'F')]
        k=self.index(np.array(range(szB[0])),np.array(range(1,szB[1])),np.array(range(szB[2])),szB)
        self.cB = self.B[np.unravel_index(k, self.B.shape, 'F')]

        # if four point numbers of a match four point numbers of b, i.e., all
        # a(2i-1:2i, :)==b(2j-1:2j,:), then cellsA(i) and cellsB(j) match exactly
        # along a face.
        #
        # Below, we construct a logical vector that picks out all pillar pairs
        # with ONLY MATCHING cells.
        self.sz2  = sz2=[szP[0], szP[1]//2-1, szP[2]//2];
        h=np.all(self.a==self.b,1)
        h=np.all(h.reshape((sz2),order='F'),0)
        h = np.repeat(h[np.newaxis,:, :], szP[0], axis=0)
        hvec=h.reshape((h.size),order='F')
        print("Found %d faulted stacks\n" % (np.sum(~h)//szP[0]));
        # Keep node numbers of faces that DO NOT match exactly
        self.a=self.a[np.where(~hvec)[0],:]
        self.b=self.b[np.where(~hvec)[0],:]

        # Keep stacks of cells with faces that DO NOT match exactly
        self.cA=self.cA[np.where(~hvec[::2])[0]]
        self.cB=self.cB[np.where(~hvec[::2])[0]]

        # Make artificial z-increment to separate each stack completely in
        # the next function call.
        dz   = np.max(self.nodes.coords[:,2])-np.min(self.nodes.coords[:,2])+1;
        auxz = (np.array(range(1,sz2[1]*sz2[2]+1))*dz).reshape([sz2[1],sz2[2]],order='F');
        auxz =np.repeat(auxz[:,:,np.newaxis],szP[0],axis=2)
        dZ   = auxz.transpose( [2,0,1]);
        dZ=dZ[np.unravel_index(np.where(~hvec)[0], h.shape, 'F')]
        dZ=np.repeat(dZ[:,np.newaxis],2,axis=1)

        # filter out inactive cells from c1 and c2 and faces belonging to inactive
        # cells from a and b. Also, remove entries in dZ.
        ind_a= self.actnum[np.unravel_index(np.repeat(self.cA,2), self.actnum.shape, 'F')];
        ind_b= self.actnum[np.unravel_index(np.repeat(self.cB,2), self.actnum.shape, 'F')];
        self.a = self.a[np.where(ind_a)[0],:];
        self.b = self.b[np.where(ind_b)[0],:];
        self.cA=self.cA[np.where(self.actnum[np.unravel_index(self.cA, self.actnum.shape, 'F')]==1)[0]]
        self.cB=self.cB[np.where(self.actnum[np.unravel_index(self.cB, self.actnum.shape, 'F')]==1)[0]]

        dZa = dZ[np.where(ind_a)[0], :];
        dZb = dZ[np.where(ind_b)[0], :];

        # Find z-coordinates + artificial z-increment of each point in a and b.
        z    = self.nodes.coords[:,2];
        za   = z[self.a]+dZa;
        zb   = z[self.b]+dZb;
        del z;

        # Process connectivity across pillar pair, i.e., determine which faces
        # b(2j-1:2j,:) on side B overlap with faces a(2i-1:2i, :) each face on side
        # A. C returns index pair (2i,2j) for each cell-cell match and a (2i, 2j-1)
        # of (2i-1,2j) pair for each cell-void match. These matches translate to
        # cell connectivity and pieces of faces that are outer or internal
        # boundaries.
        C=self.findConnections(za,zb);

        #  Construct face-to-cell map NEW; each row corresponds to a face, each of
        #  the two columns hold cell numbers.  Cell number zero mark outer/inner
        #  boundary face

        ka   = np.zeros((self.a.shape[0]-1,1),dtype=int)-1; ka[::2] = self.cA[:,np.newaxis];
        kb   = np.zeros((self.b.shape[0]-1,1),dtype=int)-1; kb[::2] = self.cB[:,np.newaxis];
        new  =np.column_stack((ka[C[:,0]], kb[C[:,1]]));
        ind  = np.any(new!=-1, axis=1);  # Keep rows with at least one non-zero
        new  = new[np.where(ind)[0], :];

        self.faces.neighbors=np.vstack((self.faces.neighbors,new))
        self.faces.cellTags   = np.vstack( ( self.faces.cellTags,   np.repeat(self.tag[np.newaxis,:], new.shape[0], axis=0) ) );

        # Compute new node coordinates and geometry of the newly found faces
        self.computeFaceGeometry(C[ind,:]);

        print('f')

    def process_pillar_faces(self):
        print('Processing regular faces')
        self.findFaces()
        self.findFaults()

    def processGRDECL(self):
        # Construct nodal coordinates  from pillars and ZCORN
        # -> self.X,self.Y,self.Z
        self.buildCornerPointNodes()
        # Init Grid nodes,cells attributes
        # Store indices in P to reconstruct all repeated point coords from unique points
        # Store each cart Block/Cell unique Id in array B : lexico growing i,j,k
        # -> self.nodes,self.cells,self.P,self.B
        self.createInitialGrid()

        # # Free X,Y,Z spaces
        # self.X,self.Y,self.Z=[],[],[]

        # Process faces with constant i-index
        self.tag=np.array([1,2]) #West East
        self.process_pillar_faces()

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


    # Coarsening methods
    def create_partition_indices(self, Cgrid):
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Generate partition Fine scale Ids in P0
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Matrix blocks of i*ones(nx/Nx,ny/Ny,nz/Nz) (i=1,2,...,Nx*Ny*Nz)
        Fgrid=self;
        print('[Partioning] Partitoning from Fine grid of size: [%d,%d,%d]' \
              % (Fgrid.NX, Fgrid.NY, Fgrid.NZ), \
              'to Coarse grid of size: [%d,%d,%d]' \
              % (Cgrid.NX, Cgrid.NY, Cgrid.NZ))
        P0 = np.ones((Fgrid.NX, Fgrid.NY, Fgrid.NZ), dtype='int64');
        ind_coarse = -1;
        for k in range(Cgrid.NZ):
            for j in range(Cgrid.NY):
                for i in range(Cgrid.NX):
                    ind_coarse += 1
                    i0,i1 = coarse_start_end_indices(self.Rx,i)
                    j0,j1 = coarse_start_end_indices(self.Ry,j)
                    k0,k1 = coarse_start_end_indices(self.Rz,k)
                    P0[i0:i1, j0:j1, k0:k1] = ind_coarse;

        P0 = np.reshape(P0, (Fgrid.NX * Fgrid.NY * Fgrid.NZ), order='F')  # fortran (start loop with 'i'):: for k: for j: for i:
        return P0

    def fill_coarse_grdecl_DXDYDZTOPS(self,CGrid):
        # Assign GRDECL cartesian attributes DX,DY,DZ,TOPS to coarse grid
        FGrid = self;
        CGrid.DX = np.zeros((CGrid.NX, CGrid.NY, CGrid.NZ))
        CGrid.DY = np.zeros((CGrid.NX, CGrid.NY, CGrid.NZ))
        CGrid.DZ = np.zeros((CGrid.NX, CGrid.NY, CGrid.NZ))
        CGrid.TOPS = np.zeros((CGrid.NX, CGrid.NY, CGrid.NZ))
        DX = np.reshape(FGrid.DX, (FGrid.NX, FGrid.NY, FGrid.NZ), order='F')
        DY = np.reshape(FGrid.DY, (FGrid.NX, FGrid.NY, FGrid.NZ), order='F')
        DZ = np.reshape(FGrid.DZ, (FGrid.NX, FGrid.NY, FGrid.NZ), order='F')
        for k in range(CGrid.NZ):
            for j in range(CGrid.NY):
                for i in range(CGrid.NX):
                    i0,i1 = coarse_start_end_indices(self.Rx,i)
                    j0,j1 = coarse_start_end_indices(self.Ry,j)
                    k0,k1 = coarse_start_end_indices(self.Rz,k)
                    CGrid.DX[i, j, k] = (DX[i0:i1, j0, k0]).sum()
                    CGrid.DY[i, j, k] = (DY[i0, j0:j1, k0]).sum()
                    CGrid.DZ[i, j, k] = (DZ[i0,    j0, k0:k1]).sum()
        assert(np.sum(CGrid.DX[:,0,0])==np.sum(FGrid.DX[0:self.NX])),"SUM(Coarse DX)!=SUM(Fine DX)"


        for k in range(CGrid.NZ-1):
            for j in range(CGrid.NY):
                for i in range(CGrid.NX):
                    CGrid.TOPS[i, j, k+1] = CGrid.TOPS[i, j, k] +CGrid.DZ[i, j, k+1]
        CGrid.DX = np.reshape(CGrid.DX, (-1, 1), order='F')
        CGrid.DY = np.reshape(CGrid.DY, (-1, 1), order='F')
        CGrid.DZ = np.reshape(CGrid.DZ, (-1, 1), order='F')
        CGrid.TOPS = np.reshape(CGrid.TOPS, (-1, 1), order='F')
        CGrid.fill_VOL()

    def fill_coarse_grdecl_COORD_ZCORN(self,CGrid):
        ord='F'
        self.ZCORN=self.ZCORN.reshape((2*self.NX,2*self.NY,2*self.NZ),order=ord)

        CGrid.ZCORN=np.zeros((2*CGrid.NX,2*CGrid.NY,2*CGrid.NZ))
        for K in range(CGrid.NZ):
            for J in range(CGrid.NY):
                for I in range(CGrid.NX):
                    # indices of min/max in corresponding min/max fine cells
                    i0,i1= coarse_start_end_indices(self.Rx,I)
                    iloc=[i0,i1-1]
                    j0,j1= coarse_start_end_indices(self.Ry,J)
                    jloc=[j0,j1-1]
                    k0,k1= coarse_start_end_indices(self.Rz,K)
                    kloc=[k0,k1-1]
                    for K_inc in range(2):
                        for J_inc in range(2):
                            for I_inc in range(2):
                                CGrid.ZCORN[2*I+I_inc,2*J+J_inc,2*K+K_inc]=\
                                 self.ZCORN[2*iloc[I_inc]+I_inc ,2*jloc[J_inc]+J_inc ,2*kloc[K_inc]+K_inc]

        self.COORD=self.COORD.reshape(((self.NX+1)*(self.NY+1),6))
        CGrid.COORD=np.zeros(((CGrid.NX+1)*(CGrid.NY+1),6))


        for J in range(CGrid.NY):
            for I in range(CGrid.NX):
                for J_inc in range(2):
                    for I_inc in range(2):
                        iglob_coarse= (I+I_inc) +(J+J_inc) *(CGrid.NX+1)
                        i0, i1 = coarse_start_end_indices(self.Rx, I)
                        iloc = [i0, i1 ]
                        j0, j1 = coarse_start_end_indices(self.Ry, J)
                        jloc = [j0, j1 ]
                        iglob_fine=iloc[I_inc]+jloc[J_inc]*(self.NX+1)
                        CGrid.COORD[iglob_coarse,:]=self.COORD[iglob_fine ,: ]

        self.ZCORN=self.ZCORN.reshape((8*self.N),order=ord)
        self.COORD=self.COORD.reshape(((self.NX+1)*(self.NY+1)*6))

        CGrid.ZCORN=CGrid.ZCORN.reshape((8*CGrid.N),order=ord)
        CGrid.COORD=CGrid.COORD.reshape(((CGrid.NX+1)*(CGrid.NY+1)*6))
        # print('Coarse CPG filled')

    def fill_local_grdecl_COORD_ZCORN(self,LocalGrid, ind,Glob_ind):
        ord='F'
        self.ZCORN=self.ZCORN.reshape((2*self.NX,2*self.NY,2*self.NZ),order=ord)

        coarseNx=len(self.Rx)
        coarseNy=len(self.Ry)
        coarseNz=len(self.Rz)
        I,J,K=getI_J_K(ind, coarseNx, coarseNy, coarseNz)

        [rx, ry, rz] =  [self.Rx[I], self.Ry[J], self.Rz[K]]
        LocalGrid.ZCORN=np.zeros((2*LocalGrid.NX,2*LocalGrid.NY,2*LocalGrid.NZ))
        for K in range(LocalGrid.NZ):
            for J in range(LocalGrid.NY):
                for I in range(LocalGrid.NX):
                    linear_index=I+LocalGrid.NX*(J+LocalGrid.NY*K)
                    ijk=Glob_ind[linear_index]
                    i,j,k=getI_J_K(ijk,self.NX,self.NY,self.NZ)
                    for K_inc in range(2):
                        for J_inc in range(2):
                            for I_inc in range(2):
                                LocalGrid.ZCORN[2*I+I_inc,2*J+J_inc,2*K+K_inc]=\
                                 self.ZCORN[2*i+I_inc,2*j+J_inc,2*k+K_inc]
        del I_inc,J_inc,K_inc
        self.COORD=self.COORD.reshape(((self.NX+1)*(self.NY+1),6))
        LocalGrid.COORD=np.zeros(((LocalGrid.NX+1)*(LocalGrid.NY+1),6))
        pillarIds=[]
        nx=self.NX+1
        ny=self.NY+1
        nz=self.NZ+1
        for J in range(LocalGrid.NY):
            for I in range(LocalGrid.NX):
                linear_index=I+J*(LocalGrid.NX)
                ijk=Glob_ind[linear_index]
                i, j, k = getI_J_K(ijk, self.NX, self.NY, self.NZ)
                pillarIds.append( getIJK(i  ,j  ,0,nx,ny,nz) )
                pillarIds.append( getIJK(i+1,j  ,0,nx,ny,nz) )
                pillarIds.append( getIJK(i  ,j+1,0,nx,ny,nz) )
                pillarIds.append( getIJK(i+1,j+1,0,nx,ny,nz) )
        pillarIds=sorted(list(set(pillarIds)))
        assert(len(pillarIds)==(LocalGrid.NX+1)*(LocalGrid.NY+1)),"PBM creating local mesh pillar numbers!=(nx+1)x(ny+1)"
        ind_local_pillar=0
        for ipillar in pillarIds:
            LocalGrid.COORD[ind_local_pillar,:]=self.COORD[ipillar,:]
            ind_local_pillar+=1

        self.ZCORN=self.ZCORN.reshape((8*self.N),order=ord)
        self.COORD=self.COORD.reshape(((self.NX+1)*(self.NY+1)*6))

        LocalGrid.ZCORN=LocalGrid.ZCORN.reshape((8*LocalGrid.N),order=ord)
        LocalGrid.COORD=LocalGrid.COORD.reshape(((LocalGrid.NX+1)*(LocalGrid.NY+1)*6))
        # print('Local CPG filled')

    def fill_local_grid(self,LocalMod,Partition,ind):
        LocalGrid = LocalMod.GRDECL_Data
        LocalGrid.GRID_type = self.GRID_type

        # Assign Local grid dimensions
        coarseNx=len(self.Rx)
        coarseNy=len(self.Ry)
        coarseNz=len(self.Rz)
        I,J,K=getI_J_K(ind, coarseNx, coarseNy, coarseNz)
        LocalGrid.NX = self.Rx[I]
        LocalGrid.NY = self.Ry[J]
        LocalGrid.NZ = self.Rz[K]
        LocalGrid.N = LocalGrid.NX*LocalGrid.NY*LocalGrid.NZ

        # Get local attributes
        # Global indices
        Glob_ind=np.where(Partition==ind)[0]
        assert(len(Glob_ind) ==LocalGrid.N ),"PBM number of indices not matching partition and local grid"
        # local DX DY DZ TOPS
        LocalGrid.DX=np.array(self.DX)[Glob_ind]
        LocalGrid.DY=np.array(self.DY)[Glob_ind]
        LocalGrid.DZ=np.array(self.DZ)[Glob_ind]
        LocalGrid.TOPS=np.array(self.TOPS)[Glob_ind]

        # Assign coarse COORD and ZCORN attributes
        if self.GRID_type=='CornerPoint':
            self.fill_local_grdecl_COORD_ZCORN(LocalGrid,ind,Glob_ind)

        # Build up basic spatial propertis
        # LocalGrid.CreateCellData(varname="PERMX",  val_array= \
        #     np.array(self.SpatialDatas["PERMX"])[Glob_ind])
        # LocalGrid.CreateCellData(varname="PERM",  val_array= \
        #     np.arrayself.SpatialDatas["PERMY"])[Glob_ind])
        # LocalGrid.SpatialDatas["PERMZ"] = np.ones(LocalGrid.N) * 10.0
        # LocalGrid.SpatialDatas["PORO"] = np.ones(LocalGrid.N) * 0.3
        LocalGrid.SpatialDatas["PERMX"] = np.array(self.SpatialDatas["PERMX"])[Glob_ind]
        LocalGrid.SpatialDatas["PERMY"] = np.array(self.SpatialDatas["PERMY"])[Glob_ind]
        LocalGrid.SpatialDatas["PERMZ"] = np.array(self.SpatialDatas["PERMZ"])[Glob_ind]
        LocalGrid.SpatialDatas["PORO"]  = np.array(self.SpatialDatas["PORO"])[Glob_ind]
        return Glob_ind

#############################################
#
#  Auxiliary function
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



def linear_ind_val_array(a,T):
    ncol=T.shape[1]
    out = np.ones(T.shape, dtype=int)
    for i in range(ncol):
        out[:, i] = a[np.unravel_index(T[:,i], a.shape, 'F')]
    return out


def doIntersect(za1, za2, zb1, zb2):
     #  Does cells given by points (a11,a12,a21,a22) and (b11,b12,b21,b22)
     #  given along pillar (1) and pillar (2) have an nonzero intersection?
     #  We need only the z-ccordinate to check:
     #
     #        |                 |
     # zb1(1) *                 |
     # za1(1) o-*---------------o za1(2)
     #        |. .*             |
     #        | . . *           |
     #        |. . . .*         |
     #        | . . . . *       |
     # za2(1) o-----------*-----o za2(2)
     #        |             *   |
     #        |               * |
     #        |                 * zb1(2)
     #        |                 |
     # zb2(1) * * * * * * * * * * zb2(2)
     #        |                 |
     #       (1)               (2)

    val = overlap (za1[0], za2[0], zb1[0], zb2[0]) or \
          overlap (za1[1], za2[1], zb1[1], zb2[1]) or \
          (za1[0]-zb1[0])*(za1[1]-zb1[1]) < 0     or \
          (za2[0]-zb2[0])*(za2[1]-zb2[1]) < 0
    if np.all(za1-za2==0): val=False
    if np.all(zb1-zb2==0): val=False
    return val

def overlap(min1, max1, min2, max2):
    #Math: overlap of 1D segments
    #Reference: https://stackoverflow.com/questions/16691524/calculating-the-overlap-distance-of-two-1d-line-segments?rq=1
    return max(0.0, min(max1, max2) - max(min1, min2))

# Compute coarse dimensions of hexaedric domain buildCPG
def coarse_sizes_vect(N,r,inv=False):
    assert (N>1),"No coarsening if N<=1"
    # assert(r<N),"No coarsening if r>=N"
    if (N % r >1):
        if inv:
            return [N % r]    + [r]*(N // r)
        else:
            return [r]*(N//r) + [N % r]
    elif(N%r==0):
        return [r] * (N // r)
    if inv:
        return [r+N%r] + [r]*(N//r-1)
    else:
        return [r]*(N//r-1) +[r+N%r]

def coarse_start_end_indices(R,i):
    istart=int(np.sum(R[:i]))
    iend= int(istart+(R[i]))
    return istart,iend